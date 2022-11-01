#include "voronoi_diagrams.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>

namespace fs = std::filesystem;

namespace rs_toolset {
namespace mosaicking {

GDALDatasetUniquePtr VoronoiDiagramsImpl::Run(
    const std::vector<std::string>& rasters_path,
    const std::string& output_path,
    OGRSpatialReference* spatial_ref,
    bool with_refinement,
    int low_overview_trunc,
    int high_overview_trunc,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  std::string string(
      "Running a voronoi diagrams task with {}\n"
      " - Spatial reference name: {}\n"
      " - With refinement: {}\n"
      " - Low overview trunction: {}\n"
      " - High overview trunction: {}\n"
      " - RGB bands' map: ");
  if (rgb_bands_map.empty()) {
    string.append("1,2,3");
  } else if (rgb_bands_map.size() != 3) {
    spdlog::error(
        "The size of \"rgb_bands_map\" {} must be 3", rgb_bands_map.size());
    return nullptr;
  } else {
    for (const auto& idx : rgb_bands_map)
      string.append(std::to_string(idx)).append(",");
    string.pop_back();
  }
  spdlog::info(
      string, output_path, spatial_ref->GetName(), with_refinement,
      low_overview_trunc, high_overview_trunc);
  double geotrans[6];
  for (const auto& path : rasters_path) {
    if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
            path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)); !dataset) {
      spdlog::error("Opening {} failed", path);
      return nullptr;
    } else if (dataset->GetGeoTransform(geotrans) != CE_None ||
        !dataset->GetSpatialRef()) {
      spdlog::error(
          "{} does not have the geotransform or the spatial reference. "
          "Please check whether the raster is DOM", path);
      return nullptr;
    }
  }
  auto driver(utils::GetVectorDriverByPath(output_path));
  if (!driver) {
    spdlog::error("Finding the driver for {} failed", output_path);
    return nullptr;
  }
  GDALDatasetUniquePtr output_dataset(driver->Create(
      output_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
  if (!output_dataset) {
    spdlog::error("Creating {} failed", output_path);
    return nullptr;
  }
  if (!spatial_ref) {
    spdlog::error(
        "No spatial reference is specified for the output composite table");
    return nullptr;
  }
  if (with_refinement && !mosaicking_) {
    spdlog::error("No mosaicking can be used for the seamline refinement");
    return nullptr;
  }
  if (low_overview_trunc < 0 || high_overview_trunc < 0) {
    spdlog::error(
        "The low overview trunction {} and the high overview trunction {} "
        "must be non-negative", low_overview_trunc, high_overview_trunc);
    return nullptr;
  }
  string = "The following raster(s) will be operated:\n";
  for (const auto& path : rasters_path)
    string.append(path).append("\n");
  string.append("{} task(s) in total");
  spdlog::info(string, rasters_path.size());

  auto output_layer(output_dataset->CreateLayer(
      "", spatial_ref, wkbPolygon, nullptr));
  OGRFieldDefn path_field("path", OFTString);
  output_layer->CreateField(&path_field);

  // Create borders and calculate centroids
  spdlog::debug("Creating borders and calculating centroids");
  auto rasters_count(static_cast<int>(rasters_path.size()));
  std::vector<OGRGeometryUniquePtr> borders(rasters_count);
  std::vector<cv::Point2f> centroids(rasters_count);
#pragma omp parallel for schedule(dynamic)
  for (int i(0); i < rasters_count; ++i) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        rasters_path[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    OGRPoint point;
    borders[i] = utils::CreateBorder(dataset.get(), spatial_ref);
    borders[i]->Centroid(&point);
    centroids[i] = {
        static_cast<float>(point.getX()), static_cast<float>(point.getY())};
  }
  spdlog::info("Creating borders and calculating centroids - done");

  // Run the voronoi diagrams algorithm
  spdlog::debug("Running the voronoi diagrams algorithm");
  OGREnvelope enve, union_enve;
  for (const auto& border : borders) {
    border->getEnvelope(&enve);
    union_enve.Merge(enve);
  }
  cv::Subdiv2D subdiv({
      static_cast<int>(floor(union_enve.MinX)),
      static_cast<int>(floor(union_enve.MinY)),
      static_cast<int>(ceil(union_enve.MaxX) - floor(union_enve.MinX)),
      static_cast<int>(ceil(union_enve.MaxY) - floor(union_enve.MinY))});
  subdiv.insert(centroids);
  std::vector<std::vector<cv::Point2f>> polygons_points;
  std::vector<cv::Point2f> polygons_center;
  subdiv.getVoronoiFacetList({}, polygons_points, polygons_center);

  // Create faetures within a rough rectangle
  auto rect(utils::CreateGeometryFromEnve(union_enve));
  for (int i(0); i < rasters_count; ++i) {
    auto linear_ring(new OGRLinearRing);
    for (const auto& point : polygons_points[i])
      linear_ring->addPoint(point.x, point.y);
    linear_ring->addPoint(polygons_points[i][0].x, polygons_points[i][0].y);
    OGRGeometryUniquePtr geometry(new OGRPolygon);
    geometry->toPolygon()->addRingDirectly(linear_ring);
    OGRFeatureUniquePtr feature(OGRFeature::CreateFeature(
        output_layer->GetLayerDefn()));
    feature->SetGeometryDirectly(geometry->Intersection(rect.get()));
    feature->SetField(0, rasters_path[i].c_str());
    output_layer->CreateFeature(feature.get());
  }
  spdlog::info("Running the voronoi diagrams algorithm - done");

  // Create features' neighbor indexes
  std::vector<std::vector<int>> features_neighbor_idxes(rasters_count);
  for (int i(0); i < rasters_count; ++i) {
    OGRFeatureUniquePtr feature1(output_layer->GetFeature(i));
    for (int j(i + 1); j < rasters_count; ++j) {
      if (OGRFeatureUniquePtr feature2(output_layer->GetFeature(j));
          feature1->GetGeometryRef()->Intersect(feature2->GetGeometryRef())) {
        features_neighbor_idxes[i].push_back(j);
        features_neighbor_idxes[j].push_back(i);
      }
    }
  }

  // Redefine initial mosaicking geometries points
  spdlog::debug("Redefining initial mosaicking geometries points");
  std::vector<OGRGeometryUniquePtr> dislocated_geometries;
  for (int i(0); i < rasters_count; ++i) {
    OGRFeatureUniquePtr feature1(output_layer->GetFeature(i));
    OGRGeometryUniquePtr boundary(borders[i]->Boundary());
    for (const auto& j : features_neighbor_idxes[i]) {
      OGRFeatureUniquePtr feature2(output_layer->GetFeature(j));
      OGRGeometryUniquePtr
          geometry2(feature2->StealGeometry()),
          seamline(feature1->GetGeometryRef()->Intersection(geometry2.get()));
      if (seamline->getGeometryType() == wkbMultiLineString)
        utils::JointMultiLineString(seamline);
      if (seamline->Intersect(boundary.get())) {
        OGRGeometryUniquePtr geometry1(feature1->GetGeometryRef()->clone());
        utils::RearrangeNeighborGeometries(
            seamline.get(), geometry1, geometry2);
        seamline.reset(seamline->Difference(boundary.get()));
        if (seamline->getGeometryType() == wkbMultiLineString)
          utils::JointMultiLineString(seamline);
        utils::GraftSeamline(seamline, geometry1, geometry2);
        feature2->SetGeometryDirectly(geometry2.release());
        output_layer->SetFeature(feature2.get());
      }
    }

    // Add dislocated geometries
    switch (OGRGeometryUniquePtr diff_geometry(
            feature1->GetGeometryRef()->Difference(borders[i].get()));
        diff_geometry->getGeometryType()) {
      case wkbPolygon: {
        dislocated_geometries.push_back(std::move(diff_geometry));
        break;
      }
      case wkbMultiPolygon: {
        for (const auto& geometry : diff_geometry->toMultiPolygon())
          dislocated_geometries.emplace_back(geometry->clone());
      }
    }

    feature1->SetGeometryDirectly(feature1->GetGeometryRef()->Intersection(
        borders[i].get()));
    output_layer->SetFeature(feature1.get());
  }
  spdlog::info("Redefining initial mosaicking geometries points - done");

  // Handle the internal dislocated geometries
  spdlog::debug("Handling the internal dislocated geometries");
  for (const auto& geometry : dislocated_geometries) {
    for (int i(0); i < rasters_count; ++i) {
      if (OGRFeatureUniquePtr feature(output_layer->GetFeature(i));
          geometry->Intersect(feature->GetGeometryRef()) &&
          borders[i]->Contains(geometry.get())) {
        if (OGRGeometryUniquePtr new_geometry(geometry->Union(
                feature->GetGeometryRef()));
            new_geometry->getGeometryType() == wkbPolygon) {
          feature->SetGeometryDirectly(new_geometry.release());
          output_layer->SetFeature(feature.get());
          break;
        }
      }
    }
  }
  spdlog::debug("Handling the internal dislocated geometries - done");

  // Create neighbor pairs
  std::vector<std::pair<int, int>> neighbor_pairs;
  for (int i(0); i < rasters_count; ++i) {
    OGRFeatureUniquePtr feature1(output_layer->GetFeature(i));
    for (int j(i + 1); j < rasters_count; ++j) {
      if (OGRFeatureUniquePtr feature2(output_layer->GetFeature(j));
          feature1->GetGeometryRef()->Intersect(feature2->GetGeometryRef())) {
        neighbor_pairs.push_back({i, j});
      }
    }
  }

  // Refine seamlines in the initial mosaicking network if the mosaicking is given
  if (with_refinement)
    RefineSeamlines(
        rasters_path, borders, neighbor_pairs, low_overview_trunc,
        high_overview_trunc, rgb_bands_map, color_balancing, output_layer);

  for (auto& feature : output_layer) {
    feature->SetField(
        0, fs::path(feature->GetFieldAsString(0)).filename().string().c_str());
    output_layer->SetFeature(feature.get());
  }
  spdlog::info("Running a voronoi diagrams task - done");
  return output_dataset;
}

void VoronoiDiagramsImpl::RefineSeamlines(
    const std::vector<std::string>& paths,
    const std::vector<OGRGeometryUniquePtr>& borders,
    const std::vector<std::pair<int, int>>& neighbor_pairs,
    int low_overview_trunc,
    int high_overview_trunc,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    OGRLayer* layer) {
  // Create refined seamlines
  spdlog::info("Creating refined seamlines");
  auto neighbors_count(static_cast<int>(neighbor_pairs.size()));
  std::vector<OGRGeometryUniquePtr> geometries(neighbors_count);
  std::vector<std::string> color_balancing_names;
  if (color_balancing)
    color_balancing_names = color_balancing->ExportAllRastersName();
  std::vector<int> color_balancing_idxes(neighbors_count, -1);
  for (int i(0); i < paths.size(); ++i) {
    OGRFeatureUniquePtr feature(layer->GetFeature(i));
    geometries[i].reset(feature->StealGeometry());
    if (color_balancing) {
      if (auto it(std::find(
              color_balancing_names.begin(), color_balancing_names.end(),
              fs::path(paths.back()).filename().string()));
          it != color_balancing_names.end()) {
        color_balancing_idxes[i] = (static_cast<int>(
            it - color_balancing_names.begin()));
      }
    }
  }
  std::vector<OGRGeometryUniquePtr> seamlines(neighbors_count);
  spdlog::default_logger()->set_level(spdlog::level::err);
  auto parallel_sink(spdlog::stdout_color_mt("Parallel"));
  int count(0);
#pragma omp parallel for schedule(dynamic)
  for (int i(0); i < neighbors_count; ++i) {
    auto idx1(neighbor_pairs[i].first), idx2(neighbor_pairs[i].second);
    seamlines[i] = mosaicking_->RunTaskForPair(
        paths[idx1], paths[idx2], borders[idx1].get(), borders[idx2].get(),
        geometries[idx1].get(), geometries[idx2].get(), low_overview_trunc,
        high_overview_trunc, rgb_bands_map, color_balancing,
        color_balancing_idxes[idx1], color_balancing_idxes[idx2]);
#pragma omp critical
    {
      parallel_sink->info(
          "---------- {}/{} ({} and {}) - done ----------", ++count,
          neighbors_count, fs::path(paths[idx1]).stem().string(),
          fs::path(paths[idx2]).stem().string());
    }
  }
  spdlog::default_logger()->set_level(spdlog::level::info);

  // Update features with refined seamlines
  spdlog::debug("Updating features with refined seamlines");
  for (int i(0); i < neighbors_count; ++i) {
    OGRFeatureUniquePtr
        feature1(layer->GetFeature(neighbor_pairs[i].first)),
        feature2(layer->GetFeature(neighbor_pairs[i].second));
    OGRGeometryUniquePtr
        geometry1(feature1->StealGeometry()),
        geometry2(feature2->StealGeometry()),
        seamline(geometry1->Intersection(geometry2.get()));
    if (seamline->getGeometryType() == wkbMultiLineString)
      utils::JointMultiLineString(seamline);
    utils::RearrangeNeighborGeometries(seamline.get(), geometry1, geometry2);
    utils::GraftSeamline(seamlines[i], geometry1, geometry2);
    feature1->SetGeometryDirectly(geometry1.release());
    feature2->SetGeometryDirectly(geometry2.release());
    layer->SetFeature(feature1.get());
    layer->SetFeature(feature2.get());
  }
  spdlog::info("Updating features with refined seamlines - done");
}

std::shared_ptr<VoronoiDiagrams> VoronoiDiagrams::Create(
    const std::shared_ptr<MosaickingInterface>& mosaicking) {
  return std::make_shared<VoronoiDiagramsImpl>(mosaicking);
}

}  // namespace mosaicking
}  // namespace rs_toolset