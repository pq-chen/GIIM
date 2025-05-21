#include "mosaicking_base.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gdalwarper.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>
#include <rs-toolset/color_balancing.h>

namespace fs = std::filesystem;

namespace {

using namespace rs_toolset;

bool RunTaskForSerialCheck(
    const std::string& path,
    OGRLayer* mosaicking_layer,
    OGRLayer* border_layer,
    OGRGeometry* covered_border,
    std::vector<double>& borders_area,
    int& low_overviews_trunc,
    int& high_overviews_trunc,
    double surrounded_buffer,
    double rejection_ratio,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    GDALDatasetUniquePtr& source_raster_dataset,
    std::vector<int>& new_rgb_bands_map,
    int& color_balancing_idx) {
  if (source_raster_dataset.reset(GDALDataset::Open(
          path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
      !source_raster_dataset) {
    spdlog::error("Opening {} failed", path);
    return false;
  } else if (double geotrans[6];
      source_raster_dataset->GetGeoTransform(geotrans) != CE_None ||
      !source_raster_dataset->GetSpatialRef()) {
    spdlog::error(
        "{} does not have the geotransform or the spatial reference. "
        "Please check whether the raster is DOM", path);
    return false;
  }
  if (!mosaicking_layer || !border_layer || !covered_border) {
    spdlog::error(
        "One of the \"mosaicking_layer\", the \"border_layer\" and the "
        "\"covered_border\" is empty");
    return false;
  }
  if (mosaicking_layer->GetFeatureCount() != border_layer->GetFeatureCount()) {
    spdlog::error(
        "The features' count {} in the \"mosaicking_layer\" must be the same "
        "with the features' count {} in the \"border_layer\"",
        mosaicking_layer->GetFeatureCount(), border_layer->GetFeatureCount());
    return false;
  }
  if (covered_border->getGeometryType() != wkbMultiPolygon) {
    spdlog::error(
        "The \"covered_border\" geometry type {} must be 6(wkbMultiPolygon)",
        int(covered_border->getGeometryType()));
    return false;
  }
  if (borders_area.size() != border_layer->GetFeatureCount()) {
    spdlog::error(
        "The \"borders_area\" size {} must be the same with the features' "
        "count {} in the \"border_layer\"", borders_area.size(),
        border_layer->GetFeatureCount());
    return false;
  }
  if (low_overviews_trunc < 0 || high_overviews_trunc < 0) {
    spdlog::error(
        "The \"low_overviews_trunc\" {} and the \"high_overviews_trunc\" {} "
        "must be non-negative", low_overviews_trunc, high_overviews_trunc);
    return false;
  }
  utils::CreateRasterPyra(source_raster_dataset.get());
  if (auto overviews_count(
          source_raster_dataset->GetRasterBand(1)->GetOverviewCount());
      (low_overviews_trunc + high_overviews_trunc) >= overviews_count) {
    spdlog::warn(
        "The \"low_overviews_trunc\" {} plus the \"high_overviews_trunc\" {} "
        "must be less than the overviews' count {}", low_overviews_trunc,
        high_overviews_trunc, overviews_count);
    high_overviews_trunc = std::min(high_overviews_trunc, overviews_count - 1);
    low_overviews_trunc = overviews_count - 1 - high_overviews_trunc;
    spdlog::info(
        "Reset arguments:\n"
        " - Low overviews trunction: {}\n"
        " - High overviews trunction: {}", low_overviews_trunc,
        high_overviews_trunc);
  }
  if (surrounded_buffer < 0.0) {
    spdlog::error(
        "The \"surrounded_buffer\" {} must be non-negative", surrounded_buffer);
    return false;
  }
  if (rejection_ratio < 0.0 || rejection_ratio > 1.0) {
    spdlog::error(
        "The \"rejection_ratio\" {} must be between 0.0 and 1.0",
        rejection_ratio);
    return false;
  }
  if (rgb_bands_map.empty()) {
    new_rgb_bands_map = {1, 2, 3};
  } else if (rgb_bands_map.size() == 3) {
    new_rgb_bands_map = rgb_bands_map;
  } else {
    spdlog::error(
        "The \"rgb_bands_map\" size {} must be 3 if not empty",
        rgb_bands_map.size());
    return false;
  }
  if (color_balancing) {
    auto names(color_balancing->ExportAllRastersName());
    if (auto it(std::find(
            names.begin(), names.end(), fs::path(path).filename().string()));
        it != names.end()) {
      color_balancing_idx = static_cast<int>(it - names.begin());
    }
  }
  return true;
}

void UpdateGeotrans(
    GDALDataset* source_raster_dataset,
    OGRSpatialReference* target_spatial_ref,
    double(&geotrans)[6]) {
  std::unique_ptr<void, void (*)(void*)> trans_arg(
      GDALCreateGenImgProjTransformer4(
          const_cast<OGRSpatialReference*>(
              source_raster_dataset->GetSpatialRef()), geotrans,
          target_spatial_ref, nullptr, nullptr),
      [](void* p) { GDALDestroyGenImgProjTransformer(p); });
  int x_size, y_size;
  GDALSuggestedWarpOutput(
      source_raster_dataset, GDALGenImgProjTransform, trans_arg.get(), geotrans,
      &x_size, &y_size);
}

}  // anonymous namespace

namespace rs_toolset {
namespace mosaicking {

MosaickingBase::MosaickingBase(double tol)
    : mem_driver_(GetGDALDriverManager()->GetDriverByName("MEM")),
      memory_driver_(GetGDALDriverManager()->GetDriverByName("Memory")),
      tol_(tol) {
  buffers_ = new int[10]{5, 10, 15, 20, 25, 30, 30, 30, 30, 30};
  spdlog::info("Creating a mosaicking base with\n - Tolerance: {}", tol);
}

bool MosaickingBase::RunTaskForSerial(
    const std::string& path,
    OGRLayer* mosaicking_layer,
    OGRLayer* border_layer,
    OGRGeometry* covered_border,
    std::vector<double>& borders_area,
    int low_overviews_trunc,
    int high_overviews_trunc,
    double surrounded_buffer,
    double rejection_ratio,
    bool anti_surrounded,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  GDALDatasetUniquePtr source_raster_dataset(nullptr);
  std::vector<int> new_rgb_bands_map;
  int color_balancing_idx(-1);
  if (!RunTaskForSerialCheck(
          path, mosaicking_layer, border_layer, covered_border, borders_area,
          low_overviews_trunc, high_overviews_trunc, surrounded_buffer,
          rejection_ratio, rgb_bands_map, color_balancing,
          source_raster_dataset, new_rgb_bands_map, color_balancing_idx)) {
    return false;
  }
  std::string string(
      "Running a mosaicking task for the existing mosaicking layer with {}\n"
      " - Low overviews trunction: {}\n"
      " - High overviews trunction: {}\n"
      " - Surrounded buffer: {}\n"
      " - Rejection ratio: {}\n"
      " - With anti-surrounded: {}\n"
      " - RGB bands' map: ");
  for (const auto& idx : new_rgb_bands_map)
    string.append(std::to_string(idx)).append(",");
  string.pop_back();
  spdlog::info(
      string, path, low_overviews_trunc, high_overviews_trunc,
      surrounded_buffer, rejection_ratio, anti_surrounded);

  double geotrans[6];
  source_raster_dataset->GetGeoTransform(geotrans);
  auto overviews_count(
      source_raster_dataset->GetRasterBand(1)->GetOverviewCount());
  OGRGeometryUniquePtr
      source_border(utils::CreateBorder(
          source_raster_dataset.get(), mosaicking_layer->GetSpatialRef())),
      original_source_border(source_border->clone());
  if (source_border->toPolygon()->getNumInteriorRings() != 0) {
    spdlog::warn(
        "Skipping adding the task "
        "since the source border has the internal hole(s)");
    return true;
  }
  if (OGRGeometryUniquePtr new_covered_polygon(source_border->clone());
      !source_border->Intersect(covered_border)) {
    spdlog::info(
        "Skipping updating the source border and the composite table layer "
        "in no intersection case");
    spdlog::debug("Updating the covered border");
    covered_border->toMultiPolygon()->addGeometryDirectly(
        new_covered_polygon.release());
    spdlog::info("Updating the covered border - done");
  } else {
    if (!source_raster_dataset->GetSpatialRef()->IsSame(
            mosaicking_layer->GetSpatialRef())) {
      UpdateGeotrans(
          source_raster_dataset.get(), mosaicking_layer->GetSpatialRef(),
          geotrans);
    }

    // Traverse all polygons in the covered border
    auto _covered_border(covered_border->toMultiPolygon());
    for (auto i(_covered_border->getNumGeometries() - 1); i >= 0; --i) {
      auto covered_polygon(_covered_border->getGeometryRef(i));
      if (!source_border->Intersect(covered_polygon)) continue;
      spdlog::debug("Initializing overlap geometries for the first loop");
      auto factor(pow(2, overviews_count - 1 - low_overviews_trunc));
      std::vector<OGRGeometryUniquePtr>
          covered_overlap_geometries, new_overlap_geometries, source_borders;
      Status status;
      if (!InitializeOverlapGeometries(
              factor, geotrans, surrounded_buffer, rejection_ratio,
              anti_surrounded, covered_polygon, source_border,
              covered_overlap_geometries, new_overlap_geometries,
              source_borders, status)) {
        return true;
      }
      for (auto j(0); j < source_borders.size(); ++j) {
        auto subtasks_count(1);
        if (status == Status::OVERLAP) {
          OGRGeometryUniquePtr geometry(source_borders[j]->Boundary());
          geometry.reset(geometry->Intersection(covered_polygon));
          if (geometry->getGeometryType() == wkbMultiLineString)
            subtasks_count = geometry->toMultiLineString()->getNumGeometries();
        }
        auto idx(overviews_count - 1 - low_overviews_trunc);
        factor = pow(2, idx);
        OGRGeometryUniquePtr valid_geometry(
            source_borders[j]->Intersection(covered_polygon));
        valid_geometry.reset(valid_geometry->Buffer(4 * geotrans[1]));
        GDALDatasetUniquePtr
            covered_overlap_dataset(nullptr),
            new_overlap_dataset(nullptr),
            label_raster_dataset(nullptr);
        spdlog::debug("Operating on the {} times downsampled overview", factor);
        CreateOverlapDatasets(
            factor, geotrans, new_rgb_bands_map, path, mosaicking_layer,
            covered_overlap_geometries[j].get(),
            new_overlap_geometries[j].get(), color_balancing,
            color_balancing_idx, covered_overlap_dataset, new_overlap_dataset,
            label_raster_dataset);
        double overlap_geotrans[6];
        covered_overlap_dataset->GetGeoTransform(overlap_geotrans);
        CreateSeamlines(
            overlap_geotrans, covered_overlap_dataset.get(),
            new_overlap_dataset.get(), label_raster_dataset,
            covered_overlap_geometries[j], new_overlap_geometries[j],
            true, status == Status::ANTI_SURROUNDED, subtasks_count);
        std::vector<std::pair<OGRGeometryUniquePtr, OGRGeometryUniquePtr>>
            subtasks;
        for (const auto& polygon :
                covered_overlap_geometries[j]->toMultiPolygon()) {
          subtasks.emplace_back(
              polygon->clone(), new_overlap_geometries[j]->clone());
        }
        for (auto& subtask : subtasks) {
          UpdateMediums(
              true, idx, idx == high_overviews_trunc, status == Status::OVERLAP,
              overlap_geotrans, valid_geometry.get(), subtask.first,
              subtask.second);
        }
        spdlog::info(
            "Operating on the {} times downsampled overview - done", factor);

        OGRGeometryUniquePtr
            covered_geometry(new OGRMultiPolygon),
            new_geometry(new OGRMultiPolygon);
        for (auto& subtask : subtasks) {
          spdlog::info("Adding a subtask with an intersected covered poylgon");
          idx = overviews_count - 1 - low_overviews_trunc - 1;
          factor = pow(2, idx);
          GDALDatasetUniquePtr sub_label_raster_dataset(mem_driver_->CreateCopy(
              "", label_raster_dataset.get(), true, nullptr, nullptr, nullptr));
          bool valid_task(true);
          for (; idx >= high_overviews_trunc; --idx, factor /= 2.0) {
            spdlog::debug(
                "Operating on the {} times downsampled overview", factor);
            CreateOverlapDatasets(
                factor, geotrans, new_rgb_bands_map, path, mosaicking_layer,
                subtask.first.get(), subtask.second.get(), color_balancing,
                color_balancing_idx, covered_overlap_dataset,
                new_overlap_dataset, sub_label_raster_dataset);
            covered_overlap_dataset->GetGeoTransform(overlap_geotrans);
            if (!CreateSeamlines(
                    overlap_geotrans, covered_overlap_dataset.get(),
                    new_overlap_dataset.get(), sub_label_raster_dataset,
                    subtask.first, subtask.second, false,
                    status == Status::ANTI_SURROUNDED)) {
              valid_task = false;
              break;
            }
            if (!UpdateMediums(
                    true, idx, idx == high_overviews_trunc,
                    status == Status::OVERLAP, overlap_geotrans,
                    valid_geometry.get(), subtask.first, subtask.second)) {
              UpdateMediums(
                    true, idx, true, status == Status::OVERLAP,
                    overlap_geotrans, valid_geometry.get(), subtask.first,
                    subtask.second);
              idx = high_overviews_trunc;
            }
            spdlog::info(
                "Operating on the {} times downsampled overview - done",
                factor);
          }
          if (valid_task && !covered_geometry->Intersect(subtask.first.get()) &&
              !new_geometry->Intersect(subtask.second.get())) {
            covered_geometry->toMultiPolygon()->addGeometryDirectly(
                subtask.first.release());
            new_geometry->toMultiPolygon()->addGeometryDirectly(
                subtask.second.release());
          } else {
            spdlog::info(
                "skipping running the task since "
                "the mosaicking result is invalid");
            continue;
          }
          spdlog::info(
              "Adding a subtask with an intersected covered poylgon - done");
        }
        switch (status) {
          case Status::OVERLAP: {
            UpdateResults(
                covered_geometry.get(), new_geometry.get(), covered_polygon,
                rejection_ratio, source_borders[j], mosaicking_layer,
                border_layer, borders_area);
            break;
          }
          case Status::SURROUNDED: {
            source_borders[j].reset(source_borders[j]->Difference(
                covered_geometry.get()));
            for (auto j(mosaicking_layer->GetFeatureCount() - 1); j >= 0; --j) {
              OGRFeatureUniquePtr feature(mosaicking_layer->GetFeature(j));
              if (auto geometry(feature->GetGeometryRef());
                  covered_geometry->Intersect(geometry)) {
                feature->SetGeometryDirectly(covered_geometry->Intersection(
                    geometry));
                mosaicking_layer->SetFeature(feature.get());
              } else if (source_borders[j]->Contains(geometry)) {
                auto last_fid(mosaicking_layer->GetFeatureCount() - 1);
                mosaicking_layer->DeleteFeature(j);
                border_layer->DeleteFeature(j);
                borders_area.erase(borders_area.begin() + j);
                if (j != last_fid) {
                  feature.reset(mosaicking_layer->GetFeature(last_fid));
                  feature->SetFID(j);
                  mosaicking_layer->DeleteFeature(last_fid);
                  mosaicking_layer->SetFeature(feature.get());
                  feature.reset(border_layer->GetFeature(last_fid));
                  feature->SetFID(j);
                  border_layer->DeleteFeature(last_fid);
                  border_layer->SetFeature(feature.get());
                }
              }
            }
            break;
          }
          case Status::ANTI_SURROUNDED: {
            source_borders[j] = std::move(new_geometry);
            for (auto k(mosaicking_layer->GetFeatureCount() - 1); k >= 0; --k) {
              OGRFeatureUniquePtr feature(mosaicking_layer->GetFeature(k));
              if (auto geometry(feature->GetGeometryRef());
                  covered_geometry->Intersect(geometry)) {
                feature->SetGeometryDirectly(geometry->Difference(
                    source_borders[j].get()));
                mosaicking_layer->SetFeature(feature.get());
              } else if (source_borders[j]->Contains(geometry)) {
                auto last_fid(mosaicking_layer->GetFeatureCount() - 1);
                mosaicking_layer->DeleteFeature(k);
                border_layer->DeleteFeature(k);
                borders_area.erase(borders_area.begin() + k);
                if (k != last_fid) {
                  feature.reset(mosaicking_layer->GetFeature(last_fid));
                  feature->SetFID(k);
                  mosaicking_layer->DeleteFeature(last_fid);
                  mosaicking_layer->SetFeature(feature.get());
                  feature.reset(border_layer->GetFeature(last_fid));
                  feature->SetFID(k);
                  border_layer->DeleteFeature(last_fid);
                  border_layer->SetFeature(feature.get());
                }
              }
            }
          }
        }
      }

      // Reconstruct source border
      if (source_border) {
        for (const auto& geometry : source_borders)
          source_border.reset(source_border->Union(geometry.get()));
      } else {
        source_border = std::move(source_borders[0]);
      }

      // Delete the processed covered polygon and union the new covered polygon
      new_covered_polygon.reset(new_covered_polygon->Union(covered_polygon));
      _covered_border->removeGeometry(i);
    }

    // Add the new covered polygon
    spdlog::debug("Updating the covered border");
    switch (new_covered_polygon->getGeometryType()) {
      case wkbPolygon: {
        auto _new_covered_polygon(new_covered_polygon->toPolygon());
        for (int i(_new_covered_polygon->getNumInteriorRings() - 1); i >= 0;
            --i) {
          if (_new_covered_polygon->getInteriorRing(i)->get_Area() == 0)
            _new_covered_polygon->removeRing(i + 1);
        }
        break;
      }
      case wkbMultiPolygon: {
        utils::ExtractBiggestPolygon(new_covered_polygon);
      }
    }
    _covered_border->addGeometryDirectly(new_covered_polygon.release());
    spdlog::info("Updating the covered border - done");
  }

  // Update the border layer
  spdlog::debug("Updating the border layer");
  OGRFeatureUniquePtr border_feature(OGRFeature::CreateFeature(
      border_layer->GetLayerDefn()));
  border_feature->SetGeometryDirectly(original_source_border.release());
  border_feature->SetField(0, path.c_str());
  border_feature->SetFID(border_layer->GetFeatureCount());
  border_layer->CreateFeature(border_feature.get());
  borders_area.push_back(
      border_feature->GetGeometryRef()->toPolygon()->get_Area());
  spdlog::info("Updating the border layer - done");

  OGRFeatureUniquePtr composite_table_feature(OGRFeature::CreateFeature(
      mosaicking_layer->GetLayerDefn()));
  composite_table_feature->SetField(0, path.c_str());
  if (color_balancing)
    composite_table_feature->SetField(1, color_balancing_idx);
  composite_table_feature->SetGeometryDirectly(source_border.release());
  composite_table_feature->SetFID(mosaicking_layer->GetFeatureCount());
  mosaicking_layer->CreateFeature(composite_table_feature.get());
  spdlog::info("Running a mosaicking task - done");
  return true;
}

OGRGeometryUniquePtr MosaickingBase::RunTaskForPair(
    const std::string& path1,
    const std::string& path2,
    OGRGeometry* border1,
    OGRGeometry* border2,
    OGRGeometry* geometry1,
    OGRGeometry* geometry2,
    int low_overviews_trunc,
    int high_overviews_trunc,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    int color_balancing_idx1,
    int color_balancing_idx2) {
  std::string string(
      "Running a mosaicking task for the raster pair with\n"
      " - Raster1 path: {}\n"
      " - Raster2 path: {}\n"
      " - Low overviews trunction: {}\n"
      " - High overviews trunction: {}\n"
      " - RGB bands' map: ");
  if (rgb_bands_map.empty()) {
    string.append("1,2,3");
  } else if (rgb_bands_map.size() != 3) {
    spdlog::error(
        "The size of \"rgb_bands_map\" {} must be 3", rgb_bands_map.size());
    return false;
  } else {
    for (const auto& idx : rgb_bands_map)
      string.append(std::to_string(idx)).append(",");
    string.pop_back();
  }
  spdlog::info(string,path1, path2, low_overviews_trunc, high_overviews_trunc);
  GDALDatasetUniquePtr
      dataset1(GDALDataset::Open(
          path1.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      dataset2(GDALDataset::Open(
          path2.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!dataset1) {
    spdlog::error("Opening {} failed", path1);
    return nullptr;
  }
  if (!dataset2) {
    spdlog::error("Opening {} failed", path2);
    return nullptr;
  }
  double geotrans1[6], geotrans2[6];
  if (dataset1->GetGeoTransform(geotrans1) != CE_None ||
      !dataset1->GetSpatialRef()) { 
    spdlog::error(
        "{} does not have the geotransform or the spatial reference. "
        "Please check whether the raster is DOM", path1);
    return nullptr;
  }
  if (dataset2->GetGeoTransform(geotrans2) != CE_None ||
      !dataset2->GetSpatialRef()) { 
    spdlog::error(
        "{} does not have the geotransform or the spatial reference. "
        "Please check whether the raster is DOM", path2);
    return nullptr;
  }
  if (!geometry1->getSpatialReference()) {
    spdlog::error("The \"geometry1\" does not have the spatial reference");
    return nullptr;
  }
  if (!geometry2->getSpatialReference()) {
    spdlog::error("The \"geometry2\" does not have the spatial reference");
    return nullptr;
  }
  if (!geometry1->getSpatialReference()->IsSame(
          geometry2->getSpatialReference())) {
    spdlog::error(
        "The \"geometry1\" and the \"geometry2\" do not have "
        "the same spatial reference");
    return nullptr;
  }
  std::unique_ptr<OGRSpatialReference> spatial_ref(
      geometry1->getSpatialReference()->Clone());
  utils::CreateRasterPyra(dataset1.get());
  utils::CreateRasterPyra(dataset2.get());
  auto overviews_count(std::min(
      dataset1->GetRasterBand(1)->GetOverviewCount(),
      dataset2->GetRasterBand(1)->GetOverviewCount()));
  if (low_overviews_trunc < 0 || high_overviews_trunc < 0) {
    spdlog::error(
        "The \"low_overviews_trunc\" {} and the \"high_overviews_trunc\" {} "
        "must be non-negative", low_overviews_trunc, high_overviews_trunc);
    return nullptr;
  }
  if ((low_overviews_trunc + high_overviews_trunc) >= overviews_count) {
    spdlog::warn(
        "The \"low_overviews_trunc\" {} plus the \"high_overviews_trunc\" {} "
        "must be less than the overviews' count {}", low_overviews_trunc,
        high_overviews_trunc, overviews_count);
    high_overviews_trunc = std::min(high_overviews_trunc, overviews_count - 1);
    low_overviews_trunc = overviews_count - 1 - high_overviews_trunc;
    spdlog::info(
        "Reset arguments:\n"
        " - Low overviews trunction: {}\n"
        " - High overviews trunction: {}", low_overviews_trunc,
        high_overviews_trunc);
  }

  // Update the geotransform after reprojection
  auto geotrans(geotrans1[1] > geotrans2[1] ? geotrans1 : geotrans2);
  auto dataset(geotrans1[1] > geotrans2[1] ? dataset1.get() : dataset2.get());
  std::unique_ptr<void, void (*)(void*)> trans_arg(
      GDALCreateGenImgProjTransformer4(
          const_cast<OGRSpatialReference*>(dataset->GetSpatialRef()), geotrans,
          spatial_ref.get(), nullptr, nullptr),
      [](void* p) { GDALDestroyGenImgProjTransformer(p); });
  int x_size, y_size;
  GDALSuggestedWarpOutput(
      dataset, GDALGenImgProjTransform, trans_arg.get(), geotrans, &x_size,
      &y_size);

  // Create overlap geometries for the first loop
  spdlog::debug("Initializing overlap geometries for the first loop");
  auto factor(pow(2, overviews_count - 1 - low_overviews_trunc));
  OGRGeometryUniquePtr seamline(geometry1->Intersection(geometry2));
  if (seamline->getGeometryType() == wkbMultiLineString)
    utils::JointMultiLineString(seamline);
  auto _seamline(seamline->toLineString());
  OGRPoint centroid1, centroid2, head_point, tail_point, point;
  geometry1->Centroid(&centroid1);
  geometry2->Centroid(&centroid2);
  _seamline->StartPoint(&head_point);
  _seamline->EndPoint(&tail_point);
  auto valid_linear_ring(new OGRLinearRing);
  valid_linear_ring->addPoint(&head_point);
  valid_linear_ring->addPoint(centroid1.getX(), centroid1.getY());
  valid_linear_ring->addPoint(&tail_point);
  valid_linear_ring->addPoint(centroid2.getX(), centroid2.getY());
  valid_linear_ring->addPoint(&head_point);
  OGRGeometryUniquePtr valid_geometry(new OGRPolygon);
  valid_geometry->toPolygon()->addRingDirectly(valid_linear_ring);
  OGRGeometryUniquePtr _border1(nullptr), _border2(nullptr);
  if (!border1) {
    _border1 = utils::CreateBorder(dataset1.get(), spatial_ref.get());
    border1 = _border1.get();
  }
  if (!border2) {
    _border2 = utils::CreateBorder(dataset2.get(), spatial_ref.get());
    border2 = _border2.get();
  }
  valid_geometry.reset(valid_geometry->Intersection(border1));
  valid_geometry.reset(valid_geometry->Intersection(border2));
  valid_geometry.reset(valid_geometry->Buffer(-factor * geotrans[1]));
  valid_linear_ring = valid_geometry->toPolygon()->getExteriorRing();
  std::pair closest_head(-1, DBL_MAX), closest_tail(-1, DBL_MAX);
  auto points_count(valid_linear_ring->getNumPoints());
  for (int i(0); i < points_count; ++i) {
    valid_linear_ring->getPoint(i, &point);
    if (auto dist(head_point.Distance(&point)); dist < closest_head.second)
      closest_head = {i, dist};
    if (auto dist(tail_point.Distance(&point)); dist < closest_tail.second)
      closest_tail = {i, dist};
  }
  if (closest_head.first == 0 || closest_head.first == points_count - 1) {
    valid_linear_ring->setPoint(0, &head_point);
    valid_linear_ring->setPoint(points_count - 1, &head_point);
  } else {
    valid_linear_ring->setPoint(closest_head.first, &head_point);
  }
  if (closest_tail.first == 0 || closest_tail.first == points_count - 1) {
    valid_linear_ring->setPoint(0, &tail_point);
    valid_linear_ring->setPoint(points_count - 1, &tail_point);
  } else {
    valid_linear_ring->setPoint(closest_tail.first, &tail_point);
  }
  OGRGeometryUniquePtr
      overlap1_geometry(valid_geometry->Union(geometry2)),
      overlap2_geometry(valid_geometry->Union(geometry1));
  spdlog::debug(
      "Initializing overlap geometries for the first loop - done");

  double overlap_geotrans[6];
  GDALDatasetUniquePtr
      overlap1_dataset(nullptr),
      overlap2_dataset(nullptr),
      label_raster_dataset(nullptr);
  for (int idx(overviews_count - 1 - low_overviews_trunc);
      idx >= high_overviews_trunc; --idx) {
    factor = pow(2, idx);
    spdlog::debug("Operating on the {} times downsampled overview", factor);
    CreateOverlapDatasets(
        factor, geotrans, rgb_bands_map, spatial_ref.get(), path1, path2,
        overlap1_geometry.get(), overlap2_geometry.get(), color_balancing,
        color_balancing_idx1, color_balancing_idx2, overlap1_dataset,
        overlap2_dataset, label_raster_dataset);
    overlap1_dataset->GetGeoTransform(overlap_geotrans);
    CreateSeamlines(
        overlap_geotrans, overlap1_dataset.get(), overlap2_dataset.get(),
        label_raster_dataset, overlap1_geometry, overlap2_geometry);
    UpdateMediums(
        false, idx, idx == high_overviews_trunc, true, overlap_geotrans,
        valid_geometry.get(), overlap1_geometry, overlap2_geometry);
    spdlog::info(
        "Operating on the {} times downsampled overview - done", factor);
  }
  seamline.reset(overlap1_geometry->Intersection(overlap2_geometry.get()));
  if (seamline->getGeometryType() == wkbMultiLineString)
    utils::JointMultiLineString(seamline);
  _seamline = seamline->toLineString();
  _seamline->StartPoint(&point);
  while (!valid_geometry->Contains(&point)) {
    _seamline->removePoint(0);
    _seamline->StartPoint(&point);
  }
  _seamline->EndPoint(&point);
  while (!valid_geometry->Contains(&point)) {
    _seamline->removePoint(_seamline->getNumPoints() - 1);
    _seamline->EndPoint(&point);
  }
  if (point.Distance(&head_point) > point.Distance(&tail_point)) {
    _seamline->addPoint(&tail_point);
    _seamline->reversePoints();
    _seamline->addPoint(&head_point);
  } else {
    _seamline->addPoint(&head_point);
    _seamline->reversePoints();
    _seamline->addPoint(&tail_point);
  }
  return seamline;
}

bool MosaickingBase::InitializeOverlapGeometries(
    double factor,
    double* geotrans,
    double surrounded_buffer,
    double rejection_ratio,
    bool anti_surrounded,
    OGRPolygon* covered_polygon,
    OGRGeometryUniquePtr& source_border,
    std::vector<OGRGeometryUniquePtr>& covered_overlap_geometries,
    std::vector<OGRGeometryUniquePtr>& new_overlap_geometries,
    std::vector<OGRGeometryUniquePtr>& source_borders,
    Status& status) {
  OGRGeometryUniquePtr covered_outer_polygon(new OGRPolygon);
  covered_outer_polygon->toPolygon()->addRingDirectly(
      covered_polygon->getExteriorRing()->clone());
  if (source_border->Contains(covered_outer_polygon.get())) {
    spdlog::debug("Initializing overlap geometries in the surrounded case");
    status = Status::SURROUNDED;
    covered_overlap_geometries.emplace_back(covered_polygon->clone());
    OGRGeometryUniquePtr geometry(covered_overlap_geometries[0]->Buffer(
        -surrounded_buffer * factor * geotrans[1]));
    new_overlap_geometries.emplace_back(covered_overlap_geometries[0]->Buffer(
        surrounded_buffer * factor * geotrans[1]));
    new_overlap_geometries[0].reset(source_border->Intersection(
        new_overlap_geometries[0].get()));
    new_overlap_geometries[0].reset(new_overlap_geometries[0]->Difference(
        geometry.get()));
  } else if (covered_outer_polygon->Contains(source_border.get())) {
    if (covered_polygon->Contains(source_border.get()) && !anti_surrounded) {
      spdlog::info(
          "Skipping adding the task "
          "since the covered border contains the source border");
      return false;
    }
    spdlog::debug(
        "Initializing overlap geometries in the anti-surrounded case");
    status = Status::ANTI_SURROUNDED;
    new_overlap_geometries.emplace_back(source_border->clone());
    OGRGeometryUniquePtr geometry(new_overlap_geometries[0]->Buffer(
        -surrounded_buffer * factor * geotrans[1]));
    covered_overlap_geometries.emplace_back(new_overlap_geometries[0]->Buffer(
        surrounded_buffer * factor * geotrans[1]));
    covered_overlap_geometries[0].reset(covered_polygon->Intersection(
        covered_overlap_geometries[0].get()));
    covered_overlap_geometries[0].reset(
        covered_overlap_geometries[0]->Difference(geometry.get()));
  } else {
    spdlog::debug("Initializing overlap geometries in the overlap case");
    status = Status::OVERLAP;
    switch (OGRGeometryUniquePtr geometry(source_border->Intersection(
            covered_polygon)); geometry->getGeometryType()) {
      case wkbPolygon: {
        if (double percent(geometry->toPolygon()->get_Area() /
                source_border->toPolygon()->get_Area());
            percent > (1 - rejection_ratio)) {
          spdlog::warn(
              "Skipping adding the task since the intersected area "
              "divided by the source border area is {}", percent);
          return false;
        }
        covered_overlap_geometries.emplace_back(geometry->Buffer(
            20 * factor * geotrans[1]));
        new_overlap_geometries.emplace_back(
          covered_overlap_geometries[0]->clone());
        break;
      }
      case wkbMultiPolygon: {
        for (const auto& subgeometry : geometry->toMultiPolygon()) {
          covered_overlap_geometries.emplace_back(subgeometry->Buffer(
              20 * factor * geotrans[1]));
          new_overlap_geometries.emplace_back(
              covered_overlap_geometries.back()->clone());
          source_borders.emplace_back(subgeometry->Buffer(
              0.1 * factor * geotrans[1]));
          source_border.reset(source_border->Difference(
              source_borders.back().get()));
          source_borders.back().reset(source_borders.back()->Buffer(
              0.1 * factor * geotrans[1]));
          source_borders.back().reset(source_borders.back()->Difference(
              source_border.get()));
        }
      }
    }
  }
  if (source_borders.empty())
    source_borders.push_back(std::move(source_border));
  spdlog::debug("Initializing overlap geometries for the first loop - done");
  return true;
}

void MosaickingBase::CreateOverlapDatasets(
    double factor,
    double* geotrans,
    const std::vector<int>& rgb_bands_map,
    const std::string& source_raster_path,
    OGRLayer* mosaicking_layer,
    OGRGeometry* covered_overlap_geometry,
    OGRGeometry* new_overlap_geometry,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    int color_balancing_idx,
    GDALDatasetUniquePtr& covered_overlap_dataset,
    GDALDatasetUniquePtr& new_overlap_dataset,
    GDALDatasetUniquePtr& label_raster_dataset) {
  spdlog::debug("Updating overlap datasets");

  // Create the union range and geotransform
  OGRGeometryUniquePtr union_geometry(covered_overlap_geometry->Union(
      new_overlap_geometry));
  OGREnvelope enve;
  union_geometry->getEnvelope(&enve);
  double inv_geotrans[6], union_geotrans[6], x, y;
  GDALInvGeoTransform(geotrans, inv_geotrans);
  memset(union_geotrans, 0, 6 * sizeof(double));
  union_geotrans[1] = factor * geotrans[1];
  union_geotrans[5] = factor * geotrans[5];
  int range[4];
  GDALApplyGeoTransform(inv_geotrans, enve.MinX, enve.MaxY, &x, &y);
  range[0] = static_cast<int>(floor(x));
  range[1] = static_cast<int>(floor(y));
  GDALApplyGeoTransform(inv_geotrans, enve.MaxX, enve.MinY, &x, &y);
  range[2] = static_cast<int>(ceil((x - range[0]) / factor));
  range[3] = static_cast<int>(ceil((y - range[1]) / factor));
  GDALApplyGeoTransform(
      geotrans, range[0], range[1], union_geotrans, union_geotrans + 3);
  auto spatial_ref(mosaicking_layer->GetSpatialRef());

  covered_overlap_dataset.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  covered_overlap_dataset->SetGeoTransform(union_geotrans);
  covered_overlap_dataset->SetSpatialRef(spatial_ref);
  new_overlap_dataset.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  new_overlap_dataset->SetGeoTransform(union_geotrans);
  new_overlap_dataset->SetSpatialRef(spatial_ref);

  // Update the covered overlap dataset
  std::vector<int> color_balanced_rasters_idx;
  std::vector<std::string> original_rasters_path;
  std::vector<OGRGeometry*>
      color_balancing_intersecting_geometries,
      original_intersecting_geometries;
  for (const auto& feature : mosaicking_layer) {
    if (feature->GetGeometryRef()->Intersect(covered_overlap_geometry)) {
      auto geoemtry(feature->GetGeometryRef()->Intersection(
          covered_overlap_geometry));
      geoemtry->assignSpatialReference(spatial_ref);
      if (color_balancing && feature->GetFieldAsInteger(1) != -1) {
        color_balanced_rasters_idx.push_back(feature->GetFieldAsInteger(1));
        color_balancing_intersecting_geometries.push_back(geoemtry);
      } else {
        original_rasters_path.push_back(feature->GetFieldAsString(0));
        original_intersecting_geometries.push_back(geoemtry);
      }
    }
  }
  if (!color_balanced_rasters_idx.empty()) {
    color_balancing->WarpByGeometry(
        color_balanced_rasters_idx, color_balancing_intersecting_geometries,
        covered_overlap_dataset, rgb_bands_map);
  }
  for (auto& geoemtry : color_balancing_intersecting_geometries)
    OGRGeometryFactory::destroyGeometry(geoemtry);
  if (!original_rasters_path.empty()) {
    utils::WarpByGeometry(
        original_rasters_path, original_intersecting_geometries,
        covered_overlap_dataset, rgb_bands_map);
  }
  for (auto& geometry : original_intersecting_geometries)
    OGRGeometryFactory::destroyGeometry(geometry);
  spdlog::debug("Updating the covered overlap dataset - done");

  // Update the new overlap dataset
  new_overlap_geometry->assignSpatialReference(spatial_ref);
  if (color_balancing && color_balancing_idx != -1) {
    color_balancing->WarpByGeometry(
        {color_balancing_idx}, {new_overlap_geometry}, new_overlap_dataset,
        rgb_bands_map);
  } else {
    utils::WarpByGeometry(
        {source_raster_path}, {new_overlap_geometry}, new_overlap_dataset,
        rgb_bands_map);
  }
  spdlog::debug("Updating the new overlap dataset - done");

  // Update the label raster dataset if necessary
  if (label_raster_dataset) {
    GDALDatasetUniquePtr dataset(mem_driver_->Create(
        "", range[2], range[3], 1, GDT_Byte, nullptr));
    dataset->SetGeoTransform(union_geotrans);
    dataset->SetSpatialRef(spatial_ref);
    utils::WarpByGeometry(
        {label_raster_dataset.get()}, {nullptr}, dataset, {},
        GRA_NearestNeighbour);
    label_raster_dataset.swap(dataset);
    spdlog::debug("Updating the label raster dataset - done");
  }
  spdlog::debug("Updating overlap datasets - done");
}

void MosaickingBase::CreateOverlapDatasets(
    double factor,
    double* geotrans,
    const std::vector<int>& rgb_bands_map,
    OGRSpatialReference* spatial_ref,
    const std::string& raster_path1,
    const std::string& raster_path2,
    OGRGeometry* overlap_geometry1,
    OGRGeometry* overlap_geometry2,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    int color_balancing_idx1,
    int color_balancing_idx2,
    GDALDatasetUniquePtr& overlap_dataset1,
    GDALDatasetUniquePtr& overlap_dataset2,
    GDALDatasetUniquePtr& label_raster_dataset) {
  spdlog::debug("Updating overlap datasets");

  // Create the union range and geotransform
  OGRGeometryUniquePtr union_geometry(overlap_geometry1->Union(
      overlap_geometry2));
  OGREnvelope enve;
  union_geometry->getEnvelope(&enve);
  double inv_geotrans[6], union_geotrans[6], x, y;
  GDALInvGeoTransform(geotrans, inv_geotrans);
  memset(union_geotrans, 0, 6 * sizeof(double));
  union_geotrans[1] = factor * geotrans[1];
  union_geotrans[5] = factor * geotrans[5];
  int range[4];
  GDALApplyGeoTransform(inv_geotrans, enve.MinX, enve.MaxY, &x, &y);
  range[0] = static_cast<int>(floor(x));
  range[1] = static_cast<int>(floor(y));
  GDALApplyGeoTransform(inv_geotrans, enve.MaxX, enve.MinY, &x, &y);
  range[2] = static_cast<int>(ceil((x - range[0]) / factor));
  range[3] = static_cast<int>(ceil((y - range[1]) / factor));
  GDALApplyGeoTransform(
      geotrans, range[0], range[1], union_geotrans, union_geotrans + 3);

  overlap_dataset1.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  overlap_dataset1->SetGeoTransform(union_geotrans);
  overlap_dataset1->SetSpatialRef(spatial_ref);
  overlap_dataset2.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  overlap_dataset2->SetGeoTransform(union_geotrans);
  overlap_dataset2->SetSpatialRef(spatial_ref);

  // Update the overlap1 dataset
  overlap_geometry1->assignSpatialReference(spatial_ref);
  if (color_balancing && color_balancing_idx1 != -1) {
    color_balancing->WarpByGeometry(
        {color_balancing_idx1}, {overlap_geometry1}, overlap_dataset1,
        rgb_bands_map);
  } else {
    utils::WarpByGeometry(
        {raster_path1}, {overlap_geometry1}, overlap_dataset1, rgb_bands_map);
  }
  spdlog::debug("Updating the overlap1 dataset - done");

  // Update the overlap2 dataset
  overlap_geometry2->assignSpatialReference(spatial_ref);
  if (color_balancing && color_balancing_idx2 != -1) {
    color_balancing->WarpByGeometry(
        {color_balancing_idx2}, {overlap_geometry2}, overlap_dataset2,
        rgb_bands_map);
  } else {
    utils::WarpByGeometry(
        {raster_path2}, {overlap_geometry2}, overlap_dataset2, rgb_bands_map);
  }
  spdlog::debug("Updating the overlap2 dataset - done");

  // Update the label raster dataset if necessary
  if (label_raster_dataset) {
    GDALDatasetUniquePtr dataset(mem_driver_->Create(
        "", range[2], range[3], 1, GDT_Byte, nullptr));
    dataset->SetGeoTransform(union_geotrans);
    dataset->SetSpatialRef(spatial_ref);
    utils::WarpByGeometry(
        {label_raster_dataset.get()}, {nullptr}, dataset, {},
        GRA_NearestNeighbour);
    label_raster_dataset.swap(dataset);
    spdlog::debug("Updating the label raster dataset - done");
  }
  spdlog::debug("Updating overlap datasets - done");
}

bool MosaickingBase::CreateSeamlines(
    double* geotrans,
    GDALDataset* covered_overlap_dataset,
    GDALDataset* new_overlap_dataset,
    GDALDatasetUniquePtr& label_raster_dataset,
    OGRGeometryUniquePtr& covered_overlap_geometry,
    OGRGeometryUniquePtr& new_overlap_geometry,
    bool connection_analysis,
    bool swap,
    int subtasks_count) {
  cv::Mat covered_mat, new_mat, label_mat;
  PrepareData(
      covered_overlap_dataset, new_overlap_dataset, label_raster_dataset.get(),
      covered_mat, new_mat, label_mat, connection_analysis);
  return ExecuteMosaicking(
      covered_mat, new_mat, label_mat, geotrans,
      const_cast<OGRSpatialReference*>(new_overlap_dataset->GetSpatialRef()),
      label_raster_dataset, covered_overlap_geometry, new_overlap_geometry,
      swap, subtasks_count);
}

bool MosaickingBase::UpdateMediums(
    bool serial,
    int idx,
    bool last,
    bool overlap,
    double* geotrans,
    OGRGeometry* valid_geometry,
    OGRGeometryUniquePtr& covered_overlap_geometry,
    OGRGeometryUniquePtr& new_overlap_geometry) {
  OGRGeometryUniquePtr seamline(
      covered_overlap_geometry->Intersection(new_overlap_geometry.get()));
  if (seamline->getGeometryType() == wkbMultiLineString)
    utils::JointMultiLineString(seamline);

  if (last) {
    if (tol_) {
      if (auto _new_overlap_geometry(new_overlap_geometry->toPolygon());
          _new_overlap_geometry->getNumInteriorRings() == 1) {
        covered_overlap_geometry.reset(
            covered_overlap_geometry->Simplify(tol_ * geotrans[1]));
        _new_overlap_geometry->removeRing(1);
        _new_overlap_geometry->addRingDirectly(
            covered_overlap_geometry->toPolygon()->getExteriorRing()->clone());
      }
      else if (auto _covered_overlap_geometry(
              covered_overlap_geometry->toPolygon());
          _covered_overlap_geometry->getNumInteriorRings() == 1) {
        new_overlap_geometry.reset(
            new_overlap_geometry->Simplify(tol_ * geotrans[1]));
        _covered_overlap_geometry->removeRing(1);
        _covered_overlap_geometry->addRingDirectly(
            new_overlap_geometry->toPolygon()->getExteriorRing()->clone());
      } else {
        utils::RearrangeNeighborGeometries(
            seamline.get(), covered_overlap_geometry, new_overlap_geometry);
        utils::GraftSeamline(
            seamline, covered_overlap_geometry, new_overlap_geometry,
            tol_ * geotrans[1]);
        if (serial) {
          OGRGeometryUniquePtr 
              buffered_covered_geometry(covered_overlap_geometry->Simplify(
                  5.0 * geotrans[1])),
              buffered_new_geometry(new_overlap_geometry->Simplify(
                  5.0 * geotrans[1]));
          buffered_covered_geometry.reset(buffered_covered_geometry->Buffer(
              5.0 * geotrans[1])),
          buffered_new_geometry.reset(buffered_new_geometry->Buffer(
              5.0 * geotrans[1]));
          buffered_covered_geometry.reset(buffered_covered_geometry->Difference(
              new_overlap_geometry.get()));
          buffered_new_geometry.reset(buffered_new_geometry->Difference(
              covered_overlap_geometry.get()));
          if (buffered_covered_geometry->getGeometryType() == wkbMultiPolygon)
            utils::ExtractBiggestPolygon(buffered_covered_geometry);
          if (buffered_new_geometry->getGeometryType() == wkbMultiPolygon)
            utils::ExtractBiggestPolygon(buffered_new_geometry);
          buffered_covered_geometry.swap(covered_overlap_geometry);
          buffered_new_geometry.swap(new_overlap_geometry);
        }
      }
    }
    return true;
  } else {
    spdlog::debug("Updating the overlap geometries");
    seamline.reset(seamline->Simplify(5 * geotrans[1]));
    double buffer_at_end(serial ? 2.0 : 0.0);
    OGRGeometryUniquePtr geometry1(nullptr), geometry2(nullptr);
    do {
      // Create the seamline polygon and the buffered seamline polygon
      OGRGeometryUniquePtr intersection_geometry(seamline->Buffer(
          (buffers_[idx] - buffer_at_end) * geotrans[1]));
      intersection_geometry.reset(intersection_geometry->Simplify(
          5 * geotrans[1]));
      intersection_geometry.reset(intersection_geometry->Intersection(
          valid_geometry));
      OGRGeometryUniquePtr union_geometry(intersection_geometry->Buffer(
          (buffers_[idx] - buffer_at_end) * geotrans[1]));

      // Update overlap geometries
      geometry1.reset(covered_overlap_geometry->Union(
          intersection_geometry.get()));
      geometry2.reset(new_overlap_geometry->Union(intersection_geometry.get()));
      geometry1.reset(geometry1->Intersection(union_geometry.get()));
      geometry2.reset(geometry2->Intersection(union_geometry.get()));
      if (buffer_at_end) {
        geometry1.reset(geometry1->Buffer(buffer_at_end * geotrans[1]));
        geometry2.reset(geometry2->Buffer(buffer_at_end * geotrans[1]));
      }
      if (!overlap || (geometry1->toPolygon()->getNumInteriorRings() == 0 &&
          geometry2->toPolygon()->getNumInteriorRings() == 0)) {
        covered_overlap_geometry.swap(geometry1);
        new_overlap_geometry.swap(geometry2);
        spdlog::debug("Updating the overlap geometries - done");
        return true;
      }
      idx--;
    } while (idx >= 0);
    spdlog::info(
        "Stopping the operation early since the overlap geometreis still have "
        "the internal hole(s) with the least buffer {}", buffers_[0]);
    return false;
  }
}

void MosaickingBase::UpdateResults(
    OGRGeometry* covered_geometry,
    OGRGeometry* new_geometry,
    OGRGeometry* covered_polygon,
    double rejection_ratio,
    OGRGeometryUniquePtr& source_border,
    OGRLayer* composite_table_layer,
    OGRLayer* border_layer,
    std::vector<double>& borders_area) {
  // Update the source border
  spdlog::debug("Updating the source border");
  switch (OGRGeometryUniquePtr diff_geometry(source_border->Difference(
          covered_geometry)); diff_geometry->getGeometryType()) {
    case wkbPolygon: {
      source_border.reset(diff_geometry.release());
      break;
    }
    case wkbMultiPolygon: {
      source_border.reset(source_border->Difference(covered_polygon));
      if (source_border->getGeometryType() == wkbMultiPolygon)
        utils::ExtractBiggestPolygon(source_border);
      for (const auto& geometry : diff_geometry->toMultiPolygon()) {
        if (geometry->Intersect(new_geometry)) {
          source_border.reset(source_border->Union(geometry));
          break;
        }
      }
    }
  }
  spdlog::info("Updating the source border - done");

  // Update the composite table
  spdlog::debug("Updating the composite table");
  std::vector<std::pair<int, OGRGeometryUniquePtr>> dislocated_infos;
  for (auto i(static_cast<int>(composite_table_layer->GetFeatureCount() - 1));
      i >= 0; --i) {
    OGRFeatureUniquePtr feature(composite_table_layer->GetFeature(i));
    if (auto geometry(feature->GetGeometryRef());
        covered_geometry->Intersect(geometry)) {
      OGRGeometryUniquePtr diff_geometry(geometry->Difference(
          source_border.get()));
      switch (diff_geometry->getGeometryType()) {
        case wkbPolygon: {
          dislocated_infos.emplace_back(i, std::move(diff_geometry));
          break;
        }
        case wkbMultiPolygon: {
          bool b(true);
          double max_area(0.0);
          for (const auto& polygon : diff_geometry->toMultiPolygon()) {
            dislocated_infos.emplace_back(i, polygon->clone());
            if (auto area(polygon->get_Area()); !b && area < max_area) {
              dislocated_infos[dislocated_infos.size() - 2].second.swap(
                  dislocated_infos.back().second);
            } else {
              max_area = area;
            }
            b = false;
          }
        }
      }
      feature->SetGeometryDirectly(dislocated_infos.back().second.release());
      composite_table_layer->SetFeature(feature.get());
      dislocated_infos.pop_back();
    } else if (source_border->Contains(geometry)) {
      auto last_fid(composite_table_layer->GetFeatureCount() - 1);
      composite_table_layer->DeleteFeature(i);
      border_layer->DeleteFeature(i);
      borders_area.erase(borders_area.begin() + i);
      if (i != last_fid) {
        feature.reset(composite_table_layer->GetFeature(last_fid));
        feature->SetFID(i);
        composite_table_layer->DeleteFeature(last_fid);
        composite_table_layer->SetFeature(feature.get());
        feature.reset(border_layer->GetFeature(last_fid));
        feature->SetFID(i);
        border_layer->DeleteFeature(last_fid);
        border_layer->SetFeature(feature.get());
      }
    }
  }
  std::vector<std::pair<int, OGRGeometryUniquePtr>> new_infos;
  for (auto& info : dislocated_infos) {
    bool b(true);
    for (int i(0); i < composite_table_layer->GetFeatureCount(); ++i) {
      if (OGRFeatureUniquePtr
              composite_table_feature(composite_table_layer->GetFeature(i)),
              border_feature(border_layer->GetFeature(i));
          info.second->Intersect(composite_table_feature->GetGeometryRef()) &&
          border_feature->GetGeometryRef()->Contains(info.second.get())) {
        if (OGRGeometryUniquePtr geometry(
                info.second->Union(composite_table_feature->GetGeometryRef()));
            geometry->getGeometryType() == wkbPolygon) {
          composite_table_feature->SetGeometryDirectly(geometry.release());
          composite_table_layer->SetFeature(composite_table_feature.get());
          b = false;
          break;
        }
      }
    }
    if (b)
      new_infos.push_back(std::move(info));
  }
  for (auto& info : new_infos) {
    if (auto area(info.second->toPolygon()->get_Area());
        area > rejection_ratio * borders_area[info.first]) {
      OGRFeatureUniquePtr
          composite_table_feature(composite_table_layer->GetFeature(
              info.first)),
          border_feature(border_layer->GetFeature(info.first));
      composite_table_feature->SetGeometryDirectly(info.second.release());
      composite_table_feature->SetFID(composite_table_layer->GetFeatureCount());
      border_feature->SetFID(border_layer->GetFeatureCount());
      composite_table_layer->CreateFeature(composite_table_feature.get());
      border_layer->CreateFeature(border_feature.get());
      borders_area.push_back(area);
    }
  }
  spdlog::info("Updating the composite table - done");
}

}  // namespace mosaicking
}  // namespace rs_toolset