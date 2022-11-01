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

namespace rs_toolset {
namespace mosaicking {

MosaickingBase::MosaickingBase(double tol)
    : mem_driver_(GetGDALDriverManager()->GetDriverByName("MEM")),
      memory_driver_(GetGDALDriverManager()->GetDriverByName("Memory")),
      tol_(tol) {
  buffers_ = new int[10]{0, 10, 15, 20, 25, 30, 30, 30, 30, 30};
  spdlog::info("Creating a mosaicking base with\n - Tolerance: {}", tol);
}

bool MosaickingBase::RunTaskForExisting(
    const std::string& path,
    OGRLayer* composite_table_layer,
    OGRLayer* border_layer,
    OGRGeometry* covered_border,
    std::vector<double>& borders_area,
    double rejection_ratio,
    int low_overview_trunc,
    int high_overview_trunc,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  std::string string(
      "Running a mosaicking task for the existing composite table with {}\n"
      " - Low overview trunction: {}\n"
      " - High overview trunction: {}\n"
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
  spdlog::info(string, path, low_overview_trunc, high_overview_trunc);
  GDALDatasetUniquePtr source_raster_dataset(GDALDataset::Open(
      path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!source_raster_dataset) {
    spdlog::error("Opening {} failed", path);
    return false;
  }
  double geotrans[6];
  if (source_raster_dataset->GetGeoTransform(geotrans) != CE_None ||
      !source_raster_dataset->GetSpatialRef()) { 
    spdlog::error(
        "{} does not have the geotransform or the spatial reference. "
        "Please check whether the raster is DOM", path);
    return false;
  }
  if (!composite_table_layer || !border_layer || !covered_border) {
    spdlog::error(
        "One of the \"composite_table_layer\", the \"border_layer\" and "
        "the \"covered_border\" is empty");
    return false;
  }
  if (covered_border->getGeometryType() != wkbMultiPolygon) {
    spdlog::error(
        "The \"covered_border\" geometry type {} must be 6(wkbMultiPolygon)",
        covered_border->getGeometryType());
    return false;
  }
  if (borders_area.size() != border_layer->GetFeatureCount()) {
    spdlog::error(
        "The size of \"borders_area\" {} must be equal to the features' count "
        "{} in the \"border_layer\"", borders_area.size(),
        border_layer->GetFeatureCount());
    return false;
  }
  int color_balancing_idx(-1);
  if (color_balancing) {
    auto names(color_balancing->ExportAllRastersName());
    if (auto it(std::find(
            names.begin(), names.end(), fs::path(path).filename().string()));
        it != names.end()) {
      color_balancing_idx = static_cast<int>(it - names.begin());
    }
  }
  utils::CreateRasterPyra(source_raster_dataset.get());
  auto overviews_count(
      source_raster_dataset->GetRasterBand(1)->GetOverviewCount());
  if (low_overview_trunc < 0 || high_overview_trunc < 0) {
    spdlog::error(
        "The low overview trunction {} and the high overview trunction {} "
        "must be non-negative", low_overview_trunc, high_overview_trunc);
    return false;
  }
  if ((low_overview_trunc + high_overview_trunc) >= overviews_count) {
    spdlog::warn(
        "The low overview trunction {} add the high overview trunction {} "
        "must be less than the overviews' count {}", low_overview_trunc,
        high_overview_trunc, overviews_count);
    high_overview_trunc = std::min(high_overview_trunc, overviews_count - 1);
    low_overview_trunc = overviews_count - 1 - high_overview_trunc;
    spdlog::info(
        "Reset arguments:\n"
        " - Low overview trunction: {}\n"
        " - High overview trunction: {}", low_overview_trunc,
        high_overview_trunc);
  }

  OGRGeometryUniquePtr
      source_border(utils::CreateBorder(
          source_raster_dataset.get(),
          composite_table_layer->GetSpatialRef())),
      original_source_border(source_border->clone());
  if (source_border->toPolygon()->getNumInteriorRings() != 0) {
    spdlog::warn(
        "Skipping adding the task "
        "since the source border has the internal hole(s)");
    return true;
  }
  if (OGRGeometryUniquePtr new_covered_polygon(source_border->clone());
      covered_border->IsEmpty() || !source_border->Intersect(covered_border)) {
    spdlog::info(
        "Skipping updating the source border and the composite table layer "
        "for the no intersection situation");
    spdlog::debug("Updating the covered border");
    covered_border->toMultiPolygon()->addGeometryDirectly(
        new_covered_polygon.release());
    spdlog::info("Updating the covered border - done");
  } else if (covered_border->Contains(source_border.get())) {
    spdlog::info(
        "Skipping adding the task "
        "since the covered border contains the source border");
    return true;
  } else {
    // Traverse all polygons in the covered border
    auto _covered_border(covered_border->toMultiPolygon());
    for (int i(_covered_border->getNumGeometries() - 1); i >= 0; i--) {
      // Continue the loop for the no intersection situation
      auto covered_polygon(_covered_border->getGeometryRef(i));
      if (!source_border->Intersect(covered_polygon)) continue;

      spdlog::info("Adding a subtask with an intersected covered poylgon");

      // Update the geotransform after reprojection
      std::unique_ptr<void, void (*)(void*)> trans_arg(
          GDALCreateGenImgProjTransformer4(
              const_cast<OGRSpatialReference*>(
                  source_raster_dataset->GetSpatialRef()), geotrans,
              composite_table_layer->GetSpatialRef(), nullptr, nullptr),
          [](void* p) { GDALDestroyGenImgProjTransformer(p); });
      int x_size, y_size;
      GDALSuggestedWarpOutput(
          source_raster_dataset.get(), GDALGenImgProjTransform,
          trans_arg.get(), geotrans, &x_size, &y_size);

      // Create overlap geometries for the first loop
      spdlog::debug("Initializing overlap geometries for the first loop");
      bool overlap_case(true);
      auto factor(pow(2, overviews_count - 1 - low_overview_trunc));
      OGRGeometryUniquePtr
          covered_overlap_geometry(nullptr),
          new_overlap_geometry(nullptr);
      if (source_border->Contains(covered_polygon)) {
        spdlog::debug(
            "Initializing overlap geometries in the all-surrounded case");
        overlap_case = false;
        covered_overlap_geometry.reset(covered_polygon->clone());
        new_overlap_geometry.reset(covered_overlap_geometry->Buffer(
            -20 * factor * geotrans[1]));
        new_overlap_geometry.reset(source_border->Difference(
            new_overlap_geometry.get()));
      } else {
        spdlog::debug("Initializing overlap geometries in the overlap case");
        covered_overlap_geometry.reset(source_border->Intersection(
            covered_polygon));
        covered_overlap_geometry.reset(covered_overlap_geometry->Buffer(
            20 * factor * geotrans[1]));
        new_overlap_geometry.reset(covered_overlap_geometry->clone());
      }
      spdlog::debug(
          "Initializing overlap geometries for the first loop - done");

      OGRGeometryUniquePtr valid_geometry(
          source_border->Intersection(covered_polygon));
      valid_geometry.reset(valid_geometry->Buffer(4 * geotrans[1]));
      GDALDatasetUniquePtr
          covered_overlap_dataset(nullptr),
          new_overlap_dataset(nullptr),
          label_raster_dataset(nullptr);
      for (int idx(overviews_count - 1 - low_overview_trunc);
          idx >= high_overview_trunc; idx--) {
        factor = pow(2, idx);
        spdlog::debug("Operating on the {} times downsampled overview", factor);
        CreateOverlapDatasets(
            factor, geotrans, rgb_bands_map, source_raster_dataset.get(),
            composite_table_layer, covered_overlap_geometry.get(),
            new_overlap_geometry.get(), color_balancing, color_balancing_idx,
            covered_overlap_dataset, new_overlap_dataset, label_raster_dataset);
        UpdateMediums(
            idx, covered_overlap_dataset.get(), new_overlap_dataset.get(),
            valid_geometry.get(), label_raster_dataset,
            covered_overlap_geometry, new_overlap_geometry,
            idx == high_overview_trunc, 2.0);
        spdlog::info(
            "Operating on the {} times downsampled overview - done", factor);
      }
      if (overlap_case) {
        UpdateResults(
            covered_overlap_geometry.get(), new_overlap_geometry.get(),
            covered_polygon, rejection_ratio, source_border,
            composite_table_layer, border_layer, borders_area);
      } else {
        source_border.reset(source_border->Difference(
            covered_overlap_geometry.get()));
        for (auto& feature : composite_table_layer) {
          if (covered_overlap_geometry->Intersect(feature->GetGeometryRef())) {
            feature->SetGeometryDirectly(covered_overlap_geometry->Intersection(
                feature->GetGeometryRef()));
            composite_table_layer->SetFeature(feature.get());
          }
        }
      }

      // Delete the processed covered polygon and union the new covered polygon
      new_covered_polygon.reset(new_covered_polygon->Union(covered_polygon));
      _covered_border->removeGeometry(i);
      spdlog::info(
          "Adding a subtask with an intersected covered poylgon - done");
    }

    // Add the new covered polygon
    spdlog::debug("Updating the covered border");
    _covered_border->addGeometryDirectly(new_covered_polygon.release());
    spdlog::info("Updating the covered border - done");
  }

  // Update the border layer
  spdlog::debug("Updating the border layer");
  OGRFeatureUniquePtr border_feature(OGRFeature::CreateFeature(
      border_layer->GetLayerDefn()));
  border_feature->SetGeometryDirectly(original_source_border.release());
  border_layer->CreateFeature(border_feature.get());
  borders_area.push_back(
      border_feature->GetGeometryRef()->toPolygon()->get_Area());
  spdlog::info("Updating the border layer - done");

  OGRFeatureUniquePtr composite_table_feature(OGRFeature::CreateFeature(
      composite_table_layer->GetLayerDefn()));
  composite_table_feature->SetField(0, path.c_str());
  if (color_balancing)
    composite_table_feature->SetField(1, color_balancing_idx);
  composite_table_feature->SetGeometryDirectly(source_border.release());
  composite_table_layer->CreateFeature(composite_table_feature.get());
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
    int low_overview_trunc,
    int high_overview_trunc,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    int color_balancing_idx1,
    int color_balancing_idx2) {
  std::string string(
      "Running a mosaicking task for the raster pair with\n"
      " - Raster1 path: {}\n"
      " - Raster2 path: {}\n"
      " - Low overview trunction: {}\n"
      " - High overview trunction: {}\n"
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
  spdlog::info(string,path1, path2, low_overview_trunc, high_overview_trunc);
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
    spdlog::error("The geometry1 does not have the spatial reference");
    return nullptr;
  }
  if (!geometry2->getSpatialReference()) {
    spdlog::error("The geometry2 does not have the spatial reference");
    return nullptr;
  }
  if (!geometry1->getSpatialReference()->IsSame(
          geometry2->getSpatialReference())) {
    spdlog::error(
        "The geometry1 and the geometry2 do not have "
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
  if (low_overview_trunc < 0 || high_overview_trunc < 0) {
    spdlog::error(
        "The low overview trunction {} and the high overview trunction {} "
        "must be non-negative", low_overview_trunc, high_overview_trunc);
    return nullptr;
  }
  if ((low_overview_trunc + high_overview_trunc) >= overviews_count) {
    spdlog::warn(
        "The low overview trunction {} add the high overview trunction {} "
        "must be less than the overviews' count {}", low_overview_trunc,
        high_overview_trunc, overviews_count);
    high_overview_trunc = std::min(high_overview_trunc, overviews_count - 1);
    low_overview_trunc = overviews_count - 1 - high_overview_trunc;
    spdlog::info(
        "Reset arguments:\n"
        " - Low overview trunction: {}\n"
        " - High overview trunction: {}", low_overview_trunc,
        high_overview_trunc);
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
  auto factor(pow(2, overviews_count - 1 - low_overview_trunc));
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

  GDALDatasetUniquePtr
      overlap1_dataset(nullptr),
      overlap2_dataset(nullptr),
      label_raster_dataset(nullptr);
  for (int idx(overviews_count - 1 - low_overview_trunc);
      idx >= high_overview_trunc; idx--) {
    factor = pow(2, idx);
    spdlog::debug("Operating on the {} times downsampled overview", factor);
    CreateOverlapDatasets(
        factor, geotrans, rgb_bands_map, spatial_ref.get(), dataset1.get(),
        dataset2.get(), overlap1_geometry.get(), overlap2_geometry.get(),
        color_balancing, color_balancing_idx1, color_balancing_idx2,
        overlap1_dataset, overlap2_dataset, label_raster_dataset);
    UpdateMediums(
        idx, overlap1_dataset.get(), overlap2_dataset.get(),
        valid_geometry.get(), label_raster_dataset, overlap1_geometry,
        overlap2_geometry, idx == high_overview_trunc, 0.0);
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

void MosaickingBase::CreateOverlapDatasets(
    double factor,
    double* geotrans,
    const std::vector<int>& rgb_bands_map,
    GDALDataset* source_raster_dataset,
    OGRLayer* composite_table_layer,
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
  auto spatial_ref(composite_table_layer->GetSpatialRef());

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
  std::vector<GDALDataset*> original_rasters_dataset;
  std::vector<OGRGeometry*>
      color_balancing_intersecting_geometries,
      original_intersecting_geometries;
  for (const auto& feature : composite_table_layer) {
    if (feature->GetGeometryRef()->Intersect(covered_overlap_geometry)) {
      auto geoemtry(feature->GetGeometryRef()->Intersection(
          covered_overlap_geometry));
      geoemtry->assignSpatialReference(spatial_ref);
      if (color_balancing && feature->GetFieldAsInteger(1) != -1) {
        color_balanced_rasters_idx.push_back(feature->GetFieldAsInteger(1));
        color_balancing_intersecting_geometries.push_back(geoemtry);
      } else {
        original_rasters_dataset.push_back(GDALDataset::Open(
            feature->GetFieldAsString(0), GDAL_OF_RASTER | GDAL_OF_READONLY));
        original_intersecting_geometries.push_back(geoemtry);
      }
    }
  }
  if (!color_balanced_rasters_idx.empty()) {
    color_balancing->WarpByGeometry(
        color_balanced_rasters_idx, color_balancing_intersecting_geometries,
        covered_overlap_dataset, rgb_bands_map);
  }
  if (!original_rasters_dataset.empty()) {
    utils::WarpByGeometry(
        original_rasters_dataset, original_intersecting_geometries,
        covered_overlap_dataset, rgb_bands_map);
  }
  for (auto& geoemtry : color_balancing_intersecting_geometries)
    OGRGeometryFactory::destroyGeometry(geoemtry);
  for (int i(0); i < original_rasters_dataset.size(); ++i) {
    GDALClose(original_rasters_dataset[i]);
    OGRGeometryFactory::destroyGeometry(original_intersecting_geometries[i]);
  }
  spdlog::debug("Updating the covered overlap dataset - done");

  // Update the new overlap dataset
  new_overlap_geometry->assignSpatialReference(spatial_ref);
  if (color_balancing && color_balancing_idx != -1) {
    color_balancing->WarpByGeometry(
        {color_balancing_idx}, {new_overlap_geometry}, new_overlap_dataset,
        rgb_bands_map);
  } else {
    utils::WarpByGeometry(
        {source_raster_dataset}, {new_overlap_geometry}, new_overlap_dataset,
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
    GDALDataset* raster_dataset1,
    GDALDataset* raster_dataset2,
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
        {raster_dataset1}, {overlap_geometry1}, overlap_dataset1,
        rgb_bands_map);
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
        {raster_dataset2}, {overlap_geometry2}, overlap_dataset2,
        rgb_bands_map);
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

void MosaickingBase::UpdateMediums(
    int idx,
    GDALDataset* covered_overlap_dataset,
    GDALDataset* new_overlap_dataset,
    OGRGeometry* valid_geometry,
    GDALDatasetUniquePtr& label_raster_dataset,
    OGRGeometryUniquePtr& covered_overlap_geometry,
    OGRGeometryUniquePtr& new_overlap_geometry,
    bool last,
    double buffer_at_end) {
  // Execute the mosaicking algorithm on the overlap datasets
  double geotrans[6];
  new_overlap_dataset->GetGeoTransform(geotrans);
  cv::Mat covered_mat, new_mat, label_mat;
  PrepareData(
      covered_overlap_dataset, new_overlap_dataset, label_raster_dataset.get(),
      covered_mat, new_mat, label_mat);
  ExecuteMosaicking(
      covered_mat, new_mat, label_mat, geotrans,
      const_cast<OGRSpatialReference*>(new_overlap_dataset->GetSpatialRef()),
      label_raster_dataset, covered_overlap_geometry, new_overlap_geometry);
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
      } else {
        utils::RearrangeNeighborGeometries(
            seamline.get(), covered_overlap_geometry, new_overlap_geometry);
        utils::GraftSeamline(
            seamline, covered_overlap_geometry, new_overlap_geometry,
            tol_ * geotrans[1]);
      }
    }
  } else {
    spdlog::debug("Updating the overlap geometries");

    // Create the seamline polygon and the buffered seamline polygon
    seamline.reset(seamline->Simplify(5 * geotrans[1]));
    seamline.reset(seamline->Buffer(
        (buffers_[idx] - buffer_at_end) * geotrans[1]));
    seamline.reset(seamline->Simplify(5 * geotrans[1]));
    seamline.reset(seamline->Intersection(valid_geometry));
    OGRGeometryUniquePtr buffered_overlap_geometry(seamline->Buffer(
        (buffers_[idx] - buffer_at_end) * geotrans[1]));

    // Update overlap geometries
    covered_overlap_geometry.reset(covered_overlap_geometry->Union(
        seamline.get()));
    covered_overlap_geometry.reset(covered_overlap_geometry->Intersection(
        buffered_overlap_geometry.get()));
    covered_overlap_geometry.reset(covered_overlap_geometry->Buffer(
        buffer_at_end * geotrans[1]));
    new_overlap_geometry.reset(new_overlap_geometry->Union(seamline.get()));
    new_overlap_geometry.reset(new_overlap_geometry->Intersection(
        buffered_overlap_geometry.get()));
    new_overlap_geometry.reset(new_overlap_geometry->Buffer(
        buffer_at_end * geotrans[1]));
    spdlog::debug("Updating the overlap geometries - done");
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
  for (int i(0); i < composite_table_layer->GetFeatureCount(); ++i) {
    if (OGRFeatureUniquePtr feature(composite_table_layer->GetFeature(i));
        covered_geometry->Intersect(feature->GetGeometryRef())) {
      OGRGeometryUniquePtr diff_geometry(feature->GetGeometryRef()->Difference(
          source_border.get()));
      switch (diff_geometry->getGeometryType()) {
        case wkbPolygon: {
          dislocated_infos.emplace_back(i, std::move(diff_geometry));
          break;
        }
        case wkbMultiPolygon: {
          bool b(true);
          double max_area(0.0);
          for (const auto& geometry : diff_geometry->toMultiPolygon()) {
            dislocated_infos.emplace_back(i, geometry->clone());
            if (auto area(geometry->get_Area()); !b && area < max_area) {
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
      composite_table_layer->CreateFeature(composite_table_feature.get());
      border_layer->CreateFeature(border_feature.get());
      borders_area.push_back(area);
    }
  }
  spdlog::info("Updating the composite table - done");
}

GDALDatasetUniquePtr CreateMosaickingRaster(
    const std::string& path,
    OGRLayer* layer,
    double reso,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    double blend_dist) {
  auto driver(utils::GetRasterDriverByPath(path));
  if (!driver) {
    spdlog::error("Finding the driver for {} failed", path);
    return nullptr;
  }
  OGREnvelope enve;
  layer->GetExtent(&enve);
  auto x_size(static_cast<int>(ceil((enve.MaxX - enve.MinX) / reso))),
      y_size(static_cast<int>(ceil((enve.MaxY - enve.MinY) / reso)));
  spdlog::info(
      "Creating a {}x{} mosaicking raster to {}", x_size, y_size, path);
  GDALDatasetUniquePtr output_dataset(driver->Create(
      path.c_str(), x_size, y_size, 3, GDT_Byte, nullptr));
  if (!output_dataset) {
    spdlog::error("Creating {} failed", path);
    return nullptr;
  }

  double geotrans[6]{enve.MinX, reso, 0.0, enve.MaxY, 0.0, -reso};
  output_dataset->SetGeoTransform(geotrans);
  output_dataset->SetSpatialRef(layer->GetSpatialRef());

  // Warp rasters by the corresponding geometry to the mosaicking raster
  for (int i(0); i < layer->GetFeatureCount(); ++i) {
    OGRFeatureUniquePtr feature(layer->GetFeature(i));
    std::string source_path(feature->GetFieldAsString(0));
    spdlog::info(
        "Warping {} by the corresponding geometry to the mosaicking raster",
        source_path);
    GDALDatasetUniquePtr source_dataset(GDALDataset::Open(
        source_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    OGRGeometryUniquePtr geometry(feature->StealGeometry());
    geometry->assignSpatialReference(layer->GetSpatialRef());
    if (color_balancing && feature->GetFieldAsInteger(1) != -1) {
      color_balancing->WarpByGeometry(
          {feature->GetFieldAsInteger(1)}, {geometry.get()}, output_dataset, {},
          GRA_Bilinear, blend_dist);
    } else {
      utils::WarpByGeometry(
          {source_dataset.get()}, {geometry.get()}, output_dataset, {},
          GRA_Bilinear, blend_dist);
    }
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1, layer->GetFeatureCount());
  }
  spdlog::info("Creating a mosaicking raster - done");
  return output_dataset;
}

}  // namespace mosaicking
}  // namespace rs_toolset