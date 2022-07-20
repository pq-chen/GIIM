#include "mosaicking_base.h"

#include <cmath>

#include <algorithm>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gdal_alg.h>
#include <gdal_priv.h>
#include <nlohmann/json.hpp>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace mosaicking {

MosaickingBase::MosaickingBase(double tol)
    : mem_driver_(GetGDALDriverManager()->GetDriverByName("MEM")),
      memory_driver_(GetGDALDriverManager()->GetDriverByName("Memory")),
      tol_(tol) {
  buffers_per_unit_ = new int[8]{ 0, 10, 15, 20, 25, 30, 30, 30 };
  spdlog::info(
      "Creating a mosaicking base with\n"
      " - Tolerance for simplifying the seamline: {}", tol);
}

bool MosaickingBase::Run(
    const std::string& raster_path,
    OGRLayer* composite_table_layer,
    OGRLayer* border_layer,
    OGRGeometry* covered_border,
    int last_overview_idx,
    bool use_seamline) {
  spdlog::info(
      "Running a mosaicking task from {}\n"
      " - Last overview index: {}\n"
      " - Use seamline: {}", raster_path, last_overview_idx, use_seamline);
  GDALDatasetUniquePtr source_raster_dataset(GDALDataset::Open(
      raster_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!source_raster_dataset) {
    spdlog::warn("Opening {} failed", raster_path);
    return false;
  }
  double geotrans[6];
  if (source_raster_dataset->GetGeoTransform(geotrans) != CE_None ||
      !source_raster_dataset->GetSpatialRef()) { 
    spdlog::warn(
        "{} does not have the geotransform or the spatial reference. "
        "Please check whether the raster is DOM", raster_path);
    return false;
  }
  if (!composite_table_layer || !covered_border) {
    spdlog::warn(
        "The \"composite_table_layer\" or the \"covered_border\" argument is "
        "empty");
    return false;
  }
  utils::CreateRasterPyra(source_raster_dataset.get());
  int overviews_count(
      source_raster_dataset->GetRasterBand(1)->GetOverviewCount());
  if (last_overview_idx <= 0 || last_overview_idx > overviews_count) {
    spdlog::warn(
        "The \"last_overview_idx\" argument must be positive and less equal "
        "than overviews' count");
    return false;
  }

  // Create a border for the source raster and
  int ratio(4);
  OGRGeometryUniquePtr source_border(utils::CreateBorderGeometry(
      source_raster_dataset.get(), ratio));
  source_border.reset(source_border->Simplify(geotrans[1] * ratio * 1.5));
  source_border.reset(source_border->Buffer(geotrans[1] * ratio * -2));
  source_border->transformTo(composite_table_layer->GetSpatialRef());
  OGRGeometryUniquePtr
      origin_source_border(source_border->clone()),
      new_covered_polygon(source_border->clone());
  if (covered_border->IsEmpty() || !source_border->Intersect(covered_border)) {
    // Skip for the no intersection situation
    spdlog::info(
        "Skipping updating the source border and the composite table layer "
        "for the no intersection situation");
    spdlog::debug("Updating the covered border");
    covered_border->toMultiPolygon()->addGeometryDirectly(
        source_border->clone());
    spdlog::info("Updating the covered border - done");
  } else {
    auto covered_multi_polygon(covered_border->toMultiPolygon());

    // Traverse all polygons in the covered border
    for (int i = covered_multi_polygon->getNumGeometries() - 1; i >= 0; i--) {
      auto covered_polygon(covered_multi_polygon->getGeometryRef(i));

      // Continue the loop for the no intersection situation
      if (!source_border->Intersect(covered_polygon)) continue;

      // Exit for the covered situation 
      if (covered_polygon->Contains(source_border.get())) {
        spdlog::info(
            "Skipping adding the task "
            "since the covered border contains the source border");
        return true;
      }

      // Add a subtask
      spdlog::info("Adding a subtask with an intersected covered poylgon");
      if (use_seamline) {
        spdlog::debug("Updating the source border with the seamline");

        // Update the geotransform after reprojection
        char* spatial_ref_wkt(nullptr);
        composite_table_layer->GetSpatialRef()->exportToWkt(&spatial_ref_wkt);
        std::unique_ptr<char[]> _spatial_ref_wkt(spatial_ref_wkt);
        std::unique_ptr<void, void(*)(void*)> trans_arg(
            GDALCreateGenImgProjTransformer3(
                source_raster_dataset->GetProjectionRef(), geotrans,
                spatial_ref_wkt, nullptr),
            [](void* p) { GDALDestroyGenImgProjTransformer(p); });
        int x_size, y_size;
        GDALSuggestedWarpOutput(
            source_raster_dataset.get(), GDALGenImgProjTransform,
            trans_arg.get(), geotrans, &x_size, &y_size);

        GDALDatasetUniquePtr
            covered_overlap_dataset(nullptr),
            new_overlap_dataset(nullptr),
            label_raster_dataset(nullptr);
        OGRGeometryUniquePtr
            valid_geometry(source_border->Intersection(covered_polygon)),
            covered_overlap_geometry(nullptr),
            new_overlap_geometry(nullptr);
        valid_geometry.reset(valid_geometry->Buffer(geotrans[1] * ratio));
        int overviews_count(
            source_raster_dataset->GetRasterBand(1)->GetOverviewCount());
        for (int i = overviews_count - last_overview_idx; i >= 0; i--) {
          double downsample_factor(pow(2, i));
          spdlog::debug(
              "Operating on the {} times downsampled overview",
              downsample_factor);

          // Create overlap geometries for the first loop
          if (!covered_overlap_geometry) {
            spdlog::debug(
                "Initializing overlap geometries for the first loop");
            if (source_border->Contains(covered_polygon)) {
              spdlog::debug(
                  "Initializing overlap geometries in the surrounded case");
              covered_overlap_geometry.reset(covered_polygon->clone());
              new_overlap_geometry.reset(covered_overlap_geometry->Buffer(
                  -5 * downsample_factor * pow(2, last_overview_idx - 1) *
                  geotrans[1]));
              new_overlap_geometry.reset(source_border->Difference(
                  new_overlap_geometry.get()));
            } else {
              spdlog::info(
                  "Initializing overlap geometries in the overlap case");
              covered_overlap_geometry.reset(
                  source_border->Intersection(covered_polygon));
              covered_overlap_geometry.reset(covered_overlap_geometry->Buffer(
                  5 * downsample_factor * pow(2, last_overview_idx - 1) *
                  geotrans[1]));
              new_overlap_geometry.reset(covered_overlap_geometry->clone());
            }
            spdlog::info(
                "Initializing overlap geometries for the first loop - done");
          }
          UpdateOverlapDatasets(
              downsample_factor, geotrans, composite_table_layer,
              source_raster_dataset.get(), covered_overlap_geometry.get(),
              new_overlap_geometry.get(), covered_overlap_dataset,
              new_overlap_dataset, label_raster_dataset);
          UpdateMediums(
              i, covered_overlap_dataset.get(), new_overlap_dataset.get(),
              valid_geometry.get(), label_raster_dataset,
              covered_overlap_geometry, new_overlap_geometry);
          spdlog::info(
              "Operating on the {} times downsampled overview - done",
              downsample_factor);
        }
        UpdateResults(
            covered_overlap_geometry.get(), new_overlap_geometry.get(),
            covered_polygon, border_layer, source_border,
            composite_table_layer);
        spdlog::debug("Updating the source border with the seamline - done");
      } else {
        spdlog::debug("Updating the source border without the seamline");
        for (const auto& feature : composite_table_layer) {
          OGRGeometry* geometry(feature->GetGeometryRef());
          if (source_border->Intersect(geometry))
            source_border.reset(source_border->Difference(geometry));
        }
        spdlog::debug(
            "Updating the source border without the seamline - done");
      }

      // Delete the processed covered polygon and union the new covered polygon
      new_covered_polygon.reset(new_covered_polygon->Union(covered_polygon));
      covered_multi_polygon->removeGeometry(i);
    }

    // Add the new covered polygon
    spdlog::debug("Updating the covered border");
    covered_multi_polygon->addGeometryDirectly(new_covered_polygon.release());
    spdlog::info("Updating the covered border - done");
    spdlog::info(
        "Adding a subtask with an intersected covered poylgon - done");
  }

  // Update the border layer
  OGRFeatureUniquePtr border_feature(OGRFeature::CreateFeature(
      border_layer->GetLayerDefn()));
  border_feature->SetGeometryDirectly(origin_source_border.release());
  border_layer->CreateFeature(border_feature.get());

  // Create a feature for the source border
  spdlog::debug("Create a feature for the source border");
  OGRFeatureUniquePtr composite_table_feature(OGRFeature::CreateFeature(
      composite_table_layer->GetLayerDefn()));
  composite_table_feature->SetField(0, raster_path.c_str());
  composite_table_feature->SetGeometryDirectly(source_border.release());
  composite_table_layer->CreateFeature(composite_table_feature.get());
  spdlog::debug("Create a feature for the source border - done");
  spdlog::info("Running the mosaicking task - done");
  return true;
}

void MosaickingBase::UpdateOverlapDatasets(
    double downsample_factor,
    double* geotrans,
    OGRLayer* composite_table_layer,
    GDALDataset* source_raster_dataset,
    OGRGeometry* covered_overlap_geometry,
    OGRGeometry* new_overlap_geometry,
    GDALDatasetUniquePtr& covered_overlap_dataset,
    GDALDatasetUniquePtr& new_overlap_dataset,
    GDALDatasetUniquePtr& label_raster_dataset) {
  spdlog::debug("Updating overlap datasets");

  // Create the union range and geotransform
  OGRGeometryUniquePtr union_geometry(
      covered_overlap_geometry->Union(new_overlap_geometry));
  OGREnvelope enve;
  union_geometry->getEnvelope(&enve);
  double inv_geotrans[6], union_geotrans[6], x, y;
  GDALInvGeoTransform(geotrans, inv_geotrans);
  int range[4];
  GDALApplyGeoTransform(inv_geotrans, enve.MinX, enve.MaxY, &x, &y);
  range[0] = static_cast<int>(floor(x));
  range[1] = static_cast<int>(floor(y));
  GDALApplyGeoTransform(inv_geotrans, enve.MaxX, enve.MinY, &x, &y);
  range[2] = static_cast<int>(ceil((x - range[0]) / downsample_factor));
  range[3] = static_cast<int>(ceil((y - range[1]) / downsample_factor));
  GDALApplyGeoTransform(
      geotrans, range[0], range[1], union_geotrans, union_geotrans + 3);
  union_geotrans[1] = downsample_factor * geotrans[1];
  union_geotrans[2] = 0;
  union_geotrans[4] = 0;
  union_geotrans[5] = downsample_factor * geotrans[5];
  auto spatial_ref(composite_table_layer->GetSpatialRef());
  spdlog::debug("Creating the union range and geotransform - done");

  // Update the covered overlap dataset
  covered_overlap_dataset.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  covered_overlap_dataset->SetGeoTransform(union_geotrans);
  covered_overlap_dataset->SetSpatialRef(spatial_ref);
  std::vector<GDALDataset*> rasters_dataset;
  std::vector<OGRGeometry*> intersected_geometries;
  for (const auto& feature : composite_table_layer)
    if (feature->GetGeometryRef()->Intersect(covered_overlap_geometry)) {
      rasters_dataset.push_back(GDALDataset::Open(
          feature->GetFieldAsString(0), GDAL_OF_RASTER | GDAL_OF_READONLY));
      intersected_geometries.push_back(
          feature->GetGeometryRef()->Intersection(covered_overlap_geometry));
    }
  utils::WarpByGeometry(
      rasters_dataset, intersected_geometries, covered_overlap_dataset);
  for (int i = 0; i < rasters_dataset.size(); i++) {
    GDALClose(rasters_dataset[i]);
    OGRGeometryFactory::destroyGeometry(intersected_geometries[i]);
  }
  spdlog::debug("Updating the covered overlap dataset - done");

  // Update the new overlap dataset
  new_overlap_dataset.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  new_overlap_dataset->SetGeoTransform(union_geotrans);
  new_overlap_dataset->SetSpatialRef(spatial_ref);
  utils::WarpByGeometry(
      { source_raster_dataset }, { new_overlap_geometry },
      new_overlap_dataset);
  spdlog::debug("Updating the new overlap dataset - done");

  // Update the label raster dataset if necessary
  if (label_raster_dataset) {
    GDALDatasetUniquePtr _label_raster_dataset(mem_driver_->Create(
        "", range[2], range[3], 1, GDT_Byte, nullptr));
    _label_raster_dataset->SetGeoTransform(union_geotrans);
    _label_raster_dataset->SetSpatialRef(spatial_ref);
    utils::WarpByGeometry(
        { label_raster_dataset.get() }, { nullptr }, _label_raster_dataset,
        GRA_NearestNeighbour);
    label_raster_dataset.swap(_label_raster_dataset);
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
    OGRGeometryUniquePtr& new_overlap_geometry) {
  // Execute the mosaicking algorithm on the overlap datasets
  double geotrans[6];
  new_overlap_dataset->GetGeoTransform(geotrans);
  cv::Mat covered_mat, new_mat, label_mat;
  PrepareData(
      covered_overlap_dataset, new_overlap_dataset, label_raster_dataset.get(),
      covered_mat, new_mat, label_mat);
  ExecuteMosaicking(
      covered_mat, new_mat, label_mat, geotrans,
      const_cast<OGRSpatialReference*>(
          covered_overlap_dataset->GetSpatialRef()),
      label_raster_dataset, covered_overlap_geometry, new_overlap_geometry);

  // Create the seamline line string
  OGRGeometryUniquePtr seamline_geometry(
      covered_overlap_geometry->Intersection(new_overlap_geometry.get()));
  std::vector<std::pair<double, double>> seamline_coors;
  CreateSeamlineLineString(seamline_geometry, seamline_coors);
  covered_overlap_geometry.reset(OGRGeometryFactory::forceToPolygon(
      covered_overlap_geometry.release()));
  new_overlap_geometry.reset(OGRGeometryFactory::forceToPolygon(
      new_overlap_geometry.release()));

  if (idx == 0) {
    if (tol_) {
      // Rearrange points in label geometries
      RearrangeLabelGeometries(
          seamline_coors, covered_overlap_geometry, new_overlap_geometry);

      // Simplify label geometries on the last operation
      SimplifyLabelGeometries(
          tol_ * geotrans[1], seamline_geometry, covered_overlap_geometry,
          new_overlap_geometry);
    }
  } else {
    spdlog::debug("Updating the overlap geometries");

    // Create the seamline polygon and the buffered seamline polygon
    seamline_geometry.reset(seamline_geometry->Simplify(5 * geotrans[1]));
    seamline_geometry.reset(seamline_geometry->Buffer(
        (buffers_per_unit_[idx] - 2) * geotrans[1]));
    seamline_geometry.reset(seamline_geometry->Intersection(
        valid_geometry));
    OGRGeometryUniquePtr buffered_overlap_geometry(
        seamline_geometry->Buffer((buffers_per_unit_[idx] - 2) * geotrans[1]));

    // Update overlap geometries
    covered_overlap_geometry.reset(
        covered_overlap_geometry->Union(seamline_geometry.get()));
    covered_overlap_geometry.reset(covered_overlap_geometry->Intersection(
        buffered_overlap_geometry.get()));
    covered_overlap_geometry.reset(
        covered_overlap_geometry->Buffer(2 * geotrans[1]));
    new_overlap_geometry.reset(
        new_overlap_geometry->Union(seamline_geometry.get()));
    new_overlap_geometry.reset(new_overlap_geometry->Intersection(
        buffered_overlap_geometry.get()));
    new_overlap_geometry.reset(new_overlap_geometry->Buffer(2 * geotrans[1]));
    spdlog::debug("Updating the overlap geometries - done");
  }
}

void MosaickingBase::CreateSeamlineLineString(
    OGRGeometryUniquePtr& seamline_geometry,
    std::vector<std::pair<double, double>>& seamline_coors) {
  seamline_coors.resize(0);
  auto seamline_geometry_json(nlohmann::json::parse(
      seamline_geometry->exportToJson()));
  auto it(seamline_geometry_json.at("coordinates").begin());
  std::vector<std::pair<double, double>> coors{ *it->begin() };
  while (it != seamline_geometry_json.at("coordinates").end() &&
      (it->rbegin()->at(0) == coors.back().first ||
      it->rbegin()->at(1) == coors.back().second)) {
    coors.push_back(*it->rbegin());
    it++;
  }
  while (it != seamline_geometry_json.at("coordinates").end()) {
    seamline_coors.push_back(*it->begin());
    it++;
  }
  seamline_coors.insert(seamline_coors.end(), coors.begin(), coors.end());
  nlohmann::json seamline_line_string_json;
  seamline_line_string_json["type"] = "LineString";
  seamline_line_string_json["coordinates"] = seamline_coors;
  seamline_geometry.reset(OGRGeometryFactory::createFromGeoJson(
      seamline_line_string_json.dump().c_str()));
}

void MosaickingBase::RearrangeLabelGeometries(
    const std::vector<std::pair<double, double>>& seamline_coors,
    OGRGeometryUniquePtr& label0_geometry,
    OGRGeometryUniquePtr& label1_geometry) {
  auto label0_geometry_json(nlohmann::json::parse(
          label0_geometry->exportToJson())),
      label1_geometry_json(nlohmann::json::parse(
          label1_geometry->exportToJson()));
  if (std::find(
          seamline_coors.begin(), seamline_coors.end(),
          label0_geometry_json.at("coordinates")[0][0]
          .get<std::pair<double, double>>()) != seamline_coors.end()) {
    std::vector<std::pair<double, double>> old_coors, new_coors;
    label0_geometry_json.at("coordinates")[0].get_to(old_coors);
    int it1_trunc(-1), it2_trunc(-1);
    decltype(old_coors.begin()) it1, it2;
    do {
      it1_trunc++;
      it1 = std::find(
          old_coors.begin(), old_coors.end(), 
          *(seamline_coors.begin() + it1_trunc));
    } while (it1 == old_coors.end());
    do {
      it2_trunc++;
      it2 = std::find(
          old_coors.begin(), old_coors.end(),
          *(seamline_coors.rbegin() + it2_trunc));
    } while (it2 == old_coors.end());
    if (auto it(std::find(it1 + 1, old_coors.end(), *it1));
        it != old_coors.end() &&
        it - it2 + 1 + it1_trunc + it2_trunc != seamline_coors.size()) {
        it1 = it;
    } else if (auto it(std::find(it2 + 1, old_coors.end(), *it2));
        it != old_coors.end() &&
        it - it1 + 1 + it1_trunc + it2_trunc != seamline_coors.size()) {
        it2 = it;
    }
    if (it2 > it1) {
      new_coors.insert(new_coors.end(), it1 + (it2 - it1) / 2, it2);
      new_coors.insert(
          new_coors.end(), seamline_coors.rbegin(), seamline_coors.rend());
      new_coors.insert(new_coors.end(), it1 + 1, it1 + (it2 - it1) / 2 + 1);
    } else {
      new_coors.insert(new_coors.end(), it2 + (it1 - it2) / 2, it1);
      new_coors.insert(
          new_coors.end(), seamline_coors.begin(), seamline_coors.end());
      new_coors.insert(new_coors.end(), it2 + 1, it2 + (it1 - it2) / 2 + 1);
    }
    label0_geometry_json.at("coordinates")[0] = new_coors;
    label0_geometry.reset(OGRGeometryFactory::createFromGeoJson(
        label0_geometry_json.dump().c_str()));
  }
  if (std::find(
          seamline_coors.begin(), seamline_coors.end(),
          label1_geometry_json.at("coordinates")[0][0]
          .get<std::pair<double, double>>()) != seamline_coors.end()) {
    std::vector<std::pair<double, double>> old_coors, new_coors;
    label1_geometry_json.at("coordinates")[0].get_to(old_coors);
    int it1_trunc(-1), it2_trunc(-1);
    decltype(old_coors.begin()) it1, it2;
    do {
      it1_trunc++;
      it1 = std::find(
          old_coors.begin(), old_coors.end(), 
          *(seamline_coors.begin() + it1_trunc));
    } while (it1 == old_coors.end());
    do {
      it2_trunc++;
      it2 = std::find(
          old_coors.begin(), old_coors.end(),
          *(seamline_coors.rbegin() + it2_trunc));
    } while (it2 == old_coors.end());
    if (auto it(std::find(it1 + 1, old_coors.end(), *it1));
        it != old_coors.end() &&
        it - it2 + 1 + it1_trunc + it2_trunc == seamline_coors.size()) {
        it1 = it;
    } else if (auto it(std::find(it2 + 1, old_coors.end(), *it2));
        it != old_coors.end() &&
        it - it1 + 1 + it1_trunc + it2_trunc == seamline_coors.size()) {
        it2 = it;
    }
    if (it2 > it1) {
      new_coors.insert(new_coors.end(), it1 + (it2 - it1) / 2, it2);
      new_coors.insert(
          new_coors.end(), seamline_coors.rbegin(), seamline_coors.rend());
      new_coors.insert(new_coors.end(), it1 + 1, it1 + (it2 - it1) / 2 + 1);
    } else {
      new_coors.insert(new_coors.end(), it2 + (it1 - it2) / 2, it1);
      new_coors.insert(
          new_coors.end(), seamline_coors.begin(), seamline_coors.end());
      new_coors.insert(new_coors.end(), it2 + 1, it2 + (it1 - it2) / 2 + 1);
    }
    label1_geometry_json.at("coordinates")[0] = new_coors;
    label1_geometry.reset(OGRGeometryFactory::createFromGeoJson(
        label1_geometry_json.dump().c_str()));
  }
}

void MosaickingBase::SimplifyLabelGeometries(
    double tol,
    OGRGeometryUniquePtr& seamline_line_string,
    OGRGeometryUniquePtr& label0_geometry,
    OGRGeometryUniquePtr& label1_geometry) {
  seamline_line_string.reset(seamline_line_string->Simplify(tol));
  auto label0_geometry_json(nlohmann::json::parse(
          label0_geometry->exportToJson())),
      label1_geometry_json(nlohmann::json::parse(
          label1_geometry->exportToJson())),
      seamline_line_string_json(nlohmann::json::parse(
          seamline_line_string->exportToJson()));

  // Simplify the label0 geometry
  spdlog::debug("Simplify the label0 geometry");
  std::list<std::pair<double, double>> seamline_coors;
  std::vector<std::pair<double, double>> label_coors, label_new_coors;
  seamline_line_string_json.at("coordinates").get_to(seamline_coors);
  label0_geometry_json.at("coordinates")[0].get_to(label_coors);
  auto it1(std::find(
          label_coors.begin(), label_coors.end(), seamline_coors.front())),
      it2(std::find(
          label_coors.begin(), label_coors.end(), seamline_coors.back()));
  while (it1 == label_coors.end()) {
    seamline_coors.pop_front();
    it1 = std::find(
        label_coors.begin(), label_coors.end(), seamline_coors.front());
  }
  while (it2 == label_coors.end()) {
    seamline_coors.pop_back();
    it2 = std::find(
        label_coors.begin(), label_coors.end(), seamline_coors.back());
  }
  if (it2 > it1) {
    label_new_coors.insert(label_new_coors.end(), label_coors.begin(), it1);
    label_new_coors.insert(
        label_new_coors.end(), seamline_coors.begin(), seamline_coors.end());
    label_new_coors.insert(label_new_coors.end(), it2 + 1, label_coors.end());
  } else {
    label_new_coors.insert(label_new_coors.end(), label_coors.begin(), it2);
    label_new_coors.insert(
        label_new_coors.end(), seamline_coors.rbegin(), seamline_coors.rend());
    label_new_coors.insert(label_new_coors.end(), it1 + 1, label_coors.end());
  }
  label0_geometry_json.at("coordinates")[0] = label_new_coors;
  label0_geometry.reset(OGRGeometryFactory::createFromGeoJson(
      label0_geometry_json.dump().c_str()));
  spdlog::info("Simplify the label0 geometry - done");

  // Simplify the label1 geometry
  spdlog::debug("Simplify the label1 geometry");
  label_new_coors.resize(0);
  label1_geometry_json.at("coordinates")[0].get_to(label_coors);
  it1 = std::find(
      label_coors.begin(), label_coors.end(), seamline_coors.front());
  it2 = std::find(
      label_coors.begin(), label_coors.end(), seamline_coors.back());
  while (it1 == label_coors.end()) {
    seamline_coors.pop_front();
    it1 = std::find(
        label_coors.begin(), label_coors.end(), seamline_coors.front());
  }
  while (it2 == label_coors.end()) {
    seamline_coors.pop_back();
    it2 = std::find(
        label_coors.begin(), label_coors.end(), seamline_coors.back());
  }
  if (it2 > it1) {
    label_new_coors.insert(label_new_coors.end(), label_coors.begin(), it1);
    label_new_coors.insert(
        label_new_coors.end(), seamline_coors.begin(), seamline_coors.end());
    label_new_coors.insert(label_new_coors.end(), it2 + 1, label_coors.end());
  } else {
    label_new_coors.insert(label_new_coors.end(), label_coors.begin(), it2);
    label_new_coors.insert(
        label_new_coors.end(), seamline_coors.rbegin(),seamline_coors.rend());
    label_new_coors.insert(label_new_coors.end(), it1 + 1, label_coors.end());
  }
  label1_geometry_json.at("coordinates")[0] = label_new_coors;
  label1_geometry.reset(OGRGeometryFactory::createFromGeoJson(
      label1_geometry_json.dump().c_str()));
  spdlog::info("Simplify the label1 geometry - done");
}

void MosaickingBase::UpdateResults(
    OGRGeometry* covered_geometry,
    OGRGeometry* new_geometry,
    OGRGeometry* covered_polygon,
    OGRLayer* border_layer,
    OGRGeometryUniquePtr& source_border,
    OGRLayer* composite_table_layer) {
  // Update the source border
  spdlog::debug("Updating the source border");
  OGRGeometryUniquePtr diff_geometry(
      source_border->Difference(covered_geometry));
  source_border.reset(source_border->Difference(covered_polygon));
  for (const auto& geometry : diff_geometry->toMultiPolygon())
    if (geometry->Intersect(new_geometry)) {
      source_border.reset(source_border->Union(geometry));
      break;
    }
  spdlog::info("Updating the source border - done");

  // Update the composite table
  spdlog::debug("Updating the composite table");
  std::vector<std::pair<int, OGRGeometryUniquePtr>> dislocated_geometries;
  for (auto& feature : composite_table_layer) {
    int id(static_cast<int>(feature->GetFID()));
    if (!covered_geometry->Intersect(feature->GetGeometryRef())) continue;
    diff_geometry.reset(
        feature->GetGeometryRef()->Difference(source_border.get()));
    if (diff_geometry->getGeometryType() == wkbPolygon) {
      dislocated_geometries.emplace_back(id, diff_geometry->clone());
    } else {
      bool be_head(true);
      double max_area(0.0);
      for (const auto& geometry : diff_geometry->toMultiPolygon()) {
        dislocated_geometries.emplace_back(id, geometry->clone());
        if (!be_head && geometry->get_Area() < max_area) {
          dislocated_geometries[dislocated_geometries.size() - 2].second.swap(
              dislocated_geometries.back().second);
        } else {
          max_area = geometry->get_Area();
        }
        be_head = false;
      }
    }
    feature->SetGeometryDirectly(
        dislocated_geometries.back().second.release());
    composite_table_layer->SetFeature(feature.get());
    dislocated_geometries.pop_back();
  }
  decltype(dislocated_geometries) new_geometries;
  for (auto& info : dislocated_geometries) {
    bool b(false);
    for (auto& feature : composite_table_layer)
      if (info.second->Intersects(feature->GetGeometryRef()) &&
          border_layer->GetFeature(feature->GetFID())->GetGeometryRef()
              ->Contains(info.second.get())) {
        feature->SetGeometryDirectly(
            info.second->Union(feature->GetGeometryRef()));
        composite_table_layer->SetFeature(feature.get());
        b = true;
        break;
      }
    if (!b)
      new_geometries.push_back(std::move(info));
  }
  for (auto& info : new_geometries) {
    OGRFeatureUniquePtr composite_table_feature(OGRFeature::CreateFeature(
        composite_table_layer->GetLayerDefn()));
    composite_table_feature->SetField(
        0, composite_table_layer->GetFeature(info.first)->GetFieldAsString(0));
    composite_table_feature->SetGeometryDirectly(info.second.release());
    composite_table_layer->CreateFeature(composite_table_feature.get());
    OGRFeatureUniquePtr border_feature(OGRFeature::CreateFeature(
        border_layer->GetLayerDefn()));
    border_feature->SetGeometryDirectly(
        border_layer->GetFeature(info.first)->GetGeometryRef()->clone());
    border_layer->CreateFeature(border_feature.get());
  }
  spdlog::info("Updating the composite table - done");
}

} // mosaicking
} // rs_toolset