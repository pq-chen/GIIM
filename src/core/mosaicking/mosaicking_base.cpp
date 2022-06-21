#include "mosaicking_base.h"

#include <string>
#include <utility>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace mosaicking {

MosaickingBase::MosaickingBase()
    : mem_driver_(GetGDALDriverManager()->GetDriverByName("MEM")),
      memory_driver_(GetGDALDriverManager()->GetDriverByName("Memory")) {}

bool MosaickingBase::Run(
    OGRLayer* composite_table_layer,
    OGRGeometry* covered_border,
    const std::string& raster_path,
    bool with_seamline) {
  spdlog::info("Adding a task from {}", raster_path);
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
  utils::CreateRasterPyra(source_raster_dataset.get(), "NONE");
  OGRGeometryUniquePtr 
      source_border(utils::CreateBorderGeometry(source_raster_dataset.get())),
      new_covered_polygon(source_border->clone());
  if (covered_border->IsEmpty() || !source_border->Intersect(covered_border)) {
    spdlog::info(
        "Skipping updating the source border and the composite table layer"
        "for no intersection situation");
    covered_border->toMultiPolygon()->addGeometryDirectly(
        source_border->clone());
    spdlog::debug("Updating the covered border - done");
  } else {
    OGRMultiPolygon* _covered_border_(covered_border->toMultiPolygon());
    for (int i = _covered_border_->getNumGeometries() - 1; i >= 0; i--) {
      OGRGeometry* covered_polygon(_covered_border_->getGeometryRef(i));
      if (!source_border->Intersect(covered_polygon)) continue;
      if (covered_polygon->Contains(source_border.get())) {
        spdlog::info(
            "Skipping adding a task "
            "since the covered border contains the source border");
        return true;
      }
      spdlog::info("Adding a subtask from the intersected border");
      if (with_seamline) {
        spdlog::debug("Creating the seamline vector");
        spdlog::debug("Initializing the seamline layer");
        GDALDatasetUniquePtr seamline_vector_dataset(
            memory_driver_->Create("", 0, 0, 0, GDT_Unknown, nullptr));
        OGRLayer* seamline_layer(seamline_vector_dataset->CreateLayer(
            "", const_cast<OGRSpatialReference*>(
                source_raster_dataset->GetSpatialRef()), wkbPolygon));
        OGRFieldDefn label_field("label", OFTInteger);
        seamline_layer->CreateField(&label_field);
        spdlog::debug("Initializing the seamline layer - done");
        GDALDatasetUniquePtr
            covered_overlap_dataset(nullptr),
            new_overlap_dataset(nullptr),
            label_raster_dataset(nullptr);
        OGRGeometryUniquePtr 
            covered_overlap_geometry(nullptr),
            new_overlap_geometry(nullptr);
        int overviews_count(
            source_raster_dataset->GetRasterBand(1)->GetOverviewCount());
        for (int i = overviews_count - 1; i >= 0; i--) {
          double downsample_factor(pow(2, i));
          spdlog::debug(
              "Operating on the {} times downsampled overview", 
              downsample_factor);
          if (!covered_overlap_geometry) {
            spdlog::debug(
                "Initializing overlap geometries for the last overview");
            if (source_border->Contains(covered_polygon)) {
              spdlog::debug(
                  "Initializing overlap geometries in the encirclement case");
              covered_overlap_geometry.reset(covered_polygon->clone());
              new_overlap_geometry.reset(covered_overlap_geometry->Buffer(
                  -5 * downsample_factor * geotrans[1]));
              new_overlap_geometry.reset(source_border->Difference(
                  new_overlap_geometry.get()));
            } else {
              covered_overlap_geometry.reset(
                  source_border->Intersection(covered_polygon));
              covered_overlap_geometry.reset(covered_overlap_geometry->Buffer(
                  5 * downsample_factor * geotrans[1]));
              new_overlap_geometry.reset(covered_overlap_geometry->clone());
            }
            spdlog::debug(
                "Initializing overlap geometries for the last overview - done");
          }
          CreateOverlapDatasets(
              composite_table_layer, source_raster_dataset.get(),
              covered_overlap_geometry.get(), new_overlap_geometry.get(),
              covered_overlap_dataset, new_overlap_dataset,
              label_raster_dataset, downsample_factor);
          UpdateOverlapGeometries(
              covered_overlap_dataset.get(), new_overlap_dataset.get(), 
              label_raster_dataset, covered_overlap_geometry,
              new_overlap_geometry, i ? nullptr : seamline_layer, i);
          spdlog::info(
              "Operating on the {} times downsampled overview - done",
              downsample_factor);
        }
        spdlog::debug("Creating the seamline vector - done");
        AddSeamline(composite_table_layer, seamline_layer, source_border);
      } else {
        spdlog::debug("Updating the source border without the seamline");
        for (const auto& composite_table_feature : composite_table_layer) {
          OGRGeometry* composite_table_geometry(
              composite_table_feature->GetGeometryRef());
          if (source_border->Intersect(composite_table_geometry))
            source_border.reset(
                source_border->Difference(composite_table_geometry));
        }
        spdlog::debug("Updating the source border without the seamline - done");
      }
      new_covered_polygon.reset(new_covered_polygon->Union(covered_polygon));
      _covered_border_->removeGeometry(i);
    }
    _covered_border_->addGeometryDirectly(new_covered_polygon.get());
    new_covered_polygon.release();
    spdlog::debug("Updating the covered border - done");
    spdlog::info("Adding a subtask from the intersected border - done");
  }
  OGRFeatureUniquePtr composite_table_feature(OGRFeature::CreateFeature(
      composite_table_layer->GetLayerDefn()));
  if (composite_table_layer->GetSpatialRef())
    source_border->transformTo(composite_table_layer->GetSpatialRef());
  composite_table_feature->SetField(0, raster_path.c_str());
  composite_table_feature->SetGeometryDirectly(source_border.get());
  composite_table_layer->CreateFeature(composite_table_feature.get());
  source_border.release();
  spdlog::debug("Creating the feature with the source border - done");
  spdlog::info("Adding a task from {} - done", raster_path);
  return true;
}

void MosaickingBase::CreateOverlapDatasets(
    OGRLayer* composite_table_layer,
    GDALDataset* source_raster_dataset,
    OGRGeometry* covered_overlap_geometry,
    OGRGeometry* new_overlap_geometry,
    GDALDatasetUniquePtr& covered_overlap_dataset,
    GDALDatasetUniquePtr& new_overlap_dataset,
    GDALDatasetUniquePtr& label_raster_dataset,
    double downsample_factor) {
  spdlog::debug("Creating overlap datasets");
  double source_geotrans[6], source_inv_geotrans[6], union_geotrans[6], x, y;
  source_raster_dataset->GetGeoTransform(source_geotrans);
  GDALInvGeoTransform(source_geotrans, source_inv_geotrans);
  OGRGeometryUniquePtr union_geometry(
      covered_overlap_geometry->Union(new_overlap_geometry));
  OGREnvelope enve;
  union_geometry->getEnvelope(&enve);
  int range[4];
  GDALApplyGeoTransform(source_inv_geotrans, enve.MinX, enve.MaxY, &x, &y);
  range[0] = static_cast<int>(floor(x));
  range[1] = static_cast<int>(floor(y));
  GDALApplyGeoTransform(source_inv_geotrans, enve.MaxX, enve.MinY, &x, &y);
  range[2] = static_cast<int>(ceil((x - range[0]) / downsample_factor));
  range[3] = static_cast<int>(ceil((y - range[1]) / downsample_factor));
  GDALApplyGeoTransform(
      source_geotrans, range[0], range[1], union_geotrans,
      union_geotrans + 3);
  union_geotrans[1] = downsample_factor * source_geotrans[1];
  union_geotrans[2] = 0;
  union_geotrans[4] = 0;
  union_geotrans[5] = downsample_factor * source_geotrans[5];
  spdlog::debug("Creating overlap datasets' range and geotransform - done");
  covered_overlap_dataset.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  covered_overlap_dataset->SetGeoTransform(union_geotrans);
  covered_overlap_dataset->SetSpatialRef(
      source_raster_dataset->GetSpatialRef());
  std::vector<GDALDataset*> rasters_dataset;
  for (const auto& composite_table_feature : composite_table_layer)
    if (composite_table_feature->GetGeometryRef()->Intersect(
        covered_overlap_geometry))
      rasters_dataset.push_back(GDALDataset::Open(
          composite_table_feature->GetFieldAsString(0), 
          GDAL_OF_RASTER | GDAL_OF_READONLY));
  utils::WarpByGeometry(
      rasters_dataset, covered_overlap_dataset.get(), covered_overlap_geometry);
  for (const auto& raster_dataset : rasters_dataset)
    GDALClose(raster_dataset);
  spdlog::debug("Creating the covered overlap dataset - done");
  new_overlap_dataset.reset(mem_driver_->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  new_overlap_dataset->SetGeoTransform(union_geotrans);
  new_overlap_dataset->SetSpatialRef(source_raster_dataset->GetSpatialRef());
  utils::WarpByGeometry(
      { source_raster_dataset }, new_overlap_dataset.get(),
      new_overlap_geometry);
  spdlog::debug("Creating the new overlap dataset - done");
  if (label_raster_dataset) {
    GDALDataset* _label_raster_dataset(mem_driver_->Create(
        "", range[2], range[3], 1, GDT_Byte, nullptr));
    _label_raster_dataset->SetGeoTransform(union_geotrans);
    _label_raster_dataset->SetSpatialRef(label_raster_dataset->GetSpatialRef());
    utils::WarpByGeometry(
        { label_raster_dataset.get() }, _label_raster_dataset, nullptr,
        GRA_NearestNeighbour);
    label_raster_dataset.reset(_label_raster_dataset);
    spdlog::debug("Creating the label raster dataset - done");
  }
  spdlog::debug("Creating overlap datasets - done");
}

void MosaickingBase::UpdateOverlapGeometries(
    GDALDataset* covered_raster_dataset,
    GDALDataset* new_raster_dataset,
    GDALDatasetUniquePtr& label_raster_dataset,
    OGRGeometryUniquePtr& covered_overlap_geometry,
    OGRGeometryUniquePtr& new_overlap_geometry,
    OGRLayer* seamline_layer,
    int overview_idx) {
  spdlog::debug("Updating overlap geometries");
  cv::Mat covered_mat, new_mat, label_mat;
  PrepareData(
      covered_raster_dataset, new_raster_dataset, label_raster_dataset.get(),
      covered_mat, new_mat, label_mat);
  OGRGeometryUniquePtr label0_geometry(nullptr), label1_geometry(nullptr);
  CreateSeamlineGeometries(
      new_raster_dataset, covered_mat, new_mat, label_mat, label_raster_dataset,
      label0_geometry, label1_geometry, seamline_layer);
  if (!seamline_layer) {
    covered_overlap_geometry.reset(OGRGeometryFactory::forceToPolygon(
        label0_geometry->Boundary()));
    new_overlap_geometry.reset(OGRGeometryFactory::forceToPolygon(
        label1_geometry->Boundary()));
    OGRGeometryUniquePtr overlap_geometry(
        covered_overlap_geometry->Intersection(new_overlap_geometry.get()));
    if (overlap_geometry->getGeometryType() == wkbMultiLineString) {
      nlohmann::json 
          overlap_geometry_json(nlohmann::json::parse(
              overlap_geometry->exportToJson())),
          overlap_line_string_json;
      overlap_line_string_json["type"] = "LineString";
      overlap_line_string_json["coordinates"] = {};
      for (const auto& line_string : overlap_geometry_json.at("coordinates"))
        overlap_line_string_json["coordinates"].push_back(*line_string.begin());
      overlap_geometry.reset(OGRGeometryFactory::createFromGeoJson(
          overlap_line_string_json.dump().c_str()));
    }
    double geotrans[6], buffer(fmax(fmin(30, 5 * (overview_idx + 1)), 10));
    new_raster_dataset->GetGeoTransform(geotrans);
    overlap_geometry.reset(overlap_geometry->Simplify(5 * geotrans[1]));
    overlap_geometry.reset(overlap_geometry->Buffer(buffer * geotrans[1]));
    OGRGeometryUniquePtr buffered_overlap_geometry(
        overlap_geometry->Buffer(buffer * geotrans[1]));
    covered_overlap_geometry.reset(
        label0_geometry->Union(overlap_geometry.get()));
    covered_overlap_geometry.reset(covered_overlap_geometry->Intersection(
        buffered_overlap_geometry.get()));
    new_overlap_geometry.reset(label1_geometry->Union(overlap_geometry.get()));
    new_overlap_geometry.reset(new_overlap_geometry->Intersection(
        buffered_overlap_geometry.get()));
  }
  spdlog::debug("Updating overlap geometries - done");
}

void MosaickingBase::AddSeamline(
    OGRLayer* composite_table_layer,
    OGRLayer* seamline_layer,
    OGRGeometryUniquePtr& source_border) {
  spdlog::debug(
      "Updating the source border and the composite table layer "
      "with the seamline");
  OGRGeometryUniquePtr covered_geometry(nullptr), new_geometry(nullptr);
  for (const auto& seamline_feature : seamline_layer) {
    OGRGeometry* seamline_geometry(seamline_feature->GetGeometryRef());
    switch (seamline_feature->GetFieldAsInteger(0)) {
      case 100: {
        covered_geometry.reset(seamline_geometry->clone());
        break;
      }
      case 200: {
        new_geometry.reset(seamline_geometry->clone());
      }
    }
  }
  OGRGeometryUniquePtr diff_geometry(
      source_border->Difference(covered_geometry.get()));
  for (const auto& geometry : diff_geometry->toMultiPolygon())
    if (geometry->Intersect(new_geometry.get())) {
      source_border.reset(geometry->clone());
      break;
    }
  spdlog::debug("Updating the source border with the seamline - done");
  std::vector<std::pair<int, OGRGeometryUniquePtr>> dislocated_geometries;
  for (auto& composite_table_feature : composite_table_layer) {
    int id(static_cast<int>(composite_table_feature->GetFID()));
    OGRGeometry* composite_table_geometry(
        composite_table_feature->GetGeometryRef());
    if (covered_geometry->Intersect(composite_table_geometry)) {
      diff_geometry.reset(
          composite_table_geometry->Difference(source_border.get()));
      if (diff_geometry->getGeometryType() == wkbPolygon) {
        dislocated_geometries.emplace_back(id, diff_geometry->clone());
      } else {
        double max_area(0.0);
        bool be_head(true);
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
      composite_table_feature->SetGeometry(
          dislocated_geometries.back().second.get());
      composite_table_layer->SetFeature(composite_table_feature.get());
      dislocated_geometries.pop_back();
    }
  }
  for (const auto& info : dislocated_geometries)
    for (auto& composite_table_feature : composite_table_layer)
      if (info.first != composite_table_feature->GetFID() &&
          info.second->Touches(composite_table_feature->GetGeometryRef())) {
        composite_table_feature->SetGeometry(
            info.second->Union(composite_table_feature->GetGeometryRef()));
        composite_table_layer->SetFeature(composite_table_feature.get());
      }
  spdlog::debug("Updating the composite table layer with the seamline - done");
  spdlog::debug(
      "Updating the source border and the composite table layer "
      "with the seamline - done");
}

} // mosaicking
} // rs_toolset