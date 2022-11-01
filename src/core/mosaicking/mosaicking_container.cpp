#include "mosaicking_container.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/mosaicking.h>
#include <rs-toolset/utils.hpp>

namespace fs = std::filesystem;

namespace rs_toolset {
namespace mosaicking {

bool SortByReso(const RasterInfo& info1, const RasterInfo& info2) {
  if (info1.date && info2.date && info1.date != info2.date &&
      (fabs(info1.reso - info2.reso) / fmax(info1.reso, info2.reso)) < 0.01) {
    return info1.date > info2.date;
  } else {
    return info1.reso < info2.reso;
  }
}

MosaickingContainerImpl::MosaickingContainerImpl(
    const std::shared_ptr<MosaickingInterface>& mosaicking,
    OGRSpatialReference* spatial_ref,
    double rejection_ratio,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing)
    : mosaicking_(mosaicking),
      covered_border_(new OGRMultiPolygon),
      rejection_ratio_(rejection_ratio),
      color_balancing_(color_balancing) {
  spdlog::info(
      "Creating a mosaicking container with\n"
      " - Spatial reference name: {}\n"
      " - Rejection ratio: {}\n"
      " - With color balancing: {}", spatial_ref->GetName(), rejection_ratio,
      bool(color_balancing));
  auto driver(GetGDALDriverManager()->GetDriverByName("Memory"));
  composite_table_dataset_ = driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
  composite_table_layer_ = composite_table_dataset_->CreateLayer(
      "", spatial_ref, wkbPolygon);
  OGRFieldDefn path_field("path", OFTString);
  composite_table_layer_->CreateField(&path_field);
  if (color_balancing_) {
    OGRFieldDefn cb_index_field("CB index", OFTInteger);
    composite_table_layer_->CreateField(&cb_index_field);
  }
  border_dataset_ = driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
  border_layer_ = border_dataset_->CreateLayer("", spatial_ref, wkbPolygon);
}

MosaickingContainerImpl::~MosaickingContainerImpl() {
  GDALClose(composite_table_dataset_);
  GDALClose(border_dataset_);
  OGRGeometryFactory::destroyGeometry(covered_border_);
}

void MosaickingContainerImpl::InitializeByExt(
    OGRLayer* composite_table_layer,
    const std::string& rasters_dir) {
  spdlog::info(
      "Initializing the mosaicking container by the external composite table");
  auto with_cb_index(
      composite_table_layer->FindFieldIndex("CB index", true) == 1);
  for (const auto& ext_feature : composite_table_layer) {
    auto path(rasters_dir + "/" + ext_feature->GetFieldAsString(0));
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    if (!dataset) {
      spdlog::error("Opening {} failed", path);
      continue;
    }

    spdlog::debug("Adding the mosaicking information from {}", path);
    OGRFeatureUniquePtr int_feature(OGRFeature::CreateFeature(
        composite_table_layer_->GetLayerDefn()));
    OGRGeometryUniquePtr ext_geometry(ext_feature->StealGeometry());
    int_feature->SetGeometry(ext_geometry.get());
    int_feature->SetField(0, path.c_str());
    if (color_balancing_) {
      int_feature->SetField(
          1, with_cb_index ? ext_feature->GetFieldAsInteger(1) : -1);
    }
    composite_table_layer_->CreateFeature(int_feature.get());

    // Update the border layer
    spdlog::debug("Updating the border layer");
    OGRFeatureUniquePtr border_feature(OGRFeature::CreateFeature(
        border_layer_->GetLayerDefn()));
    border_feature->SetGeometryDirectly(utils::CreateBorder(
        dataset.get(), composite_table_layer_->GetSpatialRef()).release());
    border_layer_->CreateFeature(border_feature.get());
    borders_area_.push_back(
        border_feature->GetGeometryRef()->toPolygon()->get_Area());
    spdlog::debug("Updating the border layer - done");

    // Update the covered border
    spdlog::debug("Updating the covered border");
    if (auto _covered_border_(covered_border_->toMultiPolygon());
        covered_border_->IsEmpty() ||
        !ext_geometry->Intersect(_covered_border_)) {
      _covered_border_->addGeometryDirectly(ext_geometry.release());
    } else {
      for (int i(_covered_border_->getNumGeometries() - 1); i >= 0; i--) {
        if (auto covered_polygon(_covered_border_->getGeometryRef(i));
            ext_geometry->Intersect(covered_polygon)) {
          ext_geometry.reset(ext_geometry->Union(covered_polygon));
          _covered_border_->removeGeometry(i);
        }
      }
      _covered_border_->addGeometryDirectly(ext_geometry.release());
    }
    spdlog::debug("Updating the covered border - done");
    spdlog::info("Adding the mosaicking information from {} - done", path);
  }
  spdlog::info(
      "Initializing the mosaicking container "
      "by the external composite table - done");
}

bool MosaickingContainerImpl::SortRasters(
    std::vector<std::string>& paths,
    const SortFunc& sort_func) {
  spdlog::info("Sorting rasters");
  std::vector<RasterInfo> rasters_info;
  for (const auto& path : paths) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    if (!dataset) {
      spdlog::error("Opening {} failed", path);
      return false;
    }
    double geotrans[6];
    if (dataset->GetGeoTransform(geotrans) != CE_None ||
        !dataset->GetSpatialRef()) {
      spdlog::error(
          "{} does not have the geotransform or the spatial reference. "
          "Please check whether the raster is DOM", path);
      return false;
    }

    rasters_info.push_back({
        path, geotrans[1],
        utils::DateMatching(fs::path(path).stem().string())});
  }
  std::sort(rasters_info.begin(), rasters_info.end(), sort_func);
  for (int i(0); i < paths.size(); ++i)
    paths[i] = rasters_info[i].path;
  spdlog::info("Sorting rasters - done");
  return true;
}

bool MosaickingContainerImpl::AddTask(
    const std::string& path,
    int low_overview_trunc,
    int high_overview_trunc,
    const std::vector<int>& rgb_bands_map) {
  return mosaicking_->RunTaskForExisting(
      path, composite_table_layer_, border_layer_, covered_border_,
      borders_area_, rejection_ratio_, low_overview_trunc, high_overview_trunc,
      rgb_bands_map, color_balancing_);
}

GDALDatasetUniquePtr MosaickingContainerImpl::ExportCompositeTable(
    const std::string& output_path,
    const std::string& query_path,
    const std::string& query_rasters_name_field_name,
    bool with_extension) {
  spdlog::info("Exporting the internal composite table to {}", output_path);
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
  if (covered_border_->IsEmpty()) {
    spdlog::warn("No mosaicking information can be exported");
    return nullptr;
  }
  GDALDatasetUniquePtr query_dataset(nullptr);
  OGRLayer* query_layer(nullptr);
  if (!query_path.empty()) {
    query_dataset.reset(GDALDataset::Open(
        query_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
    if (!query_dataset) {
      spdlog::error("Opening {} failed", query_path);
      return nullptr;
    } else if (query_dataset->GetLayer(0)->FindFieldIndex(
            query_rasters_name_field_name.c_str(), true) == -1) {
      spdlog::error(
          "The query compostite table doesn't contain the \"{}\" field name",
          query_rasters_name_field_name);
      return nullptr;
    } else {
      query_layer = query_dataset->GetLayer(0);
      spdlog::info("Using the query composite table from {}", query_path);
    }
  }

  char** options(nullptr);
  options = CSLSetNameValue(options, "ENCODING", "LDID/77");
  if (auto ext_layer(output_dataset->CreateLayer(
          "", composite_table_layer_->GetSpatialRef(), wkbPolygon, options));
      query_layer) {
    // Create fields from the query layer definition
    auto fields_count(query_layer->GetLayerDefn()->GetFieldCount());
    for (int i(0); i < fields_count; ++i)
      ext_layer->CreateField(query_layer->GetLayerDefn()->GetFieldDefn(i));
    std::map<std::string, int> query_raster_name_to_idx;
    for (int i(0); i < query_layer->GetFeatureCount(); ++i) {
      OGRFeatureUniquePtr feature(query_layer->GetFeature(i));
      query_raster_name_to_idx[feature->GetFieldAsString(
          query_rasters_name_field_name.c_str())] = i;
    }

    for (const auto& int_feature : composite_table_layer_) {
      OGRFeatureUniquePtr ext_feature(OGRFeature::CreateFeature(
          ext_layer->GetLayerDefn()));
      ext_feature->SetGeometryDirectly(int_feature->StealGeometry());
      auto filename(fs::path(int_feature->GetFieldAsString(0)).filename());
      if (auto it(query_raster_name_to_idx.find(filename.string()));
          it != query_raster_name_to_idx.end()) {
        OGRFeatureUniquePtr feature(query_layer->GetFeature(it->second));
        for (int i(0); i < fields_count; ++i)
          ext_feature->SetField(i, feature->GetRawFieldRef(i));
      } else if (with_extension) {
        ext_feature->SetField(
            query_rasters_name_field_name.c_str(), filename.string().c_str());
      } else {
        ext_feature->SetField(
            query_rasters_name_field_name.c_str(),
            filename.stem().string().c_str());
      }
      ext_layer->CreateFeature(ext_feature.get());
    }
  } else {
    OGRFieldDefn name_field("name", OFTString);
    ext_layer->CreateField(&name_field);
    if (color_balancing_) {
      ext_layer->CreateField(
          composite_table_layer_->GetLayerDefn()->GetFieldDefn(1));
    }

    for (const auto& int_feature : composite_table_layer_) {
      OGRFeatureUniquePtr ext_feature(OGRFeature::CreateFeature(
          ext_layer->GetLayerDefn()));
      ext_feature->SetGeometryDirectly(int_feature->StealGeometry());
      fs::path path(int_feature->GetFieldAsString(0));
      ext_feature->SetField(
          0,
          with_extension ? path.filename().string().c_str()
              : path.stem().string().c_str());
      if (color_balancing_)
        ext_feature->SetField(1, int_feature->GetFieldAsInteger(1));
      ext_layer->CreateFeature(ext_feature.get());
    }
  }
  CSLDestroy(options);
  spdlog::info("Exporting the internal composite table - done");
  return output_dataset;
}

std::vector<std::string> MosaickingContainerImpl::ExportAllRastersName() {
  std::vector<std::string> names;
  for (const auto& feature : composite_table_layer_)
    names.push_back(fs::path(feature->GetFieldAsString(0)).filename().string());
  return names;
}

std::shared_ptr<MosaickingContainer> MosaickingContainer::Create(
    const std::shared_ptr<MosaickingInterface>& mosaicking,
    OGRSpatialReference* spatial_ref,
    double rejection_ratio,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  if (!spatial_ref) {
    spdlog::error("The spatial reference can not be empty");
    return nullptr;
  }

  return std::make_shared<MosaickingContainerImpl>(
      mosaicking, spatial_ref, rejection_ratio, color_balancing);
}

std::shared_ptr<MosaickingContainer> MosaickingContainer::Create(
    const std::shared_ptr<MosaickingInterface>& mosaicking,
    const std::string& composite_tabel_path,
    const std::string& rasters_dir,
    double rejection_ratio,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  GDALDatasetUniquePtr dataset(GDALDataset::Open(
      composite_tabel_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
  if (!dataset) {
    spdlog::error("Opening {} failed", composite_tabel_path);
    return nullptr;
  }
  auto layer(dataset->GetLayer(0));
  if (!layer->GetSpatialRef()) {
    spdlog::error(
        "{} does not have the spatial reference", composite_tabel_path);
    return nullptr;
  }

  auto mosaicking_container(std::make_shared<MosaickingContainerImpl>(
      mosaicking, layer->GetSpatialRef(), rejection_ratio, color_balancing));
  mosaicking_container->InitializeByExt(layer, rasters_dir);
  return mosaicking_container;
}

}  // namespace mosaicking
}  // namespace rs_toolset