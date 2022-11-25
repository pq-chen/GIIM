#include "serial_container.h"

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
      (fabs(info1.reso - info2.reso) / fmax(info1.reso, info2.reso)) < 0.1) {
    return info1.date > info2.date;
  } else {
    return info1.reso < info2.reso;
  }
}

SerialContainerImpl::SerialContainerImpl(
    const std::shared_ptr<MosaickingInterface>& mosaicking,
    OGRSpatialReference* spatial_ref,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing)
    : mosaicking_(mosaicking),
      covered_border_(new OGRMultiPolygon),
      color_balancing_(color_balancing) {
  spdlog::info(
      "Creating a mosaicking container with\n"
      " - Spatial reference name: {}\n"
      " - With color balancing: {}", spatial_ref->GetName(),
      bool(color_balancing));
  auto driver(GetGDALDriverManager()->GetDriverByName("Memory"));
  composite_table_dataset_ = driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
  composite_table_layer_ = composite_table_dataset_->CreateLayer(
      "", spatial_ref, wkbPolygon);
  OGRFieldDefn path_field("path", OFTString);
  composite_table_layer_->CreateField(&path_field);
  if (color_balancing_) {
    OGRFieldDefn cb_idx_field("CB_idx", OFTInteger);
    composite_table_layer_->CreateField(&cb_idx_field);
  }
  border_dataset_ = driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
  border_layer_ = border_dataset_->CreateLayer("", spatial_ref, wkbPolygon);
  border_layer_->CreateField(&path_field);
}

SerialContainerImpl::~SerialContainerImpl() {
  GDALClose(composite_table_dataset_);
  GDALClose(border_dataset_);
  OGRGeometryFactory::destroyGeometry(covered_border_);
}

void SerialContainerImpl::InitializeByExt(
    OGRLayer* composite_table_layer,
    const std::string& rasters_dir) {
  spdlog::info(
      "Initializing the mosaicking container by the external composite table");
  auto with_cb_idx(
      composite_table_layer->FindFieldIndex("CB_idx", true) == 1);
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
          1, with_cb_idx ? ext_feature->GetFieldAsInteger(1) : -1);
    }
    int_feature->SetFID(composite_table_layer_->GetFeatureCount());
    composite_table_layer_->CreateFeature(int_feature.get());

    // Update the border layer
    spdlog::debug("Updating the border layer");
    OGRFeatureUniquePtr border_feature(OGRFeature::CreateFeature(
        border_layer_->GetLayerDefn()));
    border_feature->SetGeometryDirectly(utils::CreateBorder(
        dataset.get(), composite_table_layer_->GetSpatialRef()).release());
    border_feature->SetField(0, path.c_str());
    border_feature->SetFID(border_layer_->GetFeatureCount());
    border_layer_->CreateFeature(border_feature.get());
    borders_area_.push_back(
        border_feature->GetGeometryRef()->toPolygon()->get_Area());
    spdlog::debug("Updating the border layer - done");

    // Update the covered border
    spdlog::debug("Updating the covered border");
    if (auto _covered_border_(covered_border_->toMultiPolygon());
        !ext_geometry->Intersect(_covered_border_)) {
      _covered_border_->addGeometryDirectly(ext_geometry.release());
    } else {
      for (int i(_covered_border_->getNumGeometries() - 1); i >= 0; --i) {
        if (auto covered_polygon(_covered_border_->getGeometryRef(i));
            ext_geometry->Intersect(covered_polygon)) {
          ext_geometry.reset(ext_geometry->Union(covered_polygon));
          _covered_border_->removeGeometry(i);
        }
      }
      switch (ext_geometry->getGeometryType()) {
        case wkbPolygon: {
          auto _ext_geometry(ext_geometry->toPolygon());
          for (int i(_ext_geometry->getNumInteriorRings() - 1); i >= 0; --i) {
            if (_ext_geometry->getInteriorRing(i)->get_Area() == 0)
              _ext_geometry->removeRing(i + 1);
          }
          break;
        }
        case wkbMultiPolygon: {
          utils::ExtractBiggestPolygon(ext_geometry);
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

bool SerialContainerImpl::SortRasters(
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

bool SerialContainerImpl::AddTask(
    const std::string& path,
    int low_overviews_trunc,
    int high_overviews_trunc,
    double surrounded_buffer,
    double rejection_ratio,
    bool anti_surrounded,
    const std::vector<int>& rgb_bands_map) {
  return mosaicking_->RunTaskForSerial(
      path, composite_table_layer_, border_layer_, covered_border_,
      borders_area_, low_overviews_trunc, high_overviews_trunc,
      surrounded_buffer, rejection_ratio, anti_surrounded, rgb_bands_map,
      color_balancing_);
}

GDALDatasetUniquePtr SerialContainerImpl::ExportMosaickingVector(
    const std::string& output_path,
    const std::string& query_path,
    const std::string& query_raster_name_field_name,
    bool extension) {
  spdlog::info(
      "Exporting the internal composite table to {}\n - Extension: {}",
      output_path, extension);
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
  if (composite_table_layer_->GetFeatureCount() == 0) {
    spdlog::warn("No mosaicking information can be exported");
    return nullptr;
  }
  GDALDatasetUniquePtr query_dataset(nullptr);
  OGRLayer* query_layer(nullptr);
  int query_rasters_name_field_idx;
  if (!query_path.empty()) {
    query_dataset.reset(GDALDataset::Open(
        query_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
    if (!query_dataset) {
      spdlog::error("Opening {} failed", query_path);
      return nullptr;
    } else {
      query_layer = query_dataset->GetLayer(0);
      auto recoded_field_name(utils::CreateCplString(CPLRecode(
          query_raster_name_field_name.c_str(), "CP936",
          CSLFetchNameValue(
              query_layer->GetMetadata("SHAPEFILE"), "SOURCE_ENCODING"))));
      if (query_rasters_name_field_idx = query_layer->FindFieldIndex(
              recoded_field_name.get(), true);
          query_rasters_name_field_idx == -1) {
        spdlog::error(
            "The query compostite table doesn't contain the \"{}\" field name",
            query_raster_name_field_name);
        return nullptr;
      }
      spdlog::info("Using the query composite table from {}", query_path);
    }
  }

  auto options(utils::CreateCslStringList(CSLSetNameValue(
      nullptr, "ENCODING", "CP936")));
  auto ext_layer(output_dataset->CreateLayer(
      "", composite_table_layer_->GetSpatialRef(), wkbPolygon, options.get()));
  if (query_layer) {
    // Create fields from the query layer definition
    auto fields_count(query_layer->GetLayerDefn()->GetFieldCount());
    for (int i(0); i < fields_count; ++i)
      ext_layer->CreateField(query_layer->GetLayerDefn()->GetFieldDefn(i));
    std::map<std::string, int> query_raster_name_to_idx;
    for (int i(0); i < query_layer->GetFeatureCount(); ++i) {
      OGRFeatureUniquePtr feature(query_layer->GetFeature(i));
      query_raster_name_to_idx[feature->GetFieldAsString(
          query_rasters_name_field_idx)] = i;
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
      } else if (extension) {
        ext_feature->SetField(
            query_rasters_name_field_idx, filename.string().c_str());
      } else {
        ext_feature->SetField(
            query_rasters_name_field_idx, filename.stem().string().c_str());
      }
      ext_feature->SetFID(ext_layer->GetFeatureCount());
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
          extension ? path.filename().string().c_str()
              : path.stem().string().c_str());
      if (color_balancing_)
        ext_feature->SetField(1, int_feature->GetFieldAsInteger(1));
      ext_feature->SetFID(ext_layer->GetFeatureCount());
      ext_layer->CreateFeature(ext_feature.get());
    }
  }
  spdlog::info("Exporting the internal composite table - done");
  return output_dataset;
}

GDALDatasetUniquePtr SerialContainerImpl::ExportBorderVector(
    const std::string& output_path) {
  spdlog::info("Exporting the internal border vector to {}", output_path);
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
  if (border_layer_->GetFeatureCount() == 0) {
    spdlog::warn("No border information can be exported");
    return nullptr;
  }

  auto options(utils::CreateCslStringList(CSLSetNameValue(
      nullptr, "ENCODING", "CP936")));
  auto ext_layer(output_dataset->CreateLayer(
      "", composite_table_layer_->GetSpatialRef(), wkbPolygon, options.get()));
  OGRFieldDefn name_field("name", OFTString);
  ext_layer->CreateField(&name_field);
  for (const auto& int_feature : border_layer_) {
    OGRFeatureUniquePtr ext_feature(OGRFeature::CreateFeature(
        ext_layer->GetLayerDefn()));
    ext_feature->SetGeometryDirectly(int_feature->StealGeometry());
    ext_feature->SetField(
        0,
        fs::path(int_feature->GetFieldAsString(0)).filename().string().c_str());
    ext_feature->SetFID(ext_layer->GetFeatureCount());
    ext_layer->CreateFeature(ext_feature.get());
  }
  spdlog::info("Exporting the internal border vector - done");
  return output_dataset;
}

std::vector<std::string> SerialContainerImpl::ExportAllRastersPath() {
  std::vector<std::string> paths;
  for (const auto& feature : composite_table_layer_)
    paths.push_back(fs::path(feature->GetFieldAsString(0)).string());
  return paths;
}

std::shared_ptr<SerialContainer> SerialContainer::Create(
    const std::shared_ptr<MosaickingInterface>& mosaicking,
    OGRSpatialReference* spatial_ref,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  if (!spatial_ref) {
    spdlog::error("The spatial reference can not be empty");
    return nullptr;
  }

  return std::make_shared<SerialContainerImpl>(
      mosaicking, spatial_ref, color_balancing);
}

std::shared_ptr<SerialContainer> SerialContainer::Create(
    const std::shared_ptr<MosaickingInterface>& mosaicking,
    const std::string& mosaicking_vector_path,
    const std::string& rasters_dir,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  GDALDatasetUniquePtr dataset(GDALDataset::Open(
      mosaicking_vector_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
  if (!dataset) {
    spdlog::error("Opening {} failed", mosaicking_vector_path);
    return nullptr;
  }
  auto layer(dataset->GetLayer(0));
  if (!layer->GetSpatialRef()) {
    spdlog::error(
        "{} does not have the spatial reference", mosaicking_vector_path);
    return nullptr;
  }

  auto mosaicking_container(std::make_shared<SerialContainerImpl>(
      mosaicking, layer->GetSpatialRef(), color_balancing));
  mosaicking_container->InitializeByExt(layer, rasters_dir);
  return mosaicking_container;
}

}  // namespace mosaicking
}  // namespace rs_toolset