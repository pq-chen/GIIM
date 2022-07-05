#include "mosaicking_container.h"

#include <cmath>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include "graph_cut_mosaicking.h"
#include "mosaicking_base.h"
#include <rs-toolset/utils.hpp>


namespace fs = std::filesystem;

namespace rs_toolset {
namespace mosaicking {

bool SortByReso(
    const RasterInfo& info1,
    const RasterInfo& info2) {
  if (info1.date && info2.date && info1.date != info2.date &&
      (fabs(info1.reso - info2.reso) / fmax(info1.reso, info2.reso)) < 0.01) {
    return info1.date > info2.date;
  } else {
    return info1.reso < info2.reso;
  }
}

MosaickingContainerImpl::MosaickingContainerImpl(
    std::shared_ptr<MosaickingInterface> mosaicking,
    OGRSpatialReference* spatial_ref)
    : mosaicking_(mosaicking),
      covered_border_(OGRGeometryFactory::createGeometry(wkbMultiPolygon)) {
  spdlog::info(
      "Creating a mosaicking container with\n - Spatial reference name: {}",
      spatial_ref->GetName());
  GDALDriver* memory_driver(GetGDALDriverManager()->GetDriverByName("Memory"));
  composite_table_dataset_ = memory_driver->Create(
      "", 0, 0, 0, GDT_Unknown, nullptr);
  composite_table_layer_ = composite_table_dataset_->CreateLayer(
      "", spatial_ref, wkbPolygon);
  OGRFieldDefn path_field("path", OFTString);
  composite_table_layer_->CreateField(&path_field);
}

MosaickingContainerImpl::MosaickingContainerImpl(
    std::shared_ptr<MosaickingInterface> mosaicking,
    OGRLayer* composite_table_layer,
    const std::string& rasters_dir)
    : MosaickingContainerImpl(
          mosaicking, composite_table_layer->GetSpatialRef()) {
  spdlog::info(
      "Initializing the mosaicking container by the external composite table");
  for (const auto& external_feature : composite_table_layer) {
    auto path(fs::path(rasters_dir) / external_feature->GetFieldAsString(0));
    if (!std::filesystem::exists(path)) {
      spdlog::warn("{} does not exist", path.string());
      continue;
    }
    spdlog::debug("Adding the mosaicking information from {}", path.string());

    // Create an internal feature
    spdlog::debug("Creating an internal feature");
    OGRFeatureUniquePtr internal_feature(OGRFeature::CreateFeature(
        composite_table_layer_->GetLayerDefn()));
    OGRGeometry* external_geometry(external_feature->GetGeometryRef());
    internal_feature->SetField(0, path.string().c_str());
    internal_feature->SetGeometry(external_geometry);
    composite_table_layer_->CreateFeature(internal_feature.get());
    spdlog::debug("Creating an internal feature - done");

    // Update the covered border
    spdlog::debug("Updating the covered border");
    if (covered_border_->IsEmpty() ||
        !external_geometry->Intersect(covered_border_)) {
      covered_border_->toMultiPolygon()->addGeometryDirectly(
          external_geometry->clone());
    } else {
      OGRGeometryUniquePtr new_covered_polygon(external_geometry->clone());
      OGRMultiPolygon* covered_multi_polygon(
          covered_border_->toMultiPolygon());
      for (int i = covered_multi_polygon->getNumGeometries() - 1;
          i >= 0; i--) {
        OGRGeometry* covered_polygon(covered_multi_polygon->getGeometryRef(i));
        if (new_covered_polygon->Touches(covered_polygon)) {
          new_covered_polygon.reset(
              new_covered_polygon->Union(covered_polygon));
          covered_multi_polygon->removeGeometry(i);
        }
      }
      covered_multi_polygon->addGeometryDirectly(
          new_covered_polygon.release());
    }
    spdlog::debug("Updating the covered border - done");
    spdlog::info(
        "Adding the mosaicking information from {} - done", path.string());
  }
  spdlog::info(
      "Initializing the mosaicking container "
      "by the external composite table - done");
}

MosaickingContainerImpl::~MosaickingContainerImpl() {
  OGRGeometryFactory::destroyGeometry(covered_border_);
  GDALClose(composite_table_dataset_);
}

void MosaickingContainerImpl::SortRasters(
    std::vector<std::string>& rasters_path,
    const SortFunc& sort_func) {
  spdlog::debug("Sorting rasters");
  std::vector<RasterInfo> rasters_info;
  rasters_info.reserve(rasters_path.size());
  for (const auto& raster_path : rasters_path) {
    GDALDatasetUniquePtr raster_dataset(GDALDataset::Open(
        raster_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    if (!raster_dataset) {
      spdlog::warn("Opening {} failed", raster_path);
      continue;
    }
    double geotrans[6];
    if (raster_dataset->GetGeoTransform(geotrans) != CE_None ||
        !raster_dataset->GetSpatialRef()) {
      spdlog::warn(
          "{} does not have the geotransform or the spatial reference. "
          "Please check whether the raster is DOM", raster_path);
      continue;
    }
    rasters_info.push_back({
        raster_path , geotrans[1],
        utils::DateMatching(fs::path(raster_path).stem().string()) });
  }
  std::sort(rasters_info.begin(), rasters_info.end(), sort_func);
  for (int i = 0; i < rasters_path.size(); i++)
    rasters_path[i] = rasters_info[i].path;
  spdlog::info("Sorting rasters - done");
}

bool MosaickingContainerImpl::AddTask(
    const std::string& raster_path,
    int last_over_view_idx,
    bool use_seamline) {
  return mosaicking_->Run(
      raster_path, composite_table_layer_, covered_border_, last_over_view_idx,
      use_seamline);
}

bool MosaickingContainerImpl::ExportCompositeTableVector(
    const std::string& composit_table_path,
    double buffer,
    double tol) {
  spdlog::debug(
      "Export the internal composite table to {}", composit_table_path);
  if (buffer > 0) {
    spdlog::warn("Buffer: {} is not accepted as a positive number");
    return false;
  }
  if (covered_border_->IsEmpty()) {
    spdlog::warn("No mosaicking information need to be exported");
    return false;
  }

  // Create an external composite table dataset
  spdlog::debug("Creating an external composite table dataset");
  GDALDriver* esri_shapefile_driver(
      GetGDALDriverManager()->GetDriverByName("ESRI Shapefile"));
  GDALDatasetUniquePtr composite_table_vector(esri_shapefile_driver->Create(
      composit_table_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
  if (!composite_table_vector) {
    spdlog::warn(
        "Creating an external composite table dataset from {} failed",
        composit_table_path);
    return false;
  }
  OGRLayer* external_layer(composite_table_vector->CreateLayer(
      "", composite_table_layer_->GetSpatialRef(), wkbPolygon));
  OGRFieldDefn name_field("name", OFTString);
  external_layer->CreateField(&name_field);
  spdlog::debug("Creating an external composite table - done");

  // Create the buffered or simplified covered border if needed
  OGRGeometryUniquePtr _covered_border(covered_border_->clone());
  if (buffer)
    _covered_border.reset(_covered_border->Buffer(buffer));
  if (tol)
    _covered_border.reset(_covered_border->Simplify(tol));

  // Create an external feature
  spdlog::debug("Creating an external feature");
  for (const auto& internal_feature : composite_table_layer_) {
    OGRFeatureUniquePtr external_feature(OGRFeature::CreateFeature(
        external_layer->GetLayerDefn()));
    fs::path path(internal_feature->GetFieldAsString(0));
    external_feature->SetField(0, path.filename().string().c_str());
    if (buffer || tol) {
      external_feature->SetGeometryDirectly(_covered_border->Intersection(
          internal_feature->GetGeometryRef()));
    } else {
      external_feature->SetGeometry(internal_feature->GetGeometryRef());
    }
    external_layer->CreateFeature(external_feature.get());
  }
  spdlog::debug("Creating an external feature - done");
  spdlog::info(
      "Export the internal composite table to {} - done", composit_table_path);
  return true;
}

std::vector<std::string> MosaickingContainerImpl::ExportAllRastersName() {
  std::vector<std::string> rasters_name;
  if (composite_table_layer_->GetFeatureCount() != 0) {
    spdlog::debug("Export all rasters' name");
    rasters_name.reserve(composite_table_layer_->GetFeatureCount());
    for (const auto& feature : composite_table_layer_) {
      fs::path path(feature->GetFieldAsString(0));
      rasters_name.push_back(path.filename().string());
    }
    spdlog::info("Export all rasters' name - done");
  }
  return rasters_name;
}

bool MosaickingContainerImpl::CreateMosaickingRaster(
    const std::string& mosaicking_raster_path,
    const std::string& composit_table_path,
    const std::string& rasters_dir,
    double reso,
    double blend_dist) {
  spdlog::debug(
      "Creating a mosaicking raster from {}",mosaicking_raster_path);
  GDALDatasetUniquePtr _composite_table_dataset(nullptr);
  OGRLayer* _composite_table_layer;
  if (composit_table_path.empty()) {
    _composite_table_layer = composite_table_layer_;
    spdlog::info("Using the internal composite table");
  } else {
    _composite_table_dataset.reset(GDALDataset::Open(
        composit_table_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
    if (!_composite_table_dataset) {
      spdlog::warn("Opening {} failed", composit_table_path);
      return false;
    }
    _composite_table_layer = _composite_table_dataset->GetLayer(0);
    spdlog::info(
        "Using the external composite table from {}", composit_table_path);
  }

  // Create a mosaicking raster dataset
  spdlog::debug("Creating a mosaicking raster dataset");
  OGREnvelope enve;
  covered_border_->getEnvelope(&enve);
  auto x_size(static_cast<int>(ceil((enve.MaxX - enve.MinX) / reso))),
      y_size(static_cast<int>(ceil((enve.MaxY - enve.MinY) / reso)));
  double geotrans[6]{ enve.MinX , reso, 0.0, enve.MaxY, 0.0, -reso };
  GDALDriver* gtiff_driver(GetGDALDriverManager()->GetDriverByName("GTiff"));
  GDALDatasetUniquePtr mosaicking_raster_dataset(gtiff_driver->Create(
      mosaicking_raster_path.c_str(), x_size, y_size, 3, GDT_Byte, nullptr));
  if (!mosaicking_raster_dataset) {
    spdlog::warn(
        "Creating a mosaicking raster dataset from {} failed",
        mosaicking_raster_path);
    return false;
  }
  mosaicking_raster_dataset->SetGeoTransform(geotrans);
  mosaicking_raster_dataset->SetSpatialRef(
      _composite_table_layer->GetSpatialRef());
  spdlog::info(
      "Creating a {}x{} mosaicking raster dataest - done", x_size, y_size);

  // Warp rasters by the corresponding geometry to the mosaicking raster
  for (const auto& feature : _composite_table_layer) {
    fs::path path;
    if (_composite_table_layer == composite_table_layer_) {
      path = feature->GetFieldAsString(0);
    } else {
      path = fs::path(rasters_dir) / feature->GetFieldAsString(0);
    }
    if (!fs::exists(path)) {
      spdlog::warn("{} does not exist", path.string());
      continue;
    }
    spdlog::info(
        "Warping {} by the corresponding geometry to the mosaicking raster",
        path.string());
    GDALDatasetUniquePtr source_raster_dataset(GDALDataset::Open(
        path.string().c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    utils::WarpByGeometry(
        { source_raster_dataset.get() }, { feature->GetGeometryRef() },
        mosaicking_raster_dataset, GRA_Bilinear, blend_dist);
    spdlog::info(
        "Warping {} by the corresponding geometry to the mosaicking raster"
        " - done", path.string());
    spdlog::info(
        "----------- {}/{} - done ----------",
        feature->GetFID() + 1, _composite_table_layer->GetFeatureCount());
  }
  spdlog::info(
      "Creating a mosaicking raster from {} - done", mosaicking_raster_path);
  mosaicking_raster_dataset.reset(nullptr);
  mosaicking_raster_dataset.reset(GDALDataset::Open(
      mosaicking_raster_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  utils::CreateRasterPyra(mosaicking_raster_dataset.get());
  return true;
}

std::shared_ptr<MosaickingContainer> MosaickingContainer::Create(
    std::shared_ptr<MosaickingInterface> mosaicking,
    OGRSpatialReference* spatial_ref) {
  if (!spatial_ref) {
    spdlog::warn("No spatial reference is specified for the seamline layer");
    return nullptr;
  }
  return std::make_shared<MosaickingContainerImpl>(mosaicking, spatial_ref);
}

std::shared_ptr<MosaickingContainer> MosaickingContainer::Create(
    std::shared_ptr<MosaickingInterface> mosaicking,
    const std::string& composite_tabel_path,
    const std::string& rasters_dir) {
  GDALDatasetUniquePtr composite_tabel_dataset(GDALDataset::Open(
      composite_tabel_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
  if (!composite_tabel_dataset) {
    spdlog::warn("Opening {} failed", composite_tabel_path);
    return nullptr;
  }
  OGRLayer* composite_table_layer(composite_tabel_dataset->GetLayer(0));
  if (!composite_table_layer->GetSpatialRef()) {
    spdlog::warn("No spatial reference is specified for the seamline layer");
    return nullptr;
  }
  return std::make_shared<MosaickingContainerImpl>(
      mosaicking, composite_table_layer, rasters_dir);
}

} // mosaicking
} // rs_toolset