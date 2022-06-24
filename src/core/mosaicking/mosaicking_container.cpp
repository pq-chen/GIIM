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


namespace rs_toolset {
namespace mosaicking {

bool SortByReso(
    const RasterInfo& info1,
    const RasterInfo& info2) {
  if (!info1.date || !info2.date || 
      (fabs(info1.reso - info2.reso) / fmax(info1.reso, info2.reso)) > 0.01) {
    return info1.reso < info2.reso;
  } else {
    return info1.date > info2.date;
  }
}

MosaickingContainerImpl::MosaickingContainerImpl(
    std::shared_ptr<MosaickingInterface> mosaicking,
    OGRSpatialReference* spatial_ref)
    : mosaicking_(mosaicking),
      covered_border_(OGRGeometryFactory::createGeometry(wkbMultiPolygon)) {
  spdlog::info(
      "Creating the mosaicking container with\n- Spatial reference name: {}",
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
  for (const auto& _composite_table_feature : composite_table_layer) {
    auto path(
        std::filesystem::path(rasters_dir) / 
        _composite_table_feature->GetFieldAsString(0));
    if (!std::filesystem::exists(path)) {
      spdlog::warn("{} does not exist", path.string());
      continue;
    }
    spdlog::debug("Adding the mosaicking info for {}", path.string());

    // Creating the feature for the internal composite table
    OGRFeatureUniquePtr composite_table_feature(OGRFeature::CreateFeature(
        composite_table_layer_->GetLayerDefn()));
    OGRGeometry* _composite_table_geometry(
        _composite_table_feature->GetGeometryRef());
    composite_table_feature->SetField(0, path.string().c_str());
    composite_table_feature->SetGeometry(_composite_table_geometry);
    composite_table_layer_->CreateFeature(composite_table_feature.get());
    spdlog::debug(
        "Creating the feature for the internal composite table - done");

    // Updating the covered border
    if (covered_border_->IsEmpty()) {
      covered_border_->toMultiPolygon()->addGeometryDirectly(
          _composite_table_geometry->clone());
    } else {
      OGRGeometryUniquePtr new_covered_polygon(
          _composite_table_geometry->clone());
      OGRMultiPolygon* _covered_border_(covered_border_->toMultiPolygon());
      for (int i = _covered_border_->getNumGeometries() - 1; i >= 0; i--) {
        OGRGeometry* covered_polygon(_covered_border_->getGeometryRef(i));
        if (new_covered_polygon->Touches(covered_polygon)) {
          new_covered_polygon.reset(
              new_covered_polygon->Union(covered_polygon));
          _covered_border_->removeGeometry(i);
        }
      }
      _covered_border_->addGeometryDirectly(new_covered_polygon.get());
      new_covered_polygon.release();
    }
    spdlog::debug("Updating the covered border - done");
    spdlog::info("Adding the mosaicking info for {} - done", path.string());
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
  spdlog::debug("Sorting the given rasters");
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
        raster_path , geotrans[1], utils::DateMatching(
            std::filesystem::path(raster_path).stem().string()) });
  }
  std::sort(rasters_info.begin(), rasters_info.end(), sort_func);
  for (int i = 0; i < rasters_path.size(); i++)
    rasters_path[i] = rasters_info[i].path;
  spdlog::info("Sorting the given rasters - done");
}

bool MosaickingContainerImpl::AddTask(
    const std::string& raster_path,
    bool use_seamline) {
  return mosaicking_->Run(
      raster_path, composite_table_layer_, covered_border_, use_seamline);
}

bool MosaickingContainerImpl::ExportCompositeTableVector(
    const std::string& composit_table_path,
    double unit,
    double buffer,
    double tol) {
  spdlog::debug(
      "Exporting the composite table vector from {}", composit_table_path);
  if (covered_border_->IsEmpty()) {
    spdlog::warn("No mosaicking info need to be exported");
    return false;
  }

  // Create the external composite table dataset
  GDALDriver* esri_shapefile_driver(
      GetGDALDriverManager()->GetDriverByName("ESRI Shapefile"));
  GDALDatasetUniquePtr composite_table_vector(esri_shapefile_driver->Create(
      composit_table_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
  if (!composite_table_vector) {
    spdlog::warn(
        "Creating the composite table from {} failed", composit_table_path);
    return false;
  }
  OGRLayer* composite_table_layer(composite_table_vector->CreateLayer(
      "", composite_table_layer_->GetSpatialRef(), wkbPolygon));
  OGRFieldDefn name_field("name", OFTString);
  composite_table_layer->CreateField(&name_field);

  // Create the buffered or simplified covered border if needed
  OGRGeometryUniquePtr _covered_border(nullptr);
  if (buffer || tol) {
    _covered_border.reset(covered_border_->Buffer(buffer * unit));
    _covered_border.reset(_covered_border->Simplify(tol * unit));
  }

  // Create the feature for the external composite table
  for (const auto& composite_table_feature : composite_table_layer_) {
    OGRFeatureUniquePtr feature(OGRFeature::CreateFeature(
        composite_table_layer->GetLayerDefn()));
    std::filesystem::path path(composite_table_feature->GetFieldAsString(0));
    feature->SetField(0, path.filename().string().c_str());
    if (buffer || tol) {
      feature->SetGeometryDirectly(_covered_border->Intersection(
          composite_table_feature->GetGeometryRef()));
    } else {
      feature->SetGeometry(composite_table_feature->GetGeometryRef());
    }
    composite_table_layer->CreateFeature(feature.get());
  }
  spdlog::info(
      "Exporting the composite table vector from {} - done", 
      composit_table_path);
  return true;
}

std::vector<std::string> MosaickingContainerImpl::ExportAllRastersName() {
  std::vector<std::string> rasters_name;
  rasters_name.reserve(composite_table_layer_->GetFeatureCount());
  for (const auto& composite_table_feature : composite_table_layer_) {
    std::filesystem::path path(composite_table_feature->GetFieldAsString(0));
    rasters_name.push_back(path.filename().string());
  }
  return rasters_name;
}

bool MosaickingContainerImpl::CreateMosaickingRaster(
    const std::string& mosaicking_raster_path,
    const std::string& composit_table_path,
    const std::string& rasters_dir,
    double reso) {
  spdlog::debug(
      "Creating the mosaicking raster from {}", mosaicking_raster_path);
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

  // Create the mosaicking raster dataset
  OGREnvelope enve;
  covered_border_->getEnvelope(&enve);
  int x_size(static_cast<int>(ceil((enve.MaxX - enve.MinX) / reso))),
      y_size(static_cast<int>(ceil((enve.MaxY - enve.MinY) / reso)));
  double geotrans[6]{ enve.MinX , reso, 0, enve.MaxY, 0, -reso };
  GDALDriver* gtiff_driver(GetGDALDriverManager()->GetDriverByName("GTiff"));
  GDALDatasetUniquePtr mosaicking_raster_dataset(gtiff_driver->Create(
      mosaicking_raster_path.c_str(), x_size, y_size, 3, GDT_Byte, nullptr));
  if (mosaicking_raster_dataset) {
    spdlog::info(
        "Creating the {}x{} mosaicking raster from {} - done", 
        x_size, y_size, mosaicking_raster_path);
  } else {
    spdlog::warn(
        "Creating the mosaicking raster from {} failed",
        mosaicking_raster_path);
    return false;
  }
  mosaicking_raster_dataset->SetGeoTransform(geotrans);
  mosaicking_raster_dataset->SetSpatialRef(
      _composite_table_layer->GetSpatialRef());

  // Warp rasters by the corresponding geometry to the mosaicking raster
  for (const auto& composite_table_feature : _composite_table_layer) {
    std::filesystem::path path;
    if (_composite_table_layer == composite_table_layer_) {
      path = composite_table_feature->GetFieldAsString(0);
    } else {
      path = std::filesystem::path(rasters_dir) / 
          composite_table_feature->GetFieldAsString(0);
    }
    if (std::filesystem::exists(path)) {
      spdlog::info(
          "Warping {} by the corresponding geometry to the mosaicking raster", 
          path.string());
      GDALDatasetUniquePtr source_raster_dataset(GDALDataset::Open(
          path.string().c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
      utils::WarpByGeometry(
          { source_raster_dataset.get() }, mosaicking_raster_dataset.get(),
          composite_table_feature->GetGeometryRef());
      spdlog::info(
          "Warping {} by the corresponding geometry to the mosaicking raster"
          " - done", path.string());
      spdlog::info(
          "----------- {}/{} - done ----------", 
          composite_table_feature->GetFID() + 1, 
          _composite_table_layer->GetFeatureCount());
    } else {
      spdlog::warn("{} does not exist", path.string());
    }
  }
  spdlog::info(
      "Creating the mosaicking raster from {} - done", mosaicking_raster_path);
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