#include <rs-toolset/mosaicking.h>

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gdalwarper.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/utils.hpp>

namespace fs = std::filesystem;

namespace {

using namespace rs_toolset;

bool CreateMosaickingRasterCheck(
    OGRLayer* layer,
    const std::vector<std::string>& input_paths,
    const std::string& output_path,
    double reso,
    const std::string& compression,
    double blend_dist,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing,
    GDALDatasetUniquePtr& output_dataset) {
  if (!layer) {
    spdlog::error("The \"layer\" is empty");
    return false;
  }
  if (layer->GetFeatureCount() != input_paths.size()) {
    spdlog::error(
        "The features' count {} in the \"layer\" must be the same with the "
        "\"input_paths\" size {}", layer->GetFeatureCount(),
        input_paths.size());
    return false;
  }
  for (const auto& path : input_paths) {
    if (!path.empty()) {
      if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
              path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)); !dataset) {
        spdlog::error("Opening {} failed", path);
        return false;
      }
    }
  }
  auto driver(utils::GetRasterDriverByPath(output_path, compression));
  if (!driver) {
    spdlog::error("Finding the driver for {} failed", output_path);
    return false;
  }
  if (reso <= 0) {
    spdlog::error("The \"reso\" {} must be positive", reso);
    return false;
  }
  if (blend_dist < 0) {
    spdlog::error("The \"blend_dist\" {} must be non-negative", blend_dist);
    return false;
  }
  if (color_balancing && (layer->GetLayerDefn()->GetFieldCount() != 2 ||
          layer->GetLayerDefn()->GetFieldDefn(1)->GetNameRef() != "CB_idx")) {
    spdlog::error("The \"layer\" is not qualified for color balancing");
    return false;
  }
  OGREnvelope enve;
  layer->GetExtent(&enve);
  auto x_size(static_cast<int>(ceil((enve.MaxX - enve.MinX) / reso))),
      y_size(static_cast<int>(ceil((enve.MaxY - enve.MinY) / reso)));
  spdlog::info(
      "Creating a {}x{} mosaicking raster with the {} compression method to {}",
      x_size, y_size, compression.empty() ? "no" : compression, output_path);
  std::unique_ptr<char*, void(*)(char**)> options(nullptr, nullptr);
  if (!output_path.empty() && !compression.empty()) {
    options = utils::CreateCslStringList(CSLSetNameValue(
        nullptr, "COMPRESS", compression.c_str()));
  }
  if (output_dataset.reset(driver->Create(
          output_path.c_str(), x_size, y_size, 3, GDT_Byte, options.get()));
      !output_dataset) {
    spdlog::error("Creating {} failed", output_path);
    return false;
  }
  double geotrans[6]{enve.MinX, reso, 0.0, enve.MaxY, 0.0, -reso};
  output_dataset->SetGeoTransform(geotrans);
  output_dataset->SetSpatialRef(layer->GetSpatialRef());
  return true;
}

bool CreateRastersCutCheck(
    OGRLayer* layer,
    const std::vector<std::string>& input_paths,
    const std::string& output_dir,
    const std::string& raster_name_field_name,
    const std::string& compression,
    int& raster_name_field_idx) {
  if (!layer) {
    spdlog::error("The \"layer\" is empty");
    return false;
  }
  if (layer->GetFeatureCount() != input_paths.size()) {
    spdlog::error(
        "The features' count {} in the \"layer\" must be the same with the "
        "\"input_paths\" size {}", layer->GetFeatureCount(),
        input_paths.size());
    return false;
  }
  for (const auto& path : input_paths) {
    if (!path.empty()) {
      if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
              path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)); !dataset) {
        spdlog::error("Opening {} failed", path);
        return false;
      }
    }
  }
  if (!fs::is_directory(output_dir) || !fs::exists(output_dir)) {
    spdlog::error(
        "The \"output_dir\" {} is not a directory or does not exist",
        output_dir);
    return false;
  }
  if (!raster_name_field_name.empty()) {
    if (raster_name_field_idx = layer->GetLayerDefn()->GetFieldIndex(
            utils::CreateCplString(CPLRecode(
                raster_name_field_name.c_str(), "CP936", "UTF-8")).get());
        raster_name_field_idx == -1) {
      spdlog::error(
          "The \"layer\" does not contain the field name {}",
          raster_name_field_name);
      return false;
    }
  }
  return true;
}

}  // anonymous namespace

namespace rs_toolset {
namespace mosaicking {

GDALDatasetUniquePtr CreateMosaickingRaster(
    OGRLayer* layer,
    const std::vector<std::string>& input_paths,
    const std::string& output_path,
    double reso,
    const std::string& compression,
    double blend_dist,
    const std::vector<int>& rgb_bands_map,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing) {
  GDALDatasetUniquePtr output_dataset;
  if (!CreateMosaickingRasterCheck(
          layer, input_paths, output_path, reso, compression, blend_dist,
          color_balancing, output_dataset)) {
    return nullptr;
  }

  for (auto i(0); i < layer->GetFeatureCount(); ++i) {
    if (input_paths[i].empty()) {
      spdlog::info("Skipping warping since the input path is empty");
      continue;
    }
    spdlog::info(
        "Warping {} by the corresponding geometry to the mosaicking raster",
        input_paths[i]);
    OGRFeatureUniquePtr feature(layer->GetFeature(i));
    OGRGeometryUniquePtr geometry(feature->StealGeometry());
    geometry->assignSpatialReference(layer->GetSpatialRef());
    if (color_balancing && feature->GetFieldAsInteger(1) != -1) {
      color_balancing->WarpByGeometry(
          {feature->GetFieldAsInteger(1)}, {geometry.get()}, output_dataset,
          rgb_bands_map, GRA_Bilinear, blend_dist);
    } else {
      utils::WarpByGeometry(
          {input_paths[i]}, {geometry.get()}, output_dataset, rgb_bands_map,
          GRA_Bilinear, blend_dist);
    }
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1, layer->GetFeatureCount());
  }
  output_dataset.reset(nullptr);
  spdlog::info("Creating a mosaicking raster - done");
  output_dataset.reset(GDALDataset::Open(
      output_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  utils::CreateRasterPyra(output_dataset.get());
  return output_dataset;
}

bool CreateRastersCut(
    OGRLayer* layer,
    const std::vector<std::string>& input_paths,
    const std::string& output_dir,
    const std::string& raster_name_field_name,
    const std::string& compression,
    const std::vector<int>& rgb_bands_map) {
  auto idx(0);
  if (!CreateRastersCutCheck(
          layer, input_paths, output_dir, raster_name_field_name, compression,
          idx)) {
    return false;
  }

  for (auto i(0); i < input_paths.size(); ++i) {
    if (input_paths[i].empty()) {
      spdlog::info("Skipping cutting since the input path is empty");
      continue;
    }
    spdlog::info("Cutting {} by the corresponding geometry", input_paths[i]);
    GDALDatasetUniquePtr 
        source_dataset(GDALDataset::Open(
            input_paths[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
        output_dataset(nullptr);
    OGRFeatureUniquePtr feature(layer->GetFeature(i));
    auto output_path(output_dir + "\\" + utils::CreateCplString(CPLRecode(
        feature->GetFieldAsString(idx), "UTF-8", "CP936")).get());
    if (fs::exists(output_path)) {
      output_dataset.reset(GDALDataset::Open(
          output_path.c_str(), GDAL_OF_RASTER | GDAL_OF_UPDATE));
    } else {
      GDALDriver* driver(nullptr);
      std::unique_ptr<char*, void(*)(char**)> options(nullptr, nullptr);
      if (driver = utils::GetRasterDriverByPath(output_path, compression);
          driver) {
        options = utils::CreateCslStringList(CSLSetNameValue(
            nullptr, "COMPRESS", compression.c_str()));
      } else if (driver = utils::GetRasterDriverByPath(output_path, "");
          !driver) {
        spdlog::error("Finding the driver for {} failed", output_path);
        continue;
      }
      output_dataset.reset(driver->Create(
          output_path.c_str(), source_dataset->GetRasterXSize(),
          source_dataset->GetRasterYSize(), source_dataset->GetRasterCount(),
          GDT_Byte, options.get()));
      double geotrans[6];
      source_dataset->GetGeoTransform(geotrans);
      output_dataset->SetGeoTransform(geotrans);
      output_dataset->SetSpatialRef(source_dataset->GetSpatialRef());
    }
    OGRGeometryUniquePtr geometry(feature->StealGeometry());
    geometry->assignSpatialReference(const_cast<OGRSpatialReference*>(
        source_dataset->GetSpatialRef()));
    utils::WarpByGeometry(
        {source_dataset.get()}, {geometry.get()}, output_dataset,
        rgb_bands_map);
    output_dataset.reset(nullptr);
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1, input_paths.size());
    output_dataset.reset(GDALDataset::Open(
        output_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    utils::CreateRasterPyra(output_dataset.get(), "DEFLATE", true);
  }
  return true;
}

}  // namespace mosaicking
}  // namespace rs_toolset