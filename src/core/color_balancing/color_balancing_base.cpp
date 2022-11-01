#include "color_balancing_base.h"

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include <gdal_priv.h>
#include <gdalwarper.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/utils.hpp>


namespace fs = std::filesystem;

namespace rs_toolset {
namespace color_balancing {

std::vector<std::string> ColorBalancingBase::ExportAllRastersName() {
  if (rasters_path_.empty()) {
    spdlog::info("No rasters' name can be exported");
    return {};
  }
  spdlog::debug("Exporting all rasters' name");
  std::vector<std::string> names;
  for (const auto& path : rasters_path_)
    names.push_back(fs::path(path).filename().string());
  spdlog::debug("Exporting all rasters' name - done");
  return names;
}

bool ColorBalancingBase::CreateRasters(
    const std::vector<int>& idxes,
    const std::string& output_dir) {
  if (!fs::is_directory(output_dir)) {
    spdlog::warn("The output directory {} is not a directory", output_dir);
    return false;
  }
  if (!fs::exists(output_dir)) {
    spdlog::warn("The output directory {} does not exist", output_dir);
    return false;
  }
  std::vector<int> _idxes;
  if (idxes.empty()) {
    for (int i = 0; i < rasters_path_.size(); i++)
      _idxes.push_back(i);
  } else {
    for (const auto& idx : idxes)
      if (idx >= rasters_path_.size()) {
        spdlog::warn(
            "The index {} is greater than or equal to the rasters' count",
            idx);
        return false;
      }
    _idxes = idxes;
  }
  std::string string("Creating color balanced raster(s):\n");
  for (const auto& idx : _idxes)
    string.append(rasters_path_[idx]).append("\n");
  string.append("{} task(s) in total");
  spdlog::info(string, _idxes.size());

  GDALDriver* gtiff_driver(GetGDALDriverManager()->GetDriverByName("GTiff"));
  for (int i = 0; i < _idxes.size(); i++) {
    GDALDatasetUniquePtr source_dataset(GDALDataset::Open(
        rasters_path_[_idxes[i]].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    double geotrans[6];
    source_dataset->GetGeoTransform(geotrans);
    int x_size(source_dataset->GetRasterXSize()),
        y_size(source_dataset->GetRasterYSize());
    std::string output_path(
        output_dir + "/" +
        fs::path(rasters_path_[_idxes[i]]).filename().string());
    GDALDatasetUniquePtr output_dataset(gtiff_driver->Create(
        output_path.c_str(), x_size, y_size, 3, GDT_Byte, nullptr));
    output_dataset->SetGeoTransform(geotrans);
    output_dataset->SetSpatialRef(source_dataset->GetSpatialRef());
    WriteToDataset(_idxes[i], source_dataset.get(), output_dataset);
    spdlog::info("---------- {}/{} - done ----------", i + 1, _idxes.size());
    output_dataset.reset(nullptr);
    output_dataset.reset(GDALDataset::Open(
        output_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    utils::CreateRasterPyra(output_dataset.get());
  }
  spdlog::info("Creating color balanced raster(s) - done");
  return true;
}

bool ColorBalancingBase::WarpByGeometry(
    const std::vector<int>& idxes,
    const std::vector<OGRGeometry*>& geometries,
    GDALDatasetUniquePtr& output_dataset,
    const std::vector<int>& bands_map,
    GDALResampleAlg resample_arg,
    double blend_dist,
    double nodata_value) {
  if (idxes.size() != geometries.size()) return false;
  int x_size(output_dataset->GetRasterXSize());
  int y_size(output_dataset->GetRasterYSize());
  double geotrans[6];
  output_dataset->GetGeoTransform(geotrans);
  int bands_count(
      bands_map.empty() ? output_dataset->GetRasterCount()
      : static_cast<int>(bands_map.size()));
  auto _bands_map(std::make_unique<int[]>(bands_count));
  if (bands_map.empty()) {
    for (int b = 0; b < bands_count; b++)
      _bands_map[b] = b + 1;
  } else {
    for (int b = 0; b < bands_count; b++)
      _bands_map[b] = bands_map[b];
  }
  auto nodata_values(std::make_unique<double[]>(bands_count));
  for (int b = 0; b < bands_count; b++)
    nodata_values[b] = nodata_value;
  std::unique_ptr<void, void(*)(void*)> trans_arg(
      nullptr, [](void* p) { GDALDestroyGenImgProjTransformer(p); });
  std::unique_ptr<GDALWarpOptions> warp_options(GDALCreateWarpOptions());
  warp_options->eResampleAlg = resample_arg;
  warp_options->eWorkingDataType = GDT_Byte;
  warp_options->nBandCount = 3;
  warp_options->padfDstNoDataReal = nodata_values.get();
  warp_options->padfSrcNoDataReal = nodata_values.get();
  warp_options->panDstBands = _bands_map.get();
  warp_options->panSrcBands = _bands_map.get();
  warp_options->pfnTransformer = GDALGenImgProjTransform;
  warp_options->dfWarpMemoryLimit = 8.0 * 1024 * 1024 * 1024;
  warp_options->papszWarpOptions = CSLSetNameValue(
      warp_options->papszWarpOptions, "SKIP_NOSOURCE", "YES");
  warp_options->papszWarpOptions = CSLSetNameValue(
      warp_options->papszWarpOptions, "CUTLINE_ALL_TOUCHED", "TRUE");
  GDALWarpOperation warp_operation;
  for (int i = 0; i < idxes.size(); i++) {
    auto mem_driver(GetGDALDriverManager()->GetDriverByName("MEM"));
    GDALDatasetUniquePtr
        source_dataset(GDALDataset::Open(
            rasters_path_[idxes[i]].c_str(),
            GDAL_OF_RASTER | GDAL_OF_READONLY)),
        color_balanced_dataset(mem_driver->Create(
            "", x_size, y_size, 3, GDT_Byte, nullptr));
    color_balanced_dataset->SetGeoTransform(geotrans);
    color_balanced_dataset->SetSpatialRef(output_dataset->GetSpatialRef());

    // Warp the source dataset to the color balanced dataset with cutline without blend
    warp_options->dfCutlineBlendDist = 0.0;
    warp_options->hDstDS = color_balanced_dataset.get();
    trans_arg.reset(GDALCreateGenImgProjTransformer2(
        source_dataset.get(), color_balanced_dataset.get(), nullptr));
    OGRGeometryUniquePtr cutline(nullptr);
    if (geometries[i]) {
      double source_geotrans[6], inv_source_geotrans[6];
      source_dataset->GetGeoTransform(source_geotrans);
      GDALInvGeoTransform(source_geotrans, inv_source_geotrans);
      cutline.reset(geometries[i]->clone());
      if (cutline->getSpatialReference())
        cutline->transformTo(const_cast<OGRSpatialReference*>(
            source_dataset->GetSpatialRef()));
      cutline = utils::ApplyGeotransOnGeometry(
          cutline.get(), inv_source_geotrans);
      warp_options->hCutline = cutline.get();
    }
    warp_options->hSrcDS = source_dataset.get();
    warp_options->pTransformerArg = trans_arg.get();
    warp_operation.Initialize(warp_options.get());
    warp_operation.ChunkAndWarpMulti(0, 0, x_size, y_size);
    WriteToDataset(
        idxes[i], color_balanced_dataset.get(), color_balanced_dataset);

    // Warp the color balanced dataset to the output dataset without cutline with blend
    warp_options->hCutline = nullptr;
    warp_options->dfCutlineBlendDist = blend_dist;
    warp_options->hDstDS = output_dataset.get();
    trans_arg.reset(GDALCreateGenImgProjTransformer2(
        color_balanced_dataset.get(), output_dataset.get(), nullptr));
    warp_options->hSrcDS = color_balanced_dataset.get();
    warp_options->pTransformerArg = trans_arg.get();
    warp_operation.Initialize(warp_options.get());
    warp_operation.ChunkAndWarpMulti(0, 0, x_size, y_size);
  }
  return true;
}

std::shared_ptr<ColorBalancingInterface> ColorBalancingInterface::Create(
    const std::string& method,
    int block_size) {
  if (method == "global-to-local") {
    return GlobalToLocal::Create(block_size);
  } else {
    spdlog::warn("The method {} is not supported", method);
    return nullptr;
  }
}

std::shared_ptr<ColorBalancingInterface> 
ColorBalancingInterface::CreateFromArgusJson(
    const std::string& json_path,
    int block_size) {
  std::ifstream i(json_path);
  if (!i) {
    spdlog::warn("Opening {} failed", json_path);
    return nullptr;
  }

  nlohmann::json json;
  i >> json;
  try {
    auto method(json.at("method").get<std::string>());
    if (method == "global-to-local"){
      auto color_balancing(GlobalToLocal::Create(block_size));
      color_balancing->ImportArgusJson(json_path);
      return color_balancing;
    } else {
      spdlog::warn("The method {} is not supported", method);
      return nullptr;
    }
  } catch (...) {
    spdlog::warn("The json file {} has been broken", json_path);
    return nullptr;
  }
}

} // namespace color_balancing
} // namespace rs_toolset