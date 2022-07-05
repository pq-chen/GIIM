#include "component_substitution_base.h"

#include <cmath>

#include <string>
#include <vector>

#include <gdal_priv.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretch.h>
#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace pansharpening {

bool ComponentSubstitutionBase::Run(
    const std::string& pan_path,
    const std::string& ms_path,
    const std::string& pansharpened_path,
    bool use_rpc,
    bool use_stretch,
    const std::vector<int>& pansharpened_bands_map) {
  std::string string(
      "Running a pansharpening task from\n - PAN path: {}\n - MS path: {}\n"
      " - Pansharpened path: {}\n - Use RPC: {}\n - Use stretch: {}\n"
      " - Pansharpened bands' map: ");
  for (const auto& idx : pansharpened_bands_map)
    string.append(std::to_string(idx)).append(",");
  string.pop_back();
  spdlog::info(
      string, pan_path, ms_path, pansharpened_path, use_rpc, use_stretch);
  GDALDatasetUniquePtr 
      pan_dataset(GDALDataset::Open(
          pan_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      ms_dataset(GDALDataset::Open(
          ms_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!pan_dataset) {
    spdlog::warn("Opening {} failed", pan_path);
    return false;
  }
  if (!ms_dataset) {
    spdlog::warn("Opening {} failed", ms_path);
    return false;
  }
  auto pan_dataset_type(pan_dataset->GetRasterBand(1)->GetRasterDataType()),
      ms_dataset_type(ms_dataset->GetRasterBand(1)->GetRasterDataType());
  if (pan_dataset_type != ms_dataset_type) return false;
  if (pan_dataset_type != GDT_Byte && pan_dataset_type != GDT_UInt16) {
    spdlog::warn(
        "The PAN raster's data type must be unsigned 8-bit or "
        "unsigned 16-bit");
    return false;
  }
  if (ms_dataset_type != pan_dataset_type) {
    spdlog::warn(
        "The MS raster's data type must be the same as the PAN raster's");
    return false;
  }

  int bands_count(ms_dataset->GetRasterCount()),
      bytes_count(GDALGetDataTypeSizeBytes(pan_dataset_type)),
      pan_x_size(pan_dataset->GetRasterXSize()),
      pan_y_size(pan_dataset->GetRasterYSize()),
      ms_x_size(ms_dataset->GetRasterXSize()),
      ms_y_size(ms_dataset->GetRasterYSize()),
      pan_block_cols_count(static_cast<int>(ceil(static_cast<double>(
          pan_x_size) / block_size_))),
      pan_block_rows_count(static_cast<int>(ceil(static_cast<double>(
          pan_y_size) / block_size_))),
      last_pan_block_x_size(
          pan_x_size - (pan_block_cols_count - 1) * block_size_),
      last_pan_block_y_size(
          pan_y_size - (pan_block_rows_count - 1) * block_size_);
  RPCTransPtr pan_trans_arg(nullptr, nullptr), ms_trans_arg(nullptr, nullptr);
  if (use_rpc) {
    spdlog::info(
        "Using the RPC information to create the mapping relation "
        "between the PAN raster and the MS raster before orthorectification");
    pan_trans_arg = utils::CreateRPCTrans(pan_path, true),
    ms_trans_arg = utils::CreateRPCTrans(ms_path, false);
  } else {
    spdlog::info(
        "Using the geotransform to create the mapping relation "
        "between the PAN raster and the MS raster after orthorectification");
  }
  void* s(CreateStatistic(bands_count, ms_x_size, ms_y_size));

  // Update the downsample information
  if (need_downsample_info_) {
    spdlog::info(
        "Updating the downsample information in the statistic struct");
    int ms_block_cols_count(static_cast<int>(ceil(static_cast<double>(
            ms_x_size) / block_size_))),
        ms_block_rows_count(static_cast<int>(ceil(static_cast<double>(
            ms_y_size) / block_size_))),
        last_ms_block_x_size(
            ms_x_size - (ms_block_cols_count - 1) * block_size_),
        last_ms_block_y_size(
            ms_y_size - (ms_block_rows_count - 1) * block_size_);
    RPCTransPtr 
        _pan_trans_arg(nullptr, nullptr),
        _ms_trans_arg(nullptr, nullptr);
    if (use_rpc) {
      _pan_trans_arg = utils::CreateRPCTrans(pan_path, false),
      _ms_trans_arg = utils::CreateRPCTrans(ms_path, true);
    }
    spdlog::info(
        "Dividing the {}x{} size MS raster into {} {}x{} blocks",
        ms_x_size, ms_y_size, ms_block_cols_count * ms_block_rows_count,
        block_size_, block_size_);
    for (int i = 0; i < ms_block_cols_count * ms_block_rows_count; i++) {
      Data data;
      int ms_range[4];
      if (use_rpc) {
        int pan_range[4];
        CreateRangesForRPC(
            i, ms_block_cols_count, ms_block_rows_count, block_size_,
            block_size_, last_ms_block_x_size, last_ms_block_y_size,
            pan_x_size, pan_y_size, _ms_trans_arg, _pan_trans_arg, ms_range,
            pan_range);
        data = CreateDownsampledDataWithRPC(
            pan_path, ms_path, pan_range, ms_range, _pan_trans_arg, 
            _ms_trans_arg);
      } else {
        utils::CreateRange(
            i, ms_block_cols_count, ms_block_rows_count, block_size_,
            block_size_, last_ms_block_x_size, last_ms_block_y_size, ms_range);
        data = CreateDownsampledDataWithGeotrans(pan_path, ms_path, ms_range);
      }
      UpdateDownsampleInfo(data, s);
      spdlog::info(
          "---------- {}/{} - done ----------", i + 1,
          ms_block_cols_count * ms_block_rows_count);
    }
    spdlog::info(
        "Updating the downsample information in the statistic struct - done");
  }
  std::vector<double> weights(CreateWeights(s));

  // Update the upsample information
  spdlog::info("Updating the statistic struct with the upsample information");
  spdlog::info(
      "Dividing the {}x{} PAN raster into {} {}x{} blocks",
      pan_x_size, pan_y_size, pan_block_cols_count * pan_block_rows_count,
      block_size_, block_size_);
  for (int i = 0; i < pan_block_cols_count * pan_block_rows_count; i++) {
    Data data;
    int pan_range[4];
    if (use_rpc) {
      int ms_range[4];
      CreateRangesForRPC(
          i, pan_block_cols_count, pan_block_rows_count, block_size_,
          block_size_, last_pan_block_x_size, last_pan_block_y_size, ms_x_size,
          ms_y_size, pan_trans_arg, ms_trans_arg, pan_range, ms_range);
      data = CreateUpsampledDataWithRPC(
          pan_path, ms_path, pan_range, ms_range, pan_trans_arg, ms_trans_arg);
    } else {
      utils::CreateRange(
          i, pan_block_cols_count, pan_block_rows_count, block_size_,
          block_size_, last_pan_block_x_size, last_pan_block_y_size,
          pan_range);
      data = CreateUpsampledDataWithGeotrans(pan_path, ms_path, pan_range);
    }
    UpdateUpsampleInfo(data, weights, s);
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1,
        pan_block_cols_count * pan_block_rows_count);
  }
  spdlog::info(
      "Updating the statistic struct with the upsample information - done");
  std::vector<double> injection_gains(CreateInjectionGains(s));

  // Create pansharpened mats
  spdlog::info("Creating pansharpened mats");
  std::vector<cv::Mat> pansharpened_mats(bands_count);
  for (auto& pansharpened_mat : pansharpened_mats)
    pansharpened_mat = cv::Mat(pan_y_size, pan_x_size, CV_16UC1);
  std::shared_ptr<stretch::PercentClip> stretch(nullptr);
  if (use_stretch) {
    spdlog::info("Stretching on the pansharpened mats");
    stretch = stretch::PercentClip::Create();
  }
  spdlog::info(
      "Dividing {}x{} pansharpened raster into {} {}x{} blocks",
      pan_x_size, pan_y_size, pan_block_cols_count* pan_block_rows_count,
      block_size_, block_size_);
  for (int i = 0; i < pan_block_cols_count * pan_block_rows_count; i++) {
    Data data;
    int pan_range[4];
    if (use_rpc) {
      int ms_range[4];
      CreateRangesForRPC(
          i, pan_block_cols_count, pan_block_rows_count, block_size_,
          block_size_, last_pan_block_x_size, last_pan_block_y_size, ms_x_size,
          ms_y_size, pan_trans_arg, ms_trans_arg, pan_range, ms_range);
      data = CreateUpsampledDataWithRPC(
          pan_path, ms_path, pan_range, ms_range, pan_trans_arg, ms_trans_arg);
    } else {
      utils::CreateRange(
          i, pan_block_cols_count, pan_block_rows_count, block_size_,
          block_size_, last_pan_block_x_size, last_pan_block_y_size,
          pan_range);
      data = CreateUpsampledDataWithGeotrans(pan_path, ms_path, pan_range);
    }
    std::vector<cv::Mat> delta_mats(CreateDeltaMats(data, weights, s));
    cv::Rect rect(pan_range[0], pan_range[1], pan_range[2], pan_range[3]);
#pragma omp parallel for schedule(static, bands_count)
    for (int b = 0; b < bands_count; b++) {
      data.mats[b].copyTo(pansharpened_mats[b](rect));
      pansharpened_mats[b](rect) += injection_gains[b] * delta_mats[b];
      if (use_stretch)
        stretch->AddStatForSingleBlock(pansharpened_mats[b](rect), b);
    }
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1,
        pan_block_cols_count * pan_block_rows_count);
  }
  spdlog::info("Creating pansharpened mats - done");
  DestroyStatistic(s);

  // Write the pansharpened mats into the given path
  GDALDataType pansharpened_dataset_type(pan_dataset_type);
  int pansharpened_bytes_count(bytes_count);
  if (use_stretch) {
    stretch->Run(pansharpened_mats);
    stretch->ClearStat();
    pansharpened_dataset_type = GDT_Byte;
    pansharpened_bytes_count = 1;
  }
  spdlog::debug("Writing the pansharpened mats into {}", pansharpened_path);
  int output_bands_count(static_cast<int>(pansharpened_bands_map.size()));
  if (!output_bands_count)
    output_bands_count = bands_count;
  GDALDriver* gtiff_driver(GetGDALDriverManager()->GetDriverByName("GTiff"));
  GDALDatasetUniquePtr pansharpened_dataset(gtiff_driver->Create(
      pansharpened_path.c_str(), pan_x_size, pan_y_size, output_bands_count,
      pansharpened_dataset_type, nullptr));
  if (use_rpc) {
    pansharpened_dataset->SetMetadata(pan_dataset->GetMetadata("RPC"), "RPC");
  } else {
    double geotrans[6];
    pan_dataset->GetGeoTransform(geotrans);
    pansharpened_dataset->SetGeoTransform(geotrans);
    pansharpened_dataset->SetSpatialRef(pan_dataset->GetSpatialRef());
  }
#pragma omp parallel for schedule(static, output_bands_count)
  for (int b = 0; b < output_bands_count; b++) {
    int band_idx(
        pansharpened_bands_map.size() ? pansharpened_bands_map[b] - 1 : b);
    pansharpened_dataset->GetRasterBand(b + 1)->RasterIO(
        GF_Write, 0, 0, pan_x_size, pan_y_size, 
        pansharpened_mats[band_idx].data, pan_x_size, pan_y_size,
        pansharpened_dataset_type, pansharpened_bytes_count,
        pansharpened_bytes_count * pan_x_size);
  }
  spdlog::info(
      "Writing the pansharpened mats into {} - done", pansharpened_path);
  spdlog::info("Running a pansharpening task - done");
  return true;
}

} // pansharpening
} // rs_toolset