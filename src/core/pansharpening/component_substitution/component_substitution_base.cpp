#include "component_substitution_base.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <opencv2/opencv.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretching.h>
#include <rs-toolset/utils.hpp>

namespace rs_toolset {
namespace pansharpening {

bool ComponentSubstitutionBase::Run(
    const std::string& pan_path,
    const std::string& ms_path,
    const std::string& output_path,
    bool use_rpc,
    bool use_stretching,
    const std::vector<int>& bands_map) {
  std::string string(
      "Running a pansharpening task with\n"
      " - PAN path: {}\n"
      " - MS path: {}\n"
      " - Use RPC: {}\n"
      " - Use stretching: {}\n"
      " - Bands' map: ");
  if (bands_map.empty()) {
    string.append("All bands");
  } else {
    for (const auto& idx : bands_map)
      string.append(std::to_string(idx)).append(",");
    string.pop_back();
  }
  spdlog::info(string, pan_path, ms_path, use_rpc, use_stretching);
  GDALDatasetUniquePtr
      pan_dataset(GDALDataset::Open(
          pan_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      ms_dataset(GDALDataset::Open(
          ms_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!pan_dataset) {
    spdlog::error("Opening {} failed", pan_path);
    return false;
  }
  if (!ms_dataset) {
    spdlog::error("Opening {} failed", ms_path);
    return false;
  }
  auto data_type(pan_dataset->GetRasterBand(1)->GetRasterDataType());
  if (auto ms_data_type(ms_dataset->GetRasterBand(1)->GetRasterDataType());
      data_type != ms_data_type) {
    spdlog::error(
        "The PAN raster's data type {} must be the same with "
        "the MS raster's {}", data_type, ms_data_type);
    return false;
  }
  if (data_type != GDT_Byte && data_type != GDT_UInt16) {
    spdlog::error(
        "The data type {} must be unsigned 8-bit or unsigned 16-bit",
        data_type);
    return false;
  }
  auto driver(utils::GetRasterDriverByPath(output_path));
  if (!driver) {
    spdlog::error("Finding the driver for {} failed", output_path);
    return false;
  }

  auto bands_count(ms_dataset->GetRasterCount()),
      pan_x_size(pan_dataset->GetRasterXSize()),
      pan_y_size(pan_dataset->GetRasterYSize()),
      ms_x_size(ms_dataset->GetRasterXSize()),
      ms_y_size(ms_dataset->GetRasterYSize());
  auto statistics(CreateStatistics(bands_count, ms_x_size, ms_y_size));
  utils::RPCTransPtr
      pan_trans_arg(nullptr, nullptr),
      ms_trans_arg(nullptr, nullptr);
  if (use_rpc) {
    pan_trans_arg = utils::CreateRPCTrans(pan_path),
    ms_trans_arg = utils::CreateRPCTrans(ms_path);
  }
  Data data;
  int ms_range[4], pan_range[4];

  // Create weights
  spdlog::info("Creating weights");
  if (need_downsample_info_) {
    auto ms_block_cols_count(static_cast<int>(ceil(static_cast<double>(
            ms_x_size) / block_size_))),
        ms_block_rows_count(static_cast<int>(ceil(static_cast<double>(
            ms_y_size) / block_size_))),
        last_ms_block_x_size(
            ms_x_size - block_size_ * (ms_block_cols_count - 1)),
        last_ms_block_y_size(
            ms_y_size - block_size_ * (ms_block_rows_count - 1));
    spdlog::info(
        "Dividing the {}x{} MS raster into {} {}x{} blocks", ms_x_size,
        ms_y_size, ms_block_cols_count * ms_block_rows_count, block_size_,
        block_size_);
    for (int i(0); i < ms_block_cols_count * ms_block_rows_count; ++i) {
      if (use_rpc) {
        CreateRangesForRPC(
            pan_x_size, pan_y_size, i, ms_block_cols_count, ms_block_rows_count,
            block_size_, block_size_, last_ms_block_x_size,
            last_ms_block_y_size, pan_trans_arg, ms_trans_arg, pan_range,
            ms_range);
        data = CreateDownsampledDataWithRPC(
            pan_dataset.get(), ms_dataset.get(), pan_range, ms_range,
            pan_trans_arg, ms_trans_arg);
      } else {
        utils::CreateRange(
            i, ms_block_cols_count, ms_block_rows_count, block_size_,
            block_size_, last_ms_block_x_size, last_ms_block_y_size, ms_range);
        data = CreateDownsampledDataWithGeotrans(
            pan_dataset.get(), ms_dataset.get(), ms_range);
      }
      UpdateDownsampleInfo(data, statistics);
      spdlog::info(
          "---------- {}/{} - done ----------", i + 1,
          ms_block_cols_count * ms_block_rows_count);
    }
  }
  auto weights(CreateWeights(statistics));
  spdlog::info("Creating weights - done");

  // Create injection gains
  spdlog::info("Creating injection gains");
  auto pan_block_cols_count(static_cast<int>(ceil(static_cast<double>(
          pan_x_size) / block_size_))),
      pan_block_rows_count(static_cast<int>(ceil(static_cast<double>(
          pan_y_size) / block_size_))),
      last_pan_block_x_size(
          pan_x_size - block_size_ * (pan_block_cols_count - 1)),
      last_pan_block_y_size(
          pan_y_size - block_size_ * (pan_block_rows_count - 1));
  spdlog::info(
      "Dividing the {}x{} PAN raster into {} {}x{} blocks", pan_x_size,
      pan_y_size, pan_block_cols_count * pan_block_rows_count, block_size_,
      block_size_);
  for (int i(0); i < pan_block_cols_count * pan_block_rows_count; ++i) {
    if (use_rpc) {
      CreateRangesForRPC(
          ms_x_size, ms_y_size, i, pan_block_cols_count, pan_block_rows_count,
          block_size_, block_size_, last_pan_block_x_size,
          last_pan_block_y_size, ms_trans_arg, pan_trans_arg, ms_range,
          pan_range);
      data = CreateUpsampledDataWithRPC(
          pan_dataset.get(), ms_dataset.get(), pan_range, ms_range,
          pan_trans_arg, ms_trans_arg);
    } else {
      utils::CreateRange(
          i, pan_block_cols_count, pan_block_rows_count, block_size_,
          block_size_, last_pan_block_x_size, last_pan_block_y_size, pan_range);
      data = CreateUpsampledDataWithGeotrans(
          pan_dataset.get(), ms_dataset.get(), pan_range);
    }
    UpdateUpsampleInfo(data, weights, statistics);
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1,
        pan_block_cols_count * pan_block_rows_count);
  }
  auto injection_gains(CreateInjectionGains(statistics));
  DestroyStatistics(statistics);
  spdlog::info("Creating injection gains - done");

  // Create pansharpened mats
  spdlog::info("Creating pansharpened mats");
  std::vector<cv::Mat> pansharpened_mats(bands_count);
  for (auto& mat : pansharpened_mats)
    mat = cv::Mat(pan_y_size, pan_x_size, CV_16UC1);
  std::shared_ptr<stretching::PercentClip> stretching(nullptr);
  if (use_stretching)
    stretching = stretching::PercentClip::Create();
  spdlog::info(
      "Dividing {}x{} pansharpened raster into {} {}x{} blocks", pan_x_size,
      pan_y_size, pan_block_cols_count * pan_block_rows_count, block_size_,
      block_size_);
  for (int i(0); i < pan_block_cols_count * pan_block_rows_count; ++i) {
    if (use_rpc) {
      CreateRangesForRPC(
          ms_x_size, ms_y_size, i, pan_block_cols_count, pan_block_rows_count,
          block_size_, block_size_, last_pan_block_x_size,
          last_pan_block_y_size, ms_trans_arg, pan_trans_arg, ms_range,
          pan_range);
      data = CreateUpsampledDataWithRPC(
          pan_dataset.get(), ms_dataset.get(), pan_range, ms_range,
          pan_trans_arg, ms_trans_arg);
    } else {
      utils::CreateRange(
          i, pan_block_cols_count, pan_block_rows_count, block_size_,
          block_size_, last_pan_block_x_size, last_pan_block_y_size, pan_range);
      data = CreateUpsampledDataWithGeotrans(
          pan_dataset.get(), ms_dataset.get(), pan_range);
    }
    auto delta_mats(CreateDeltaMats(data, weights));
    cv::Rect rect(pan_range[0], pan_range[1], pan_range[2], pan_range[3]);
    std::vector<cv::Mat> mats(bands_count);
#pragma omp parallel for schedule(static, bands_count)
    for (int b(0); b < bands_count; ++b) {
      data.ms_mats[b].copyTo(pansharpened_mats[b](rect));
      pansharpened_mats[b](rect) += injection_gains[b] * delta_mats[b];
      mats[b] = pansharpened_mats[b](rect);
    }
    if (use_stretching)
      stretching->AccumulateStat(mats);
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1,
        pan_block_cols_count * pan_block_rows_count);
  }
  if (use_stretching) {
    auto sink(spdlog::default_logger());
    spdlog::set_default_logger(spdlog::stdout_color_mt("Stretching"));
    stretching->CreateLutMats();
    stretching->Run(pansharpened_mats);
    stretching->Clear();
    data_type = GDT_Byte;
    spdlog::set_default_logger(sink);
  }
  spdlog::info("Creating pansharpened mats - done");

  // Write the pansharpened mats into the given path
  spdlog::info("Writing the pansharpened mats into {}", output_path);
  if (!bands_map.empty())
    bands_count = static_cast<int>(bands_map.size());
  GDALDatasetUniquePtr output_dataset(driver->Create(
      output_path.c_str(), pan_x_size, pan_y_size, bands_count, data_type,
      nullptr));
  if (!output_dataset) {
    spdlog::error("Creating {} failed", output_path);
    return false;
  }
  if (use_rpc) {
    output_dataset->SetMetadata(pan_dataset->GetMetadata("RPC"), "RPC");
  } else {
    double geotrans[6];
    pan_dataset->GetGeoTransform(geotrans);
    output_dataset->SetGeoTransform(geotrans);
    output_dataset->SetSpatialRef(pan_dataset->GetSpatialRef());
  }
  auto bytes_count(GDALGetDataTypeSizeBytes(data_type));
#pragma omp parallel for schedule(static, bands_count)
  for (int b(0); b < bands_count; ++b) {
    output_dataset->GetRasterBand(b + 1)->RasterIO(
        GF_Write, 0, 0, pan_x_size, pan_y_size, 
        pansharpened_mats[bands_map.empty() ? b : bands_map[b] - 1].data,
        pan_x_size, pan_y_size, data_type, bytes_count,
        bytes_count * pan_x_size);
  }
  spdlog::info("Writing the pansharpened mats into - done");
  spdlog::info("Running a pansharpening task - done");
  return true;
}

std::vector<cv::Mat> ComponentSubstitutionBase::CreateDeltaMats(
    const Data& data,
    const std::vector<double>& weights) {
  // Minus the synthetic low resolution PAN mat
  int bands_count(static_cast<int>(data.ms_mats.size()));
  cv::Mat delta_mat;
  data.ms_mats[0].convertTo(delta_mat, CV_16SC1, -weights[0]);
  for (int i(1); i < bands_count; ++i)
    delta_mat -= weights[i] * data.ms_mats[i];
  if (weights.size() == bands_count + 1) {
    cv::Mat constant_mat;
    cv::threshold(
        data.pan_mat, constant_mat, 0, weights.back(), cv::THRESH_BINARY);
    delta_mat -= constant_mat;
  }

  // Add the histogram matched PAN mat
  delta_mat += utils::MapMatWithLut(data.pan_mat, hist_matching_mat_);
  return std::vector<cv::Mat>(bands_count, delta_mat);
}

}  // namespace pansharpening
}  // namespace rs_toolset