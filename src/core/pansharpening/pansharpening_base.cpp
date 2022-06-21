#include "pansharpening_base.h"

#include <cmath>
#include <cstdint>
#include <cstring>

#include <memory>
#include <string>
#include <vector>

#include <gdal_alg.h>
#include <gdal_priv.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace pansharpening {

void PansharpeningBase::CreateRangesForRPC(
    int target_block_idx,
    int target_block_cols_count,
    int target_block_rows_count,
    int target_block_x_size,
    int target_block_y_size,
    int last_target_block_x_size,
    int last_target_block_y_size,
    int source_x_size,
    int source_y_size,
    const RPCTransPtr& forward_trans_arg,
    const RPCTransPtr& backward_trans_arg,
    int (&target_range)[4],
    int (&source_range)[4]) {
  utils::CreateRange(
      target_block_idx, target_block_cols_count, target_block_rows_count,
      target_block_x_size, target_block_y_size, last_target_block_x_size,
      last_target_block_y_size, target_range);
  double 
      ms_xs[4]{
          static_cast<double>(target_range[0]),
          static_cast<double>(target_range[0] + target_range[2]),
          static_cast<double>(target_range[0]),
          static_cast<double>(target_range[0] + target_range[2]) },
      ms_ys[4]{
          static_cast<double>(target_range[1]),
          static_cast<double>(target_range[1]),
          static_cast<double>(target_range[1] + target_range[3]),               
          static_cast<double>(target_range[1] + target_range[3]) },
      ms_zs[4] { 0.0, 0.0, 0.0, 0.0 };                                          
  int successes[4];                                                             
  GDALRPCTransform(               
      forward_trans_arg.get(), true, 4, ms_xs, ms_ys, ms_zs, successes);
  GDALRPCTransform(
      backward_trans_arg.get(), true, 4, ms_xs, ms_ys, ms_zs, successes);
  source_range[0] = 
      static_cast<int>(floor(fmax(fmin(ms_xs[0], ms_xs[2]), 0.)));
  source_range[1] = 
      static_cast<int>(floor(fmax(fmin(ms_ys[0], ms_ys[1]), 0.)));
  source_range[2] = static_cast<int>(ceil(fmin(fmax(ms_xs[1], ms_xs[3]), 
      static_cast<double>(source_x_size)))) - source_range[0];
  source_range[3] = static_cast<int>(ceil(fmin(fmax(ms_ys[2], ms_ys[3]),
      static_cast<double>(source_y_size)))) - source_range[1];
}

PansharpeningBase::Data PansharpeningBase::CreateUpsampledDataWithRPC(
    const std::string& pan_path,
    const std::string& ms_path,
    int* pan_range,
    int* ms_range,
    const RPCTransPtr& pan_trans_arg,
    const RPCTransPtr& ms_trans_arg) {
  spdlog::debug("Creating the upsampled data with the rpc information");
  GDALDatasetUniquePtr 
      pan_dataset(GDALDataset::Open(
          pan_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      ms_dataset(GDALDataset::Open(
          ms_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  Data data;
  data.mat = utils::CreateMatFromDataset(pan_dataset.get(), pan_range);
  data.mats = CreateResampledMats(
      utils::CreateMatFromDataset(ms_dataset.get(), ms_range), data.mat,
      ms_range, pan_range, pan_trans_arg, ms_trans_arg);
  MaskData(data);
  spdlog::info("Creating the upsampled data with the rpc information - done");
  return data;
}

PansharpeningBase::Data PansharpeningBase::CreateDownsampledDataWithRPC(
    const std::string& pan_path,
    const std::string& ms_path,
    int* pan_range,
    int* ms_range,
    const RPCTransPtr& pan_trans_arg,
    const RPCTransPtr& ms_trans_arg) {
  spdlog::debug("Creating the downsampled data with the rpc information");
  GDALDatasetUniquePtr 
      pan_dataset(GDALDataset::Open(
          pan_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      ms_dataset(GDALDataset::Open(
          ms_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  Data data;
  cv::Mat ms_mat(utils::CreateMatFromDataset(ms_dataset.get(), ms_range));
  cv::split(ms_mat, data.mats);
  data.mat = CreateResampledMats(
      utils::CreateMatFromDataset(pan_dataset.get(), pan_range), ms_mat,
      pan_range, ms_range, ms_trans_arg, pan_trans_arg)[0];
  MaskData(data);
  spdlog::info("Creating the downsampled data with the rpc information - done");
  return data;
}

PansharpeningBase::Data PansharpeningBase::CreateUpsampledDataWithGeotrans(
    const std::string& pan_path,
    const std::string& ms_path,
    int* target_range) {
  spdlog::debug("Creating the upsampled data with the geotransform");
  GDALDatasetUniquePtr 
      pan_dataset(GDALDataset::Open(
          pan_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      ms_dataset(GDALDataset::Open(
          ms_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  GDALDriver* mem_driver(GetGDALDriverManager()->GetDriverByName("MEM"));
  Data data;
  double geotrans[6];
  pan_dataset->GetGeoTransform(geotrans);
  geotrans[0] += target_range[0] * geotrans[1];
  geotrans[3] += target_range[1] * geotrans[5];
  GDALDatasetUniquePtr upsampled_ms_dataset(mem_driver->Create(
      "", target_range[2], target_range[3], ms_dataset->GetRasterCount(),
      ms_dataset->GetRasterBand(1)->GetRasterDataType(), nullptr));
  upsampled_ms_dataset->SetGeoTransform(geotrans);
  upsampled_ms_dataset->SetSpatialRef(ms_dataset->GetSpatialRef());
  for (int b = 0; b < ms_dataset->GetRasterCount(); b++)
    upsampled_ms_dataset->GetRasterBand(b + 1)->SetNoDataValue(
        ms_dataset->GetRasterBand(b + 1)->GetNoDataValue());
  utils::WarpByGeometry(
      { ms_dataset.get() }, upsampled_ms_dataset.get(), nullptr);
  data.mat = utils::CreateMatFromDataset(pan_dataset.get(), target_range);
  cv::split(utils::CreateMatFromDataset(upsampled_ms_dataset.get()), data.mats);
  MaskData(data);
  spdlog::info("Creating the upsampled data with the geotransform - done");
  return data;
}

PansharpeningBase::Data PansharpeningBase::CreateDownsampledDataWithGeotrans(
    const std::string& pan_path,
    const std::string& ms_path,
    int* target_range) {
  spdlog::debug("Creating the downsampled data with the geotransform");
  GDALDatasetUniquePtr 
      pan_dataset(GDALDataset::Open(
          pan_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)),
      ms_dataset(GDALDataset::Open(
          ms_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  GDALDriver* mem_driver(GetGDALDriverManager()->GetDriverByName("MEM"));
  Data data;
  double geotrans[6];
  ms_dataset->GetGeoTransform(geotrans);
  geotrans[0] += target_range[0] * geotrans[1];
  geotrans[3] += target_range[1] * geotrans[5];
  GDALDatasetUniquePtr downsampled_pan_dataset(mem_driver->Create(
      "", target_range[2], target_range[3], pan_dataset->GetRasterCount(),
      pan_dataset->GetRasterBand(1)->GetRasterDataType(), nullptr));
  downsampled_pan_dataset->SetGeoTransform(geotrans);
  downsampled_pan_dataset->SetSpatialRef(pan_dataset->GetSpatialRef());
  for (int b = 0; b < pan_dataset->GetRasterCount(); b++)
    downsampled_pan_dataset->GetRasterBand(b + 1)->SetNoDataValue(
        pan_dataset->GetRasterBand(b + 1)->GetNoDataValue());
  utils::WarpByGeometry(
      { pan_dataset.get() }, downsampled_pan_dataset.get(), nullptr);
  cv::split(
      utils::CreateMatFromDataset(ms_dataset.get(), target_range), data.mats);
  data.mat = utils::CreateMatFromDataset(downsampled_pan_dataset.get());
  MaskData(data);
  spdlog::info("Creating the downsampled data with the geotransform - done");
  return data;
}

std::vector<cv::Mat> PansharpeningBase::CreateResampledMats(
    const cv::Mat& source_mat,
    const cv::Mat& target_mat,
    int* source_range,
    int* target_range,
    const RPCTransPtr& forward_trans_arg,
    const RPCTransPtr& backward_trans_arg) {
  int bands_count(source_mat.channels());
  uint64_t size(static_cast<uint64_t>(target_mat.rows) * target_mat.cols);
  std::unique_ptr<void, void(*)(void*)> resampled_source_data(
      nullptr, [](void* p) { delete[] p; });
  if (source_mat.depth() == CV_8U) {
    resampled_source_data.reset(static_cast<void*>(
        new uint8_t[bands_count * size]));
  } else {
    resampled_source_data.reset(static_cast<void*>(
        new uint16_t[bands_count * size]));
  }
#pragma omp parallel for schedule(dynamic)
  for (int row = 0; row < target_range[3]; row++) {
    auto xs(std::make_unique<double[]>(target_range[2])),
        ys(std::make_unique<double[]>(target_range[2])),
        zs(std::make_unique<double[]>(target_range[2]));
    auto successes(std::make_unique<int[]>(target_range[2]));
    memset(zs.get(), 0, sizeof(double) * target_range[2]);
    for (int col = 0; col < target_range[2]; col++) {
      xs[col] = col + target_range[0];
      ys[col] = row + target_range[1];
    }
    GDALRPCTransform(
        forward_trans_arg.get(), true, target_range[2], xs.get(), ys.get(),
        zs.get(), successes.get());
    GDALRPCTransform(
        backward_trans_arg.get(), true, target_range[2], xs.get(), ys.get(),
        zs.get(), successes.get());
    for (int col = 0; col < target_range[2]; col++) {
      xs[col] -= source_range[0];
      ys[col] -= source_range[1];
    }
    if (source_mat.depth() == CV_8U) {
      utils::GetInterpolation<uint8_t>(
          source_mat, target_range[2], xs.get(), ys.get(),
          static_cast<uint8_t*>(resampled_source_data.get()) + 
          static_cast<uint64_t>(bands_count) * row * target_range[2]);
    } else {
      utils::GetInterpolation<uint16_t>(
          source_mat, target_range[2], xs.get(), ys.get(),
          static_cast<uint16_t*>(resampled_source_data.get()) + 
          static_cast<uint64_t>(bands_count) * row * target_range[2]);
    }
  }
  cv::Mat resampled_source_mat(
      target_mat.size(), source_mat.type(), resampled_source_data.get());
  std::vector<cv::Mat> resampled_source_mats;
  cv::split(resampled_source_mat, resampled_source_mats);
  return resampled_source_mats;
}

void PansharpeningBase::MaskData(Data& data) {
  cv::Mat temp_mat, mask_mat;
  cv::threshold(data.mat, temp_mat, 0, 1, cv::THRESH_BINARY);
  cv::threshold(data.mats[0], mask_mat, 0, 1, cv::THRESH_BINARY);
  if (cv::countNonZero(temp_mat - mask_mat) != temp_mat.rows * temp_mat.cols) {
    spdlog::debug("Masking the upsampled data to the same area");
    cv::min(temp_mat, mask_mat, mask_mat);
    data.mat = data.mat.mul(mask_mat);
    for (auto& mat : data.mats)
      mat = mat.mul(mask_mat);
    spdlog::info("Masking the upsampled data to the same area - done");
  }
}

} // pansharpening
} // rs_toolset