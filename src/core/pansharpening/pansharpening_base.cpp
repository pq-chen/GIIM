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
    int source_x_size,
    int source_y_size,
    int target_block_idx,
    int target_block_cols_count,
    int target_block_rows_count,
    int target_block_x_size,
    int target_block_y_size,
    int last_target_block_x_size,
    int last_target_block_y_size,
    const utils::RPCTransPtr& source_trans_arg,
    const utils::RPCTransPtr& target_trans_arg,
    int(&source_range)[4],
    int(&target_range)[4]) {
  utils::CreateRange(
      target_block_idx, target_block_cols_count, target_block_rows_count,
      target_block_x_size, target_block_y_size, last_target_block_x_size,
      last_target_block_y_size, target_range);
  double
      xs[4]{
          static_cast<double>(target_range[0]),
          static_cast<double>(target_range[0] + target_range[2]),
          static_cast<double>(target_range[0]),
          static_cast<double>(target_range[0] + target_range[2])},
      ys[4]{
          static_cast<double>(target_range[1]),
          static_cast<double>(target_range[1]),
          static_cast<double>(target_range[1] + target_range[3]),
          static_cast<double>(target_range[1] + target_range[3])},
      zs[4]{0.0, 0.0, 0.0, 0.0};
  int successes[4];
  GDALRPCTransform(target_trans_arg.get(), false, 4, xs, ys, zs, successes);
  GDALRPCTransform(source_trans_arg.get(), true, 4, xs, ys, zs, successes);
  source_range[0] = static_cast<int>(floor(fmax(fmin(xs[0], xs[2]), 0.0)));
  source_range[1] = static_cast<int>(floor(fmax(fmin(ys[0], ys[1]), 0.0)));
  source_range[2] = static_cast<int>(ceil(fmin(fmax(xs[1], xs[3]),
      static_cast<double>(source_x_size)))) - source_range[0];
  source_range[3] = static_cast<int>(ceil(fmin(fmax(ys[2], ys[3]),
      static_cast<double>(source_y_size)))) - source_range[1];
}

PansharpeningBase::Data PansharpeningBase::CreateDownsampledDataWithRPC(
    GDALDataset* pan_dataset,
    GDALDataset* ms_dataset,
    int* pan_range,
    int* ms_range,
    const utils::RPCTransPtr& pan_trans_arg,
    const utils::RPCTransPtr& ms_trans_arg) {
  Data data;
  data.pan_mat = CreateResampledMatsWithRpc(
      utils::CreateMatFromDataset(pan_dataset, pan_range), pan_range, ms_range,
      pan_trans_arg, ms_trans_arg)[0];
  cv::split(utils::CreateMatFromDataset(ms_dataset, ms_range), data.ms_mats);
  MaskData(data);
  return data;
}

PansharpeningBase::Data PansharpeningBase::CreateDownsampledDataWithGeotrans(
    GDALDataset* pan_dataset,
    GDALDataset* ms_dataset,
    int* ms_range) {
  auto driver(GetGDALDriverManager()->GetDriverByName("MEM"));
  GDALDatasetUniquePtr downsampled_pan_dataset(driver->Create(
      "", ms_range[2], ms_range[3], pan_dataset->GetRasterCount(),
      pan_dataset->GetRasterBand(1)->GetRasterDataType(), nullptr));
  double geotrans[6];
  ms_dataset->GetGeoTransform(geotrans);
  geotrans[0] += ms_range[0] * geotrans[1];
  geotrans[3] += ms_range[1] * geotrans[5];
  downsampled_pan_dataset->SetGeoTransform(geotrans);
  downsampled_pan_dataset->SetSpatialRef(pan_dataset->GetSpatialRef());
  utils::WarpByGeometry({pan_dataset}, {nullptr}, downsampled_pan_dataset);
  Data data;
  data.pan_mat = utils::CreateMatFromDataset(downsampled_pan_dataset.get());
  cv::split(utils::CreateMatFromDataset(ms_dataset, ms_range), data.ms_mats);
  MaskData(data);
  return data;
}

PansharpeningBase::Data PansharpeningBase::CreateUpsampledDataWithRPC(
    GDALDataset* pan_dataset,
    GDALDataset* ms_dataset,
    int* pan_range,
    int* ms_range,
    const utils::RPCTransPtr& pan_trans_arg,
    const utils::RPCTransPtr& ms_trans_arg) {
  Data data;
  data.pan_mat = utils::CreateMatFromDataset(pan_dataset, pan_range);
  data.ms_mats = CreateResampledMatsWithRpc(
      utils::CreateMatFromDataset(ms_dataset, ms_range), ms_range, pan_range,
      ms_trans_arg, pan_trans_arg);
  MaskData(data);
  return data;
}

PansharpeningBase::Data PansharpeningBase::CreateUpsampledDataWithGeotrans(
    GDALDataset* pan_dataset,
    GDALDataset* ms_dataset,
    int* pan_range) {
  auto mem_driver(GetGDALDriverManager()->GetDriverByName("MEM"));
  GDALDatasetUniquePtr upsampled_ms_dataset(mem_driver->Create(
      "", pan_range[2], pan_range[3], ms_dataset->GetRasterCount(),
      ms_dataset->GetRasterBand(1)->GetRasterDataType(), nullptr));
  double geotrans[6];
  pan_dataset->GetGeoTransform(geotrans);
  geotrans[0] += pan_range[0] * geotrans[1];
  geotrans[3] += pan_range[1] * geotrans[5];
  upsampled_ms_dataset->SetGeoTransform(geotrans);
  upsampled_ms_dataset->SetSpatialRef(ms_dataset->GetSpatialRef());
  utils::WarpByGeometry({ms_dataset}, {nullptr}, upsampled_ms_dataset);
  Data data;
  data.pan_mat = utils::CreateMatFromDataset(pan_dataset, pan_range);
  cv::split(
      utils::CreateMatFromDataset(upsampled_ms_dataset.get()), data.ms_mats);
  MaskData(data);
  return data;
}

std::vector<cv::Mat> PansharpeningBase::CreateResampledMatsWithRpc(
    const cv::Mat& source_mat,
    int* source_range,
    int* target_range,
    const utils::RPCTransPtr& source_trans_arg,
    const utils::RPCTransPtr& target_trans_arg) {
  auto bands_count(source_mat.channels());
  auto size(static_cast<uint64_t>(target_range[2]) * target_range[3]);
  std::unique_ptr<void, void(*)(void*)> data(
      nullptr, [](void* p) { delete[] p; });
  if (source_mat.depth() == CV_8U) {
    data.reset(static_cast<void*>(new uint8_t[bands_count * size]));
  } else {
    data.reset(static_cast<void*>(new uint16_t[bands_count * size]));
  }
#pragma omp parallel for schedule(dynamic)
  for (int row(0); row < target_range[3]; ++row) {
    auto xs(std::make_unique<double[]>(target_range[2])),
        ys(std::make_unique<double[]>(target_range[2])),
        zs(std::make_unique<double[]>(target_range[2]));
    auto successes(std::make_unique<int[]>(target_range[2]));
    memset(zs.get(), 0, sizeof(double) * target_range[2]);
    for (int col(0); col < target_range[2]; ++col) {
      xs[col] = col + target_range[0];
      ys[col] = row + target_range[1];
    }
    GDALRPCTransform(
        target_trans_arg.get(), false, target_range[2], xs.get(), ys.get(),
        zs.get(), successes.get());
    GDALRPCTransform(
        source_trans_arg.get(), true, target_range[2], xs.get(), ys.get(),
        zs.get(), successes.get());
    for (int col(0); col < target_range[2]; ++col) {
      xs[col] -= source_range[0];
      ys[col] -= source_range[1];
    }

    // Calculate all rows and columns' interpolations in the source mat
    if (source_mat.depth() == CV_8U) {
      utils::CalcInterpolations(
          source_mat, target_range[2], xs.get(), ys.get(),
          static_cast<uint8_t*>(data.get()) +
          static_cast<uint64_t>(bands_count) * row * target_range[2], true);
    } else {
      utils::CalcInterpolations(
          source_mat, target_range[2], xs.get(), ys.get(),
          static_cast<uint16_t*>(data.get()) +
          static_cast<uint64_t>(bands_count) * row * target_range[2], true);
    }
  }
  std::vector<cv::Mat> output_mats;
  cv::split(
      cv::Mat(target_range[3], target_range[2], source_mat.type(), data.get()),
      output_mats);
  return output_mats;
}

void PansharpeningBase::MaskData(Data& data) {
  cv::Mat pan_mask_mat, ms_mask_mat;
  cv::threshold(data.pan_mat, pan_mask_mat, 0, 1, cv::THRESH_BINARY);
  cv::threshold(data.ms_mats[0], ms_mask_mat, 0, 1, cv::THRESH_BINARY);
  if (cv::countNonZero(cv::abs(pan_mask_mat - ms_mask_mat))) {
    cv::min(pan_mask_mat, ms_mask_mat, pan_mask_mat);
    data.pan_mat = data.pan_mat.mul(pan_mask_mat);
    for (auto& mat : data.ms_mats)
      mat = mat.mul(pan_mask_mat);
  }
}

}  // namespace pansharpening
}  // namespace rs_toolset