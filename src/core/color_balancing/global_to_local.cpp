#include "global_to_local.h"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SPQRSupport>
#include <gdal_priv.h>
#include <nlohmann/json.hpp>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/utils.hpp>


namespace fs = std::filesystem;

namespace rs_toolset {
namespace color_balancing {

bool GlobalToLocalImpl::CalcArgus(
    const std::vector<std::string>& rasters_path) {
  for (const auto& path : rasters_path)
    if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
            path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)); !dataset) {
      spdlog::warn("Opening {} failed", path);
      return false;
    }
  std::string string("Calculating arguments for the following raster(s):\n");
  for (const auto& path : rasters_path)
    string.append(path).append("\n");
  string.append("{} task(s) in total");
  spdlog::info(string, rasters_path.size());
  rasters_path_ = rasters_path;

  // Calculate the rasters' statistic
  spdlog::info("Calculating the rasters' statistic");
  int rasters_count(static_cast<int>(rasters_path.size()));
  std::vector<OGRGeometryUniquePtr> borders(rasters_count);
  std::vector<Statistic> rasters_statistic(rasters_count);
  double 
      extent[4]{ DBL_MAX, -DBL_MAX, -DBL_MAX, DBL_MAX },
      resolution(DBL_MAX);
  for (int i = 0; i < rasters_count; i++) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        rasters_path[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    borders[i] = utils::CreateBorder(dataset.get());
    int x_size(dataset->GetRasterXSize()),
        y_size(dataset->GetRasterYSize()),
        block_cols_count(static_cast<int>(
            ceil(static_cast<double>(x_size) / block_size_))),
        block_rows_count(static_cast<int>(
            ceil(static_cast<double>(y_size) / block_size_))),
        last_block_x_size(x_size - (block_cols_count - 1) * block_size_),
        last_block_y_size(y_size - (block_rows_count - 1) * block_size_);
    for (int k = 0; k < block_cols_count * block_rows_count; k++) {
      int range[4];
      utils::CreateRange(
          k, block_cols_count, block_rows_count, block_size_, block_size_,
          last_block_x_size, last_block_y_size, range);
      UpdateStatistic(
          CreateRasterData(dataset.get(), range), rasters_statistic[i]);
    }
    CalcStatistic(rasters_statistic[i]);
    spdlog::info("---------- {}/{} - done ----------", i + 1, rasters_count);
    double geotrans[6];
    dataset->GetGeoTransform(geotrans);
    double
        max_x(geotrans[0] + dataset->GetRasterXSize() * geotrans[1]),
        min_y(geotrans[3] + dataset->GetRasterYSize() * geotrans[5]);
    if (extent[0] > geotrans[0]) extent[0] = geotrans[0];
    if (extent[1] < geotrans[3]) extent[1] = geotrans[3];
    if (extent[2] < max_x) extent[2] = max_x;
    if (extent[3] > min_y) extent[3] = min_y;
    if (resolution > geotrans[1]) resolution = geotrans[1];
  }
  spdlog::info("Calculating the rasters' statistic - done");

  // Calculate the overlaps' statistic
  spdlog::info("Calculating the overlaps' statistic");
  std::vector<OverlapInfo> overlaps_info;
  std::vector<std::vector<int>> rasters_overlap_idxes(rasters_count);
  for (int i = 0; i < rasters_count; i++)
    for (int j = i + 1; j < rasters_count; j++)
      if (borders[i]->Intersect(borders[j].get())) {
        rasters_overlap_idxes[i].push_back(static_cast<int>(
            overlaps_info.size()));
        rasters_overlap_idxes[j].push_back(static_cast<int>(
            overlaps_info.size()));
        overlaps_info.push_back({ i, j, Statistic(), Statistic(), 0.0 });
      }
  for (int i = 0; i < rasters_count; i++) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        rasters_path[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    double geotrans[6], inv_geotrans[6];
    dataset->GetGeoTransform(geotrans);
    GDALInvGeoTransform(geotrans, inv_geotrans);
    for (const auto& j : rasters_overlap_idxes[i]) {
      OverlapInfo& info(overlaps_info[j]);
      auto& s(i == info.idx1 ? info.s1 : info.s2);
      OGRGeometryUniquePtr geoemtry(borders[i]->Intersection(
          borders[i == info.idx1 ? info.idx2 : info.idx1].get()));
      if (!info.weight)
        info.weight = geoemtry->toPolygon()->get_Area();
      OGREnvelope enve;
      geoemtry->getEnvelope(&enve);
      double x, y;
      int out_range[4];
      GDALApplyGeoTransform(inv_geotrans, enve.MinX, enve.MaxY, &x, &y);
      out_range[0] = static_cast<int>(ceil(x));
      out_range[1] = static_cast<int>(ceil(y));
      GDALApplyGeoTransform(inv_geotrans, enve.MaxX, enve.MinY, &x, &y);
      out_range[2] = static_cast<int>(floor(x)) - out_range[0];
      out_range[3] = static_cast<int>(floor(y)) - out_range[1];
      int block_cols_count(static_cast<int>(
              ceil(static_cast<double>(out_range[2]) / block_size_))),
          block_rows_count(static_cast<int>(
              ceil(static_cast<double>(out_range[3]) / block_size_))),
          last_block_x_size(out_range[2] - (block_cols_count - 1) *
              block_size_),
          last_block_y_size(out_range[3] - (block_rows_count - 1) *
              block_size_);
      for (int k = 0; k < block_cols_count * block_rows_count; k++) {
        int range[4];
        utils::CreateRange(
            k, block_cols_count, block_rows_count, block_size_, block_size_,
            last_block_x_size, last_block_y_size, range);
        range[0] += out_range[0];
        range[1] += out_range[1];
        UpdateStatistic(
            CreateOverlapData(dataset.get(), range, geoemtry.get()), s);
      }
      CalcStatistic(s);
    }
    spdlog::info("---------- {}/{} - done ----------", i + 1, rasters_count);
  }
  spdlog::info("Calculating the overlaps' statistic - done");

  // Calculate rasters' global arguments
  CalcRastersGlobalArgus(rasters_statistic, overlaps_info);

  // Calculate the local grid information
  spdlog::info("Calculating the local grid information");
  double grid_geotrans[6]{
      extent[0], resolution * grid_size_, 0, extent[1], 0,
      -resolution * grid_size_ };
  GDALInvGeoTransform(grid_geotrans, inv_grid_geotrans_);
  int grid_cols_count(static_cast<int>(
          ceil((extent[2] - extent[0]) / grid_geotrans[1]))),
      grid_rows_count(static_cast<int>(
          ceil((extent[3] - extent[1]) / grid_geotrans[5])));
  rasters_grid_mean_mat_.clear();
  for (int i = 0; i < rasters_count; i++)
    rasters_grid_mean_mat_.push_back(
        cv::Mat::zeros(grid_rows_count, grid_cols_count, CV_32FC3));
  cv::Mat global_grid_pixels_count_mat(
      cv::Mat::zeros(grid_rows_count, grid_cols_count, CV_32FC1));
  global_grid_mean_mat_ =
      cv::Mat::zeros(grid_rows_count, grid_cols_count, CV_32FC3);
  max_mean_ = 0.0;
  for (int i = 0; i < rasters_count; i++) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        rasters_path[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    int x_size(dataset->GetRasterXSize()),
        y_size(dataset->GetRasterYSize()),
        block_cols_count(static_cast<int>(
            ceil(static_cast<double>(x_size) / block_size_))),
        block_rows_count(static_cast<int>(
            ceil(static_cast<double>(y_size) / block_size_))),
        last_block_x_size(x_size - (block_cols_count - 1) * block_size_),
        last_block_y_size(y_size - (block_rows_count - 1) * block_size_);
    double
        geotrans[6], inv_geotrans[6], composed_geotrans[6],
        inv_composed_geotrans[6];
    dataset->GetGeoTransform(geotrans);
    GDALInvGeoTransform(geotrans, inv_geotrans);
    GDALComposeGeoTransforms(geotrans, inv_grid_geotrans_, composed_geotrans);
    GDALComposeGeoTransforms(
        grid_geotrans, inv_geotrans, inv_composed_geotrans);
    for (int k = 0; k < block_cols_count * block_rows_count; k++) {
      int raster_range[4];
      utils::CreateRange(
          k, block_cols_count, block_rows_count, block_size_, block_size_,
          last_block_x_size, last_block_y_size, raster_range);
      auto mats(CreateRasterData(dataset.get(), raster_range));
      double extent[4]{
          static_cast<double>(raster_range[0]),
          static_cast<double>(raster_range[1]),
          static_cast<double>(raster_range[2] + raster_range[0]),
          static_cast<double>(raster_range[3] + raster_range[1]) };
      GDALApplyGeoTransform(
          composed_geotrans, extent[0], extent[1], &extent[0], &extent[1]);
      GDALApplyGeoTransform(
          composed_geotrans, extent[2], extent[3], &extent[2], &extent[3]);
      int grid_range[4]{
          static_cast<int>(floor(extent[0])),
          static_cast<int>(floor(extent[1])),
          static_cast<int>(ceil(extent[2])),
          static_cast<int>(ceil(extent[3])) };
      grid_range[2] -= grid_range[0];
      grid_range[3] -= grid_range[1];
      UpdateGridInfo(
          i, mats, raster_range, grid_range, inv_composed_geotrans,
          global_grid_pixels_count_mat);
    }
    spdlog::info("---------- {}/{} - done ----------", i + 1, rasters_count);
  }
  spdlog::info("Calculating the local grid information - done");
  spdlog::info("Calculating arguments - done");
  return true;
}

bool GlobalToLocalImpl::ExportArgusJson(const std::string& json_path) {
  fs::path path(json_path);
  if (!fs::exists(path.parent_path())) {
    spdlog::warn(
        "The JSON directory {} does not exist", path.parent_path().string());
    return false;
  }
  if (path.extension().string() != ".json") {
    spdlog::warn("The JSON path {} does not end with \".json\"", json_path);
    return false;
  }
  if (rasters_path_.empty()) {
    spdlog::warn("No arguments can be exported");
    return false;
  }

  spdlog::info("Exporting arguments to the JSON file {}", json_path);
  nlohmann::json json;
  std::vector<std::vector<uint8_t>> rasters_grid_means(rasters_path_.size());
  for (int i = 0; i < rasters_path_.size(); i++)
    cv::imencode(".tif", rasters_grid_mean_mat_[i], rasters_grid_means[i]);
  std::vector<uint8_t> global_grid_mean;
  cv::imencode(".tif", global_grid_mean_mat_, global_grid_mean);
  std::vector<double> inv_grid_geotrans;
  inv_grid_geotrans.assign(inv_grid_geotrans_, inv_grid_geotrans_ + 6);
  json["method"] = "global-to-local";
  json["rasters_paths"] = rasters_path_;
  json["rasters_global_argus"] = rasters_global_argus_;
  json["rasters_grid_means"] = rasters_grid_means;
  json["global_grid_means"] = global_grid_mean;
  json["max_mean"] = max_mean_;
  json["inv_grid_geotrans"] = inv_grid_geotrans;
  std::ofstream o(json_path);
  o << json << std::endl;
  spdlog::info("Exporting arguments - done");
  return true;
}

bool GlobalToLocalImpl::ImportArgusJson(const std::string& json_path) {
  std::ifstream i(json_path);
  if (!i) {
    spdlog::warn("Opening {} failed", json_path);
    return false;
  }

  spdlog::info("Importing arguments from the JSON file {}", json_path);
  nlohmann::json json;
  i >> json;
  try {
    if (json.at("method").get<std::string>() != "global-to-local")
      return false;
    json.at("rasters_paths").get_to(rasters_path_);
    for (const auto& path : rasters_path_)
      if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
              path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY)); !dataset) {
        spdlog::warn("Opening {} failed", path);
        return false;
      }
    json.at("rasters_global_argus").get_to(rasters_global_argus_);
    rasters_grid_mean_mat_.clear();
    for (const auto& mean : json.at("rasters_grid_means"))
      rasters_grid_mean_mat_.push_back(cv::imdecode(
          mean.get<std::vector<uint8_t>>(), cv::IMREAD_UNCHANGED));
    global_grid_mean_mat_ = cv::imdecode(
        json.at("global_grid_means").get<std::vector<uint8_t>>(),
        cv::IMREAD_UNCHANGED);
    json["max_mean"].get_to(max_mean_);
    auto inv_grid_geotrans(
        json.at("inv_grid_geotrans").get<std::vector<double>>());
    memcpy(inv_grid_geotrans_, inv_grid_geotrans.data(), 6 * sizeof(double));
  } catch (...) {
    spdlog::warn("The json file {} has been broken", json_path);
    return false;
  }
  spdlog::info("Importing arguments - done");
  return true;
}

std::vector<cv::Mat> GlobalToLocalImpl::CreateOverlapData(
    GDALDataset* dataset,
    int* range,
    OGRGeometry* geometry) {
  double geotrans[6];
  dataset->GetGeoTransform(geotrans);
  geotrans[0] += range[0] * geotrans[1];
  geotrans[3] += range[1] * geotrans[5];
  GDALDriver* mem_driver(GetGDALDriverManager()->GetDriverByName("MEM"));
  GDALDatasetUniquePtr output_dataset(mem_driver->Create(
      "", range[2], range[3], 3, GDT_Byte, nullptr));
  output_dataset->SetGeoTransform(geotrans);
  output_dataset->SetSpatialRef(dataset->GetSpatialRef());
  utils::WarpByGeometry({ dataset }, { geometry }, output_dataset);
  std::vector<cv::Mat> mats;
  cv::split(utils::CreateMatFromDataset(output_dataset.get()), mats);
  return mats;
}

void GlobalToLocalImpl::UpdateStatistic(
    const std::vector<cv::Mat>& mats,
    Statistic& s) {
  int pixels_count(cv::countNonZero(mats[0]));
  if (!pixels_count) return;
  s.pixels_count_ += pixels_count;
  for (int b = 0; b < 3; b++) {
    s.sums_[b] += cv::sum(mats[b])[0];
    cv::Mat temp;
    mats[b].convertTo(temp, CV_32S);
    s.square_sums_[b] += cv::sum(temp.mul(temp))[0];
  }
}

void GlobalToLocalImpl::CalcRastersGlobalArgus(
    const std::vector<Statistic>& rasters_statistic,
    const std::vector<OverlapInfo>& overlaps_info) {
  spdlog::info("Calculating rasters' global arguments");
  int rasters_count(static_cast<int>(rasters_statistic.size())),
      overlaps_count(static_cast<int>(overlaps_info.size()));
  rasters_global_argus_ =
      std::vector<std::vector<float>>(rasters_count, std::vector<float>(6));
#pragma omp parallel for schedule(static, 3)
  for (int _b = 0; _b < 3; _b++) {
    std::vector<Eigen::Triplet<double>> A_coefs, sqrt_P_coefs;
    Eigen::VectorXd b(
        Eigen::VectorXd::Zero(2 * (overlaps_count + rasters_count)));
    double sqrt_raster_weight(0.0);
    for (int i = 0; i < overlaps_count; i++) {
      auto& info(overlaps_info[i]);
      A_coefs.emplace_back(2 * i, 2 * info.idx1, info.s1.means[_b]);
      A_coefs.emplace_back(2 * i, 2 * info.idx1 + 1, 1.0);
      A_coefs.emplace_back(2 * i, 2 * info.idx2, -info.s2.means[_b]);
      A_coefs.emplace_back(2 * i, 2 * info.idx2 + 1, -1.0);
      A_coefs.emplace_back(2 * i + 1, 2 * info.idx1, info.s1.stddevs[_b]);
      A_coefs.emplace_back(2 * i + 1, 2 * info.idx2, -info.s2.stddevs[_b]);
      sqrt_P_coefs.emplace_back(2 * i, 2 * i, sqrt(info.weight));
      sqrt_P_coefs.emplace_back(2 * i + 1, 2 * i + 1, sqrt(info.weight));
      sqrt_raster_weight += info.weight;
    }
    sqrt_raster_weight =
        sqrt_raster_weight ? sqrt(sqrt_raster_weight / (overlaps_count * 2))
        : 1.0;
    for (int i = overlaps_count; i < overlaps_count + rasters_count; i++) {
      int idx(i - overlaps_count);
      auto& s(rasters_statistic[idx]);
      A_coefs.emplace_back(2 * i, 2 * idx, s.means[_b]);
      A_coefs.emplace_back(2 * i, 2 * idx + 1, 1.0);
      A_coefs.emplace_back(2 * i + 1, 2 * idx, s.stddevs[_b]);
      sqrt_P_coefs.emplace_back(2 * i, 2 * i, sqrt_raster_weight);
      sqrt_P_coefs.emplace_back(2 * i + 1, 2 * i + 1, sqrt_raster_weight);
      b(2 * i) = s.means[_b];
      b(2 * i + 1) = s.stddevs[_b];
    }
    Eigen::SparseMatrix<double>
        A(2 * (overlaps_count + rasters_count), 2 * rasters_count),
        sqrt_P(
            2 * (overlaps_count + rasters_count),
            2 * (overlaps_count + rasters_count));
    A.setFromTriplets(A_coefs.begin(), A_coefs.end());
    sqrt_P.setFromTriplets(sqrt_P_coefs.begin(), sqrt_P_coefs.end());
    A = sqrt_P * A;
    b = sqrt_P * b;
    Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd solve(solver.solve(b));
    for (int i = 0; i < rasters_count; i++) {
      rasters_global_argus_[i][2 * _b] = static_cast<float>(solve(2 * i));
      rasters_global_argus_[i][2 * _b + 1] =
          static_cast<float>(solve(2 * i + 1));
    }
  }
  spdlog::info("Calculating rasters' global arguments - done");
}

void GlobalToLocalImpl::UpdateGridInfo(
    int idx,
    const std::vector<cv::Mat>& mats,
    int* raster_range,
    int* grid_range,
    double* composed_geotrans,
    cv::Mat& global_grid_pixels_count_mat) {
#pragma omp parallel for schedule(dynamic)
  for (int row = 0; row < grid_range[3]; row++) {
    auto raster_grid_mean(rasters_grid_mean_mat_[idx].ptr<float>(
            grid_range[1] + row, grid_range[0])),
        global_grid_pixels_count(global_grid_pixels_count_mat.ptr<float>(
            grid_range[1] + row, grid_range[0])),
        global_grid_mean(global_grid_mean_mat_.ptr<float>(
            grid_range[1] + row, grid_range[0]));
    for (int col = 0; col < grid_range[2]; col++, global_grid_pixels_count++) {
      double extent[4] {
          static_cast<double>(grid_range[0] + col),
          static_cast<double>(grid_range[1] + row),
          static_cast<double>(grid_range[0] + col + 1),
          static_cast<double>(grid_range[1] + row + 1) };
      GDALApplyGeoTransform(
          composed_geotrans, extent[0], extent[1], &extent[0], &extent[1]);
      extent[0] -= raster_range[0];
      extent[1] -= raster_range[1];
      GDALApplyGeoTransform(
          composed_geotrans, extent[2], extent[3], &extent[2], &extent[3]);
      extent[2] -= raster_range[0];
      extent[3] -= raster_range[1];
      cv::Rect rect(
          static_cast<int>(fmax(fmin(round(extent[0]), mats[0].cols - 1), 0)),
          static_cast<int>(fmax(fmin(round(extent[1]), mats[0].rows - 1), 0)),
          static_cast<int>(fmax(fmin(round(extent[2]), mats[0].cols - 1), 0)),
          static_cast<int>(fmax(fmin(round(extent[3]), mats[0].rows - 1), 0)));
      rect.width -= rect.x;
      rect.height -= rect.y;
      int cur_pixels_count(cv::countNonZero(mats[0](rect)));
      if (cur_pixels_count) {
        for (int b = 0; b < 3; b++, raster_grid_mean++, global_grid_mean++) {
          *raster_grid_mean = fmax(fmin(
              rasters_global_argus_[idx][2 * b + 1] +
              rasters_global_argus_[idx][2 * b] *
              (cv::sum(mats[b](rect))[0] / cur_pixels_count), 255.0), 1.0);
          *global_grid_mean =
              (*global_grid_mean * *global_grid_pixels_count +
              *raster_grid_mean * cur_pixels_count) /
              (*global_grid_pixels_count + cur_pixels_count);
          if (max_mean_ < *raster_grid_mean)
            max_mean_ = *raster_grid_mean;
        }
        *global_grid_pixels_count += cur_pixels_count;
      } else {
        raster_grid_mean += 3;
        global_grid_mean += 3;
      }
    }
  }
}

void GlobalToLocalImpl::WriteToDataset(
    int idx,
    GDALDataset* source_dataset,
    GDALDatasetUniquePtr& output_dataset) {
  double geotrans[6], composed_geotrans[6];
  source_dataset->GetGeoTransform(geotrans);
  GDALComposeGeoTransforms(geotrans, inv_grid_geotrans_, composed_geotrans);
  int x_size(source_dataset->GetRasterXSize()),
      y_size(source_dataset->GetRasterYSize()),
      block_cols_count(static_cast<int>(
          ceil(static_cast<double>(x_size) / block_size_))),
      block_rows_count(static_cast<int>(
          ceil(static_cast<double>(y_size) / block_size_))),
      last_block_x_size(x_size - (block_cols_count - 1) * block_size_),
      last_block_y_size(y_size - (block_rows_count - 1) * block_size_);
  for (int k = 0; k < block_cols_count * block_rows_count; k++) {
    int range[4];
    utils::CreateRange(
        k, block_cols_count, block_rows_count, block_size_, block_size_,
        last_block_x_size, last_block_y_size, range);
    auto mats(CreateRasterData(source_dataset, range));
    ExecuteColorBalancing(idx, range, composed_geotrans, mats);
#pragma omp parallel for schedule(static, 3)
    for (int b = 0; b < 3; b++)
      output_dataset->GetRasterBand(b + 1)->RasterIO(
          GF_Write, range[0], range[1], range[2], range[3], mats[b].data,
          range[2], range[3], GDT_Byte, 1, range[2]);
  }
}

void GlobalToLocalImpl::ExecuteColorBalancing(
    int idx,
    int* range,
    double* composed_geotrans,
    std::vector<cv::Mat>& mats) {
//#pragma omp parallel for schedule(dynamic)
  for (int row = 0; row < range[3]; row++) {
    auto resampled_raster_grid_data(std::make_unique<float[]>(3 * range[2])),
        resampled_global_grid_data(std::make_unique<float[]>(3 * range[2]));
    auto xs(std::make_unique<double[]>(range[2])),
        ys(std::make_unique<double[]>(range[2])),
        zs(std::make_unique<double[]>(range[2]));
    auto successes(std::make_unique<int[]>(range[2]));
    memset(zs.get(), 0, range[2] * sizeof(double));
    for (int col = 0; col < range[2]; col++)
      GDALApplyGeoTransform(
          composed_geotrans, col + range[0], row + range[1], &xs[col],
          &ys[col]);
    utils::CalcInterpolations<float>(
        rasters_grid_mean_mat_[idx], range[2], xs.get(), ys.get(),
        resampled_raster_grid_data.get());
    utils::CalcInterpolations<float>(
        global_grid_mean_mat_, range[2], xs.get(), ys.get(),
        resampled_global_grid_data.get());
    cv::Mat
        resampled_grid_mat(
            1, range[2], CV_32FC3, resampled_raster_grid_data.get()),
        resampled_global_grid_mat(
            1, range[2], CV_32FC3, resampled_global_grid_data.get());
    cv::log(resampled_grid_mat / max_mean_, resampled_grid_mat);
    cv::log(resampled_global_grid_mat / max_mean_, resampled_global_grid_mat);
    cv::Mat exponential_mat(resampled_global_grid_mat / resampled_grid_mat);
    auto ptr(exponential_mat.ptr<float>(0));
    for (int col = 0; col < range[2]; col++)
      for (int b = 0; b < 3; b++) {
        float value(static_cast<float>(mats[b].at<uint8_t>(row, col)));
        if (!value) continue;
        value = rasters_global_argus_[idx][2 * b] * value +
            rasters_global_argus_[idx][2 * b + 1];
        mats[b].at<uint8_t>(row, col) =
            value < 0 ? 1 : static_cast<uint8_t>(round(fmax(fmin(pow(
                value / max_mean_, *(ptr + 3 * col + b)) * max_mean_,
                255.0), 1.0)));
      }
  }
}

std::shared_ptr<GlobalToLocal> GlobalToLocal::Create(int block_size) {
  return std::make_shared<GlobalToLocalImpl>(block_size);
}

} // namespace color_balancing
} // namespace rs_toolset