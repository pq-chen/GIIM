#include "graph_cut_mosaicking.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>
#include <utility>

#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include "graph_cut_optimization/graph_cut_optimization.h"
#include <rs-toolset/mosaicking.h>
#include <rs-toolset/utils.hpp>

namespace rs_toolset {
namespace mosaicking {

GraphCutImpl::GraphCutImpl(
    float grad_self_low,
    float grad_self_high,
    float grad_self_exp,
    float diff_low,
    float diff_exp,
    double tol)
    : MosaickingBase(tol),
      grad_self_low_(grad_self_low),
      grad_self_high_(grad_self_high),
      grad_self_exp_(grad_self_exp),
      diff_low_(diff_low),
      diff_exp_(diff_exp) {
  spdlog::info(
      "Creating a graph cut mosaicking with\n"
      " - Gradient-self term low trunction: {}\n"
      " - Gradient-self term high trunction: {}\n"
      " - Gradient-self term exponential: {}\n"
      " - Difference term low trunction: {}\n"
      " - Difference term exponential: {}", grad_self_low, grad_self_high,
      grad_self_exp, diff_low, diff_exp);
}

void GraphCutImpl::PrepareData(
    GDALDataset* covered_overlap_dataset,
    GDALDataset* new_overlap_dataset,
    GDALDataset* label_raster_dataset,
    cv::Mat& covered_mat,
    cv::Mat& new_mat,
    cv::Mat& label_mat) {
  spdlog::debug("Preparing the data");
  covered_mat = CreateLightnessMat(covered_overlap_dataset);
  new_mat = CreateLightnessMat(new_overlap_dataset);
  if (label_raster_dataset) {
    label_mat = utils::CreateMatFromDataset(label_raster_dataset);
  } else {
    label_mat = cv::Mat::zeros(new_mat.size(), CV_8UC1);
  }
  spdlog::debug("Preparing the data - done");
}

void GraphCutImpl::ExecuteMosaicking(
    const cv::Mat& covered_mat,
    const cv::Mat& new_mat,
    const cv::Mat& label_mat,
    double* geotrans,
    OGRSpatialReference* spatial_ref,
    GDALDatasetUniquePtr& label_raster_dataset,
    OGRGeometryUniquePtr& label0_geometry,
    OGRGeometryUniquePtr& label1_geometry) {
  spdlog::debug("Execute the graph cut mosaicking algorithm");

  // Create maps between coors and valid indexes
  spdlog::debug("Creating maps between coors and valid indexes");
  std::map<std::pair<int, int>, GraphCutInfo> coor_to_idx;
  std::vector<std::pair<int, int>> idx_to_coor;
  int count(0);
  for (int row(0); row < covered_mat.rows; ++row) {
    for (int col(0); col < covered_mat.cols; ++col) {
      if (covered_mat.at<float>(row, col)) {
        if (new_mat.at<float>(row, col)) {
          coor_to_idx[{row, col}] = {count, 2};
        } else {
          coor_to_idx[{row, col}] = {count, 0};
        }
      } else {
        if (new_mat.at<float>(row, col)) {
          coor_to_idx[{row, col}] = {count, 1};
        } else {
          continue;
        }
      }
      idx_to_coor.push_back({row, col});
      ++count;
    }
  }
  spdlog::debug("Creating maps between coors and valid indexes - done");

  // Run the graph cut optimization algorithm
  spdlog::debug("Running the graph cut optimization algorithm");
  auto gco(std::make_unique<GraphCutOptimizationGeneralGraph>(count, 2));
  for (const auto& info : coor_to_idx) {
    switch (info.second.label) {
      case 0: {
        gco->SetDataCost(info.second.idx, 1, MAX_ENERGYTERM);
        gco->SetLabel(info.second.idx, 0);
        break;
      }
      case 1: {
        gco->SetDataCost(info.second.idx, 0, MAX_ENERGYTERM);
        gco->SetLabel(info.second.idx, 1);
        break;
      }
      case 2: {
        switch (label_mat.at<uint8_t>(info.first.first, info.first.second)) {
          case 100: {
            gco->SetLabel(info.second.idx, 0);
            break;
          }
          case 200: {
            gco->SetLabel(info.second.idx, 1);
          }
        }
      }
    }
    if (auto it(coor_to_idx.find({info.first.first, info.first.second - 1}));
        it != coor_to_idx.end()) {
      gco->SetNeighbors(info.second.idx, it->second.idx);
    }
    if (auto it(coor_to_idx.find({info.first.first - 1, info.first.second}));
        it != coor_to_idx.end()) {
      gco->SetNeighbors(info.second.idx, it->second.idx);
    }
  }
  cv::Mat
      x_kernel = (cv::Mat_<float>(4, 3) <<
          -0.4, 0, 0.4, -2, 0, 2, -2, 0, 2, -0.4, 0, 0.4),
      y_kernel = (cv::Mat_<float>(3, 4) <<
          -0.4, -2, -2, -0.4, 0, 0, 0, 0, 0.4, 2, 2, 0.4),
      covered_x_mat,
      covered_y_mat,
      new_x_mat,
      new_y_mat;
  cv::filter2D(covered_mat, covered_x_mat, CV_32F, x_kernel, cv::Point(1, 1));
  covered_x_mat = cv::abs(covered_x_mat);
  cv::filter2D(covered_mat, covered_y_mat, CV_32F, y_kernel, cv::Point(1, 1));
  covered_y_mat = cv::abs(covered_y_mat);
  cv::filter2D(new_mat, new_x_mat, CV_32F, x_kernel, cv::Point(1, 1));
  new_x_mat = cv::abs(new_x_mat);
  cv::filter2D(new_mat, new_y_mat, CV_32F, y_kernel, cv::Point(1, 1));
  new_y_mat = cv::abs(new_y_mat);
  Data data{
      grad_self_low_, grad_self_high_, grad_self_exp_, diff_low_, diff_exp_,
      covered_mat, covered_x_mat, covered_y_mat, new_mat, new_x_mat, new_y_mat,
      &idx_to_coor};
  gco->SetSmoothCost(SmoothCost, &data);
  gco->Expansion(2);
  spdlog::debug("Running the graph cut optimization algorithm - done");

  // Update the label raster
  spdlog::debug("Updating the label raster");
  cv::Mat new_label_mat(cv::Mat::zeros(label_mat.size(), CV_8UC1));
  for (const auto& info : coor_to_idx) {
    new_label_mat.at<uint8_t>(info.first.first, info.first.second) =
        gco->WhatLabel(info.second.idx) == 0 ? 100 : 200;
  }
  label_raster_dataset = utils::CreateDatasetFromMat(
      new_label_mat, "", geotrans, spatial_ref);
  spdlog::debug("Updating the label raster - done");

  // Polygonize the label raster
  spdlog::debug("Polygonizing the new label raster");
  GDALDatasetUniquePtr polygonized_vector_dataset(memory_driver_->Create(
      "", 0, 0, 0, GDT_Unknown, nullptr));
  auto polygonized_layer(polygonized_vector_dataset->CreateLayer(
      "", spatial_ref, wkbPolygon));
  OGRFieldDefn label_field("label", OFTInteger);
  polygonized_layer->CreateField(&label_field);
  GDALPolygonize(
      label_raster_dataset->GetRasterBand(1), nullptr, polygonized_layer, 0,
      nullptr, nullptr, nullptr);
  spdlog::debug("Polygonizing the new label raster - done");

  // Find label geometries among polygonized features
  spdlog::debug("Finding label geometries among polygonized features");
  std::pair label0_max_idx(-1, -1.0), label1_max_idx(-1, -1.0);
  for (int i(0); i < polygonized_layer->GetFeatureCount(); ++i) {
    OGRFeatureUniquePtr feature(polygonized_layer->GetFeature(i));
    switch (auto area(feature->GetGeometryRef()->toPolygon()->get_Area());
        feature->GetFieldAsInteger(0)) {
      case 100: {
        if (area > label0_max_idx.second)
          label0_max_idx = {i, area};
        break;
      }
      case 200: {
        if (area > label1_max_idx.second)
          label1_max_idx = {i, area};
      }
    }
  }
  OGRFeatureUniquePtr 
      feature0(polygonized_layer->GetFeature(label0_max_idx.first)),
      feature1(polygonized_layer->GetFeature(label1_max_idx.first));
  label0_geometry.reset(feature0->StealGeometry());
  auto _label0_geometry(label0_geometry->toPolygon());
  for (int i(0); i < _label0_geometry->getNumInteriorRings(); ++i)
    _label0_geometry->removeRing(i + 1);
  label1_geometry.reset(feature1->StealGeometry());
  label1_geometry.reset(label0_geometry->Union(label1_geometry.get()));
  auto _label1_geometry(label1_geometry->toPolygon());
  for (int i(0); i < _label1_geometry->getNumInteriorRings(); ++i)
    _label1_geometry->removeRing(i + 1);
  label1_geometry.reset(label1_geometry->Difference(label0_geometry.get()));
  spdlog::debug("Finding label geometries among polygonized features - done");
  spdlog::debug("Execute the graph cut mosaicking algorithm - done");
}

cv::Mat GraphCutImpl::CreateLightnessMat(GDALDataset* dataset) {
  int bands_map[3]{3, 2, 1};
  auto lab_mat(utils::CreateMatFromDataset(dataset, nullptr, 3, bands_map));
  lab_mat.convertTo(lab_mat, CV_32FC3, 1.0 / 255);
  cv::cvtColor(lab_mat, lab_mat, cv::COLOR_BGR2Lab);
  std::vector<cv::Mat> lab_mats;
  cv::split(lab_mat, lab_mats);
  return lab_mats[0];
}

int GraphCutImpl::SmoothCost(
    int site1,
    int site2,
    int label1,
    int label2,
    void* data) {
  if (label1 == label2) return 0;

  auto _data(static_cast<Data*>(data));
  auto row1(_data->idx_to_coor->at(site1).first),
      col1(_data->idx_to_coor->at(site1).second),
      row2(_data->idx_to_coor->at(site2).first),
      col2(_data->idx_to_coor->at(site2).second),
      min_row(std::min(row1, row2)),
      min_col(std::min(col1, col2));
  float data_value1(0.0f), data_value2(0.0f);
  if (label1 == 0) {
    data_value1 = _data->covered_mat.at<float>(row1, col1);
    data_value2 = _data->new_mat.at<float>(row2, col2);
  } else {
    data_value1 = _data->new_mat.at<float>(row1, col1);
    data_value2 = _data->covered_mat.at<float>(row2, col2);
  }
  float data_diff_term(fmax(fabs(data_value1 - data_value2), _data->diff_low)),
      grad_self_term(200.0f),
      grad_diff_term(200.0f);
  if (cv::Rect rect(
          min_col - 1, min_row - 1, 3 + (row1 == row2), 3 + (row1 != row2));
      min_row != 0 && std::max(row1, row2) != (_data->covered_mat.rows - 1) &&
      min_col != 0 && std::max(col1, col2) != (_data->covered_mat.cols - 1) &&
      cv::countNonZero(_data->covered_mat(rect)) == 12 &&
      cv::countNonZero(_data->new_mat(rect)) == 12) {
    float grad_value1, grad_value2;
    if (row1 == row2) {
      grad_value1 = _data->covered_y_mat.at<float>(min_row, min_col);
      grad_value2 = _data->new_y_mat.at<float>(min_row, min_col);
    } else {
      grad_value1 = _data->covered_x_mat.at<float>(min_row, min_col);
      grad_value2 = _data->new_x_mat.at<float>(min_row, min_col);
    }
    grad_self_term = pow(
        fmin(
            fmax(fmax(grad_value1, grad_value2), _data->grad_self_low),
            _data->grad_self_high), _data->grad_self_exp);
    grad_diff_term = fmax(
        fabs(grad_value1 - grad_value2) / 2.8f, _data->diff_low);
  }
  return static_cast<int>(round(
      100 * grad_self_term *
          pow(data_diff_term + grad_diff_term, _data->diff_exp)));
}

std::shared_ptr<GraphCut> GraphCut::Create(
    float grad_self_low,
    float grad_self_high,
    float grad_self_exp,
    float diff_low,
    float diff_exp,
    double tol) {
  return std::make_shared<GraphCutImpl>(
      grad_self_low, grad_self_high, grad_self_exp, diff_low, diff_exp, tol);
}

}  // namespace mosaicking
}  // namespace rs_toolset