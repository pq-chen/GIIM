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

#include "graph_cut_optimization/GCoptimization.h"
#include <rs-toolset/mosaicking.h>
#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace mosaicking {

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

  // Create graph cut data structures and run the graph cut algorithm
  spdlog::debug("Creating maps between coors and valid indexes");
  std::map<std::pair<int, int>, GraphCutInfo> coor_to_valid_idxes;
  std::map<int, std::pair<int, int>> valid_idx_to_coors;
  int valid_idx(0);
  for (int row = 0; row < covered_mat.rows; row++)
    for (int col = 0; col < covered_mat.cols; col++) {
      if (covered_mat.at<float>(row, col)) {
        if (new_mat.at<float>(row, col)) {
          coor_to_valid_idxes[{row, col}] = { valid_idx, 2 };
        } else {
          coor_to_valid_idxes[{row, col}] = { valid_idx, 0 };
        }
      } else {
        if (new_mat.at<float>(row, col)) {
          coor_to_valid_idxes[{row, col}] = { valid_idx, 1 };
        } else {
          continue;
        }
      }
      valid_idx_to_coors[valid_idx] = { row, col };
      valid_idx++;
    }
  spdlog::debug("Creating maps between coors and valid indexes - done");
  spdlog::debug("Creating a graph cut class and run the algorithm");
  auto gc(std::make_unique<GCoptimizationGeneralGraph>(valid_idx, 2));
  decltype(coor_to_valid_idxes.begin()) it;
  for (const auto& info : coor_to_valid_idxes) {
    switch (info.second.label) {
      case 0: {
        gc->setDataCost(info.second.idx, 1, GCO_MAX_ENERGYTERM);
        gc->setLabel(info.second.idx, 0);
        break;
      }
      case 1: {
        gc->setDataCost(info.second.idx, 0, GCO_MAX_ENERGYTERM);
        gc->setLabel(info.second.idx, 1);
        break;
      }
      case 2: {
        switch (label_mat.at<uint8_t>(info.first.first, info.first.second)){
          case 100: {
            gc->setLabel(info.second.idx, 0);
            break;
          }
          case 200: {
            gc->setLabel(info.second.idx, 1);
          }
        }
      }
    }
    it = coor_to_valid_idxes.find({ info.first.first, info.first.second - 1 });
    if (it != coor_to_valid_idxes.end())
      gc->setNeighbors(info.second.idx, it->second.idx);
    it = coor_to_valid_idxes.find({ info.first.first - 1, info.first.second });
    if (it != coor_to_valid_idxes.end())
      gc->setNeighbors(info.second.idx, it->second.idx);
  }
  cv::Mat
      x_kernel = (cv::Mat_<float>(4, 3) <<
          -0.4, 0, 0.4, -2, 0, 2, -2, 0, 2, -0.4, 0, 0.4),
      y_kernel = (cv::Mat_<float>(3, 4) <<
          -0.4, -2, -2, -0.4, 0, 0, 0, 0, 0.4, 2, 2, 0.4),
      covered_x_mat, covered_y_mat, new_x_mat, new_y_mat;
  cv::filter2D(covered_mat, covered_x_mat, CV_32F, x_kernel, cv::Point(1, 1));
  cv::filter2D(covered_mat, covered_y_mat, CV_32F, y_kernel, cv::Point(1, 1));
  cv::filter2D(new_mat, new_x_mat, CV_32F, x_kernel, cv::Point(1, 1));
  cv::filter2D(new_mat, new_y_mat, CV_32F, y_kernel, cv::Point(1, 1));
  SmoothExtraData smooth_extra_data{
      grad_term_exp_, diff_term_low_trunc_, diff_term_high_trunc_, covered_mat,
      covered_x_mat, covered_y_mat, new_mat, new_x_mat, new_y_mat,
      valid_idx_to_coors };
  gc->setSmoothCost(SmoothCost, &smooth_extra_data);
  gc->expansion(2);
  spdlog::debug("Creating a graph cut class and run the algorithm - done");

  // Update the label raster
  spdlog::debug("Updating the label raster");
  cv::Mat new_label_mat(cv::Mat::zeros(label_mat.size(), CV_8UC1));
  for (const auto& info : coor_to_valid_idxes)
    new_label_mat.at<uint8_t>(info.first.first, info.first.second) =
        gc->whatLabel(info.second.idx) == 0 ? 100 : 200;
  label_raster_dataset = utils::CreateDatasetFromMat(
      new_label_mat, "", geotrans, spatial_ref);
  spdlog::debug("Updating the label raster - done");

  // Polygonize the label raster
  spdlog::debug("Polygonizing the new label raster");
  GDALDatasetUniquePtr polygonized_vector_dataset(memory_driver_->Create(
      "", 0, 0, 0, GDT_Unknown, nullptr));
  OGRLayer* polygonized_layer(polygonized_vector_dataset->CreateLayer(
      "", spatial_ref, wkbPolygon));
  OGRFieldDefn label_field("label", OFTInteger);
  polygonized_layer->CreateField(&label_field);
  GDALPolygonize(
      label_raster_dataset->GetRasterBand(1), nullptr, polygonized_layer, 0,
      nullptr, nullptr, nullptr);
  spdlog::debug("Polygonizing the new label raster - done");

  // Find label geometries among polygonized features
  spdlog::debug("Finding label geometries among polygonized features");
  std::pair<GIntBig, double> label0_max(-1, -1.0), label1_max(-1, -1.0);
  for (const auto& polygonized_feature : polygonized_layer) {
    GIntBig id(polygonized_feature->GetFID());
    switch (polygonized_feature->GetFieldAsInteger(0)) {
      case 100: {
        double area(
            polygonized_feature->GetGeometryRef()->toPolygon()->get_Area());
        if (area > label0_max.second)
          label0_max = { id, area };
        break;
      }
      case 200: {
        double area(
            polygonized_feature->GetGeometryRef()->toPolygon()->get_Area());
        if (area > label1_max.second)
          label1_max = { id, area };
      }
    }
  }
  label0_geometry.reset(
      polygonized_layer->GetFeature(label0_max.first)->GetGeometryRef()
      ->Boundary());
  label1_geometry.reset(
      polygonized_layer->GetFeature(label1_max.first)->GetGeometryRef()
      ->Boundary());
  spdlog::debug("Finding label geometries among polygonized features - done");

  // Eliminate internal holes in the label geometries
  spdlog::debug("Eliminating internal holes in the label geometries");
  if (label0_geometry->getGeometryType() == wkbMultiLineString)
    label0_geometry.reset(
        label0_geometry->toMultiLineString()->getGeometryRef(0)->clone());
  if (label0_geometry->getGeometryType() == wkbMultiLineString)
    label1_geometry.reset(
        label1_geometry->toMultiLineString()->getGeometryRef(0)->clone());
  spdlog::debug("Eliminating internal holes in the label geometries - done");
  spdlog::debug("Execute the graph cut mosaicking algorithm - done");
}

cv::Mat GraphCutImpl::CreateLightnessMat(GDALDataset* dataset) {
  int bands_map[3]{ 3, 2, 1 };
  cv::Mat lab_mat(utils::CreateMatFromDataset(dataset, nullptr, 3, bands_map));
  std::vector<cv::Mat> lab_mats;
  lab_mat.convertTo(lab_mat, CV_32FC3, 1.0 / 255);
  cv::cvtColor(lab_mat, lab_mat, cv::COLOR_BGR2Lab);
  cv::split(lab_mat, lab_mats);
  return lab_mats[0];
}

int GraphCutImpl::SmoothCost(
      int site1,
      int site2,
      int label1,
      int label2,
      void* extra_data) {
  auto _extra_data(static_cast<SmoothExtraData*>(extra_data));
  if (label1 == label2) return 0;
  int row1(_extra_data->valid_idx_to_coors.at(site1).first),
      col1(_extra_data->valid_idx_to_coors.at(site1).second),
      row2(_extra_data->valid_idx_to_coors.at(site2).first),
      col2(_extra_data->valid_idx_to_coors.at(site2).second),
      min_row(row1 < row2 ? row1 : row2),
      min_col(col1 < col2 ? col1 : col2),
      max_row(row1 > row2 ? row1 : row2),
      max_col(col1 > col2 ? col1 : col2);
  float gradient_term(200.0);
  if (min_row != 0 && max_row != _extra_data->covered_mat.rows - 1 &&
      min_col != 0 && max_col != _extra_data->covered_mat.cols - 1) {
    if (row1 == row2) {
      cv::Rect rect(min_col - 1, min_row - 1, 4, 3);
      cv::Mat
          box_mat1(cv::Mat(_extra_data->covered_mat, rect)),
          box_mat2(cv::Mat(_extra_data->new_mat, rect));
      if (
          cv::countNonZero(box_mat1) == 12 && 
          cv::countNonZero(box_mat2) == 12) {
        gradient_term = std::abs(
            _extra_data->covered_y_mat.at<float>(min_row, min_col) -
            _extra_data->new_y_mat.at<float>(min_row, min_col));
      }
    } else {
      cv::Rect rect(min_col - 1, min_row - 1, 3, 4);
      cv::Mat
          box_mat1(cv::Mat(_extra_data->covered_mat, rect)),
          box_mat2(cv::Mat(_extra_data->new_mat, rect));
      if (cv::countNonZero(box_mat1) == 12 &&
          cv::countNonZero(box_mat2) == 12) {
        gradient_term = std::abs(
            _extra_data->covered_x_mat.at<float>(min_row, min_col) -
            _extra_data->new_x_mat.at<float>(min_row, min_col));
      }
    }
  } 
  float diff_term(0.0);
  if (label1 == 0) {
    diff_term = std::abs(
        _extra_data->covered_mat.at<float>(row1, col1) -
        _extra_data->new_mat.at<float>(row2, col2));
  } else {
    diff_term = std::abs(
        _extra_data->new_mat.at<float>(row1, col1) - 
        _extra_data->covered_mat.at<float>(row2, col2));
  }
  return static_cast<int>(round(
      100 * std::pow(gradient_term, _extra_data->grad_exp) *
      (fmax(fmin(diff_term, _extra_data->max_diff), _extra_data->min_diff))));
}

std::shared_ptr<GraphCut> GraphCut::Create(
    double grad_term_exp,
    double diff_term_low_trunc,
    double diff_term_high_trunc,
    double tol) {
  return std::make_shared<GraphCutImpl>(
      grad_term_exp, diff_term_low_trunc, diff_term_high_trunc, tol);
}

} // mosaicking
} // rs_toolset