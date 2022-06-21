#include "graph_cut_mosaicking.h"

#include <cmath>
#include <cstdint>

#include <map>
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
    GDALDataset* covered_raster_dataset,
    GDALDataset* new_raster_dataset,
    GDALDataset* seamline_raster_dataset,
    cv::Mat& covered_mat,
    cv::Mat& new_mat,
    cv::Mat& label_mat) {
  spdlog::debug("Preparing the data");
  covered_mat = CreateLumiData(covered_raster_dataset);
  new_mat = CreateLumiData(new_raster_dataset);
  if (seamline_raster_dataset) {
    label_mat = utils::CreateMatFromDataset(seamline_raster_dataset);
  } else {
    label_mat = cv::Mat::zeros(new_mat.size(), CV_8UC1);
  }
  spdlog::debug("Preparing the data - done");
}

void GraphCutImpl::CreateSeamlineGeometries(
    GDALDataset* new_raster_dataset,
    const cv::Mat& covered_mat,
    const cv::Mat& new_mat,
    const cv::Mat& label_mat,
    GDALDatasetUniquePtr& seamline_raster_dataset,
    OGRGeometryUniquePtr& label0_geometry,
    OGRGeometryUniquePtr& label1_geometry,
    OGRLayer* seamline_layer) {
  spdlog::debug("Creating seamline geometries");
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
  spdlog::debug("Creating the maps between coors and valid indexes - done");
  auto gc(std::make_unique<GCoptimizationGeneralGraph>(valid_idx, 2));
  std::map<std::pair<int, int>, GraphCutInfo>::iterator it;
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
  spdlog::debug("Setting the data energy and neighbors - done");
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
      grad_exp_, min_diff_, max_diff_, covered_mat, covered_x_mat,
      covered_y_mat, new_mat, new_x_mat, new_y_mat, valid_idx_to_coors };
  gc->setSmoothCost(SmoothCost, &smooth_extra_data);
  spdlog::debug("Setting the smooth energy - done");
  gc->expansion(2);
  spdlog::debug("Excuting the graph cut method - done");
  cv::Mat new_label_mat(cv::Mat::zeros(label_mat.size(), CV_8UC1));
  for (const auto& info : coor_to_valid_idxes)
    new_label_mat.at<uint8_t>(info.first.first, info.first.second) = 
        gc->whatLabel(info.second.idx) == 0 ? 100 : 200;
  double geotrans[6];
  new_raster_dataset->GetGeoTransform(geotrans);
  seamline_raster_dataset = utils::CreateDatasetFromMat(
      new_label_mat, "", geotrans, 
      const_cast<OGRSpatialReference*>(new_raster_dataset->GetSpatialRef()));
  GDALDatasetUniquePtr polygonized_vector_dataset(memory_driver_->Create(
      "", 0, 0, 0, GDT_Unknown, nullptr));
  OGRLayer* polygonized_layer(polygonized_vector_dataset->CreateLayer(
      "", const_cast<OGRSpatialReference*>(new_raster_dataset->GetSpatialRef()),
      wkbPolygon));
  OGRFieldDefn label_field("label", OFTInteger);
  polygonized_layer->CreateField(&label_field);
  GDALPolygonize(
      seamline_raster_dataset->GetRasterBand(1), nullptr, polygonized_layer, 0,
      nullptr, nullptr, nullptr);
  spdlog::debug("Polygonizing the new label raster - done");
  std::pair<GIntBig, double> label0_max(-1, -1.), label1_max(-1, -1.);
  for (const auto& polygonized_feature : polygonized_layer) {
    GIntBig fid = polygonized_feature->GetFID();
    switch (polygonized_feature->GetFieldAsInteger(0)) {
      case 100: {
        double area(
            polygonized_feature->GetGeometryRef()->toPolygon()->get_Area());
        if (area > label0_max.second)
          label0_max = { fid, area };
        break;
      }
      case 200: {
        double area(
            polygonized_feature->GetGeometryRef()->toPolygon()->get_Area());
        if (area > label1_max.second)
          label1_max = { fid, area };
      }
    }
  }
  spdlog::debug("Finding label geometries - done");
  label0_geometry.reset(
      polygonized_layer->GetFeature(label0_max.first)->GetGeometryRef()
      ->clone());
  OGRGeometryUniquePtr label0_boundary(label0_geometry->Boundary());
  if (label0_boundary->getGeometryType() == wkbMultiLineString) {
    OGRGeometry* label0_outer_boundary(
        label0_boundary->toMultiLineString()->getGeometryRef(0)->clone());
    label0_geometry.reset(OGRGeometryFactory::forceToPolygon(
        label0_outer_boundary));
  }
  if (seamline_layer) {
    OGRFeatureUniquePtr label0_feature(OGRFeature::CreateFeature(
        seamline_layer->GetLayerDefn()));
    label0_feature->SetGeometry(label0_geometry.get());
    label0_feature->SetField(0, 100);
    seamline_layer->CreateFeature(label0_feature.get());
    spdlog::debug("Creating the label0 feature for the seamline layer - done");
  }
  label1_geometry.reset(
      polygonized_layer->GetFeature(label1_max.first)->GetGeometryRef()
      ->clone());
  OGRGeometryUniquePtr label1_boundary(label1_geometry->Boundary());
  if (label1_boundary->getGeometryType() == wkbMultiLineString) {
    OGRGeometry* label1_outer_boundary(
        label1_boundary->toMultiLineString()->getGeometryRef(0)->clone());
    label1_geometry.reset(OGRGeometryFactory::forceToPolygon(
        label1_outer_boundary));
  }
  if (seamline_layer) {
    OGRFeatureUniquePtr label1_feature(OGRFeature::CreateFeature(
        seamline_layer->GetLayerDefn()));
    label1_feature->SetGeometry(label1_geometry.get());
    label1_feature->SetField(0, 200);
    seamline_layer->CreateFeature(label1_feature.get());
    spdlog::debug("Creating the label1 feature for the seamline layer - done");
  }
  spdlog::debug("Creating seamline geometries - done");
}

cv::Mat GraphCutImpl::CreateLumiData(GDALDataset* dataset) {
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
  SmoothExtraData* _extra_data = static_cast<SmoothExtraData*>(
      extra_data);
  if (label1 == label2) return 0;
  int row1(_extra_data->valid_idx_to_coors.at(site1).first),
      col1(_extra_data->valid_idx_to_coors.at(site1).second),
      row2(_extra_data->valid_idx_to_coors.at(site2).first),
      col2(_extra_data->valid_idx_to_coors.at(site2).second),
      min_row(row1 < row2 ? row1 : row2),
      min_col(col1 < col2 ? col1 : col2),
      max_row(row1 > row2 ? row1 : row2),
      max_col(col1 > col2 ? col1 : col2);
  float gradient_term(200.);
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
  float diff_term(0.);
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
    double grad_exp,
    double diff_min,
    double diff_max) {
  return std::make_shared<GraphCutImpl>(grad_exp, diff_min, diff_max);
}

} // mosaicking 
} // rs_toolset