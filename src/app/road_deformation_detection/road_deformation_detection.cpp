#include "road_deformation_detection.h"

#include <cmath>
#include <cstring>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gdal_alg.h>
#include <gdal_priv.h>
#include <nlohmann/json.hpp>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace road_deformation_detection {

std::shared_ptr<RoadDeformationDetection> RoadDeformationDetection::Create(
    const std::string& road_path,
    const std::string& road_layer_name,
    const std::string& dem_path,
    double ideal_sample_space_in_pixels,
    double threshold,
    int min_continual_count) {
  auto road_dataset(GDALDataset::Open(
      road_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
  if (!road_dataset) {
    spdlog::warn("Opening {} failed", road_path);
    return nullptr;
  }
  OGRLayer* road_layer(nullptr);
  if (road_layer_name.empty()) {
    road_layer = road_dataset->GetLayer(0);
  } else {
    road_layer = road_dataset->GetLayerByName(road_layer_name.c_str());
    if (!road_layer) {
      spdlog::warn("Getting road layer \"{}\" failed", road_layer_name);
      return nullptr;
    }
  }
  auto dem_dataset(GDALDataset::Open(
      dem_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!dem_dataset) {
    spdlog::warn("Opening {} failed", dem_path);
    return nullptr;
  }
  double dem_geotrans[6];
  dem_dataset->GetGeoTransform(dem_geotrans);
  auto dem_nodata_value(static_cast<float>(
      dem_dataset->GetRasterBand(1)->GetNoDataValue()));
  return std::shared_ptr<RoadDeformationDetection>(
      new RoadDeformationDetection(
          road_dataset, road_layer, dem_dataset, dem_geotrans,
          dem_nodata_value, ideal_sample_space_in_pixels, threshold,
          min_continual_count));
}

RoadDeformationDetection::~RoadDeformationDetection() {
  GDALClose(road_dataset_);
  GDALClose(dem_dataset_);
}

GDALDatasetUniquePtr RoadDeformationDetection::Run(
    const std::string& dom_path,
    const std::string& rpb_path,
    const std::string& output_path) {
  spdlog::info(
      "Running a road deformation detection task with\n"
      " - DOM path: {}\n"
      " - RPB path: {}\n"
      " - Output path: {}", dom_path, rpb_path, output_path);
  GDALDatasetUniquePtr dom_dataset(GDALDataset::Open(
      dom_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!dom_dataset) {
    spdlog::warn("Opening {} failed", dom_path);
    return nullptr;
  }
  double dom_geotrans[6];
  if (dom_dataset->GetGeoTransform(dom_geotrans) != CE_None ||
      !dom_dataset->GetSpatialRef()) {
    spdlog::warn(
        "{} does not have the geotransform or the spatial reference. "
        "Please check whether the raster is DOM", dom_path);
    return nullptr;
  }
  auto rpc_info_text(utils::LoadRPBFile(rpb_path));
  GDALRPCInfo rpc_info;
  if (!GDALExtractRPCInfo(rpc_info_text.get(), &rpc_info)) {
    spdlog::warn("Extracting RPC information failed");
    return nullptr;
  }
  auto driver(utils::GetVectorDriverByPath(output_path));
  if (!driver) {
    spdlog::warn("Finding the driver for {} failed", output_path);
    return nullptr;
  }
  GDALDatasetUniquePtr output_dataset(driver->Create(
      output_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
  if (!output_dataset) {
    spdlog::warn("Creating {} failed", output_path);
    return nullptr;
  }

  auto output_layer(output_dataset->CreateLayer(
      "", const_cast<OGRSpatialReference*>(dom_dataset->GetSpatialRef()),
      wkbPoint));
  OGRFieldDefn normal_bias_ratio_field("nbr", OFTReal);
  output_layer->CreateField(&normal_bias_ratio_field);
  OGRGeometryUniquePtr border(utils::CreateBorder(dom_dataset.get()));
  road_layer_->SetSpatialFilter(border.get());
  int range[4]{
      static_cast<int>(floor(
          (dom_geotrans[0] - dem_geotrans_[0]) / dem_geotrans_[1])),
      static_cast<int>(floor(
          (dom_geotrans[3] - dem_geotrans_[3]) / dem_geotrans_[5])),
      static_cast<int>(ceil(
          (dom_geotrans[0] - dem_geotrans_[0] + dom_dataset->GetRasterXSize() *
              dom_geotrans[1]) / dem_geotrans_[1])),
      static_cast<int>(ceil(
          (dom_geotrans[3] - dem_geotrans_[3] + dom_dataset->GetRasterYSize() *
              dom_geotrans[5]) / dem_geotrans_[5])) };
  auto dem_mat(utils::CreateMatFromDataset(dem_dataset_, range));
  double geotrans[6];
  memcpy(geotrans, dem_geotrans_, 6 * sizeof(double));
  geotrans[0] += range[0] * geotrans[1];
  geotrans[3] += range[1] * geotrans[5];
  double hori_angle, vert_shift;
  CalcLightAngles(rpc_info, hori_angle, vert_shift);
  for (const auto& feature : road_layer_) {
    auto geometry(feature->GetGeometryRef());
    switch (geometry->getGeometryType()) {
      case wkbLineString:
      case wkbLineString25D: {
        DetectLineString(
            geometry->toLineString(), border.get(), dem_mat, geotrans,
            hori_angle, vert_shift, output_layer);
        break;
      }
      case wkbMultiLineString:
      case wkbMultiLineString25D: {
        for (const auto& line_string : geometry->toMultiLineString())
          DetectLineString(
              line_string, border.get(), dem_mat, geotrans, hori_angle,
              vert_shift, output_layer);
        break;
      }
    }
  }
  spdlog::info("Running a road deformation detection task - done");
  return output_dataset;
}

RoadDeformationDetection::RoadDeformationDetection(
    GDALDataset* road_dataset,
    OGRLayer* road_layer,
    GDALDataset* dem_dataset,
    double* dem_geotrans,
    float dem_nodata_value,
    double ideal_sample_space_in_pixels,
    double threshold,
    int min_continual_count)
    : road_dataset_(road_dataset),
      road_layer_(road_layer),
      dem_dataset_(dem_dataset),
      dem_nodata_value_(dem_nodata_value),
      ideal_sample_space_in_pixels_(ideal_sample_space_in_pixels),
      threshold_(threshold),
      min_continual_count_(min_continual_count) {
  memcpy(dem_geotrans_, dem_geotrans, 6 * sizeof(double));
  spdlog::info(
      "Creating a road deformation detection with\n"
      " - Ideal sample space in pixels: {}\n"
      " - Threshold; {}\n"
      " - Minimum continual count: {}",
      ideal_sample_space_in_pixels, threshold, min_continual_count);
}

void RoadDeformationDetection::CalcLightAngles(
    const GDALRPCInfo& rpc_info,
    double& hori_angle,
    double& vert_shift) {
  spdlog::debug("Calculating the light angles from the RPC information");
  utils::RPCTransPtr forward_trans_arg(
      GDALCreateRPCTransformer(&rpc_info, false, 0.0, nullptr),
      [](void* p) { GDALDestroyRPCTransformer(p); });
  utils::RPCTransPtr back_trans_arg(
      GDALCreateRPCTransformer(&rpc_info, true, 0.0, nullptr),
      [](void* p) { GDALDestroyRPCTransformer(p); });
  double
      x(rpc_info.dfLONG_OFF),
      y(rpc_info.dfLAT_OFF),
      z(rpc_info.dfHEIGHT_OFF);
  int success;
  GDALRPCTransform(forward_trans_arg.get(), true, 1, &x, &y, &z, &success);
  GDALRPCTransform(back_trans_arg.get(), true, 1, &x, &y, &(++z), &success);
  double delta_x(x - rpc_info.dfLONG_OFF), delta_y(y - rpc_info.dfLAT_OFF);
  hori_angle = atan2(delta_y, delta_x);
  vert_shift = sqrt(delta_x * delta_x + delta_y * delta_y);
  spdlog::info("Calculating the light angles from the RPC information - done");
}

void RoadDeformationDetection::DetectLineString(
    OGRLineString* line_string,
    OGRGeometry* border,
    cv::Mat dem_mat,
    double* geotrans,
    double hori_angle,
    double vert_shift,
    OGRLayer* layer) {
  OGRPoint start_point, end_point;
  nlohmann::json json;
  std::vector<std::pair<double, double>> coors;
  for (int i = 1; i < line_string->getNumPoints(); i++) {
    line_string->getPoint(i - 1, &start_point);
    line_string->getPoint(i, &end_point);
    json["type"] = "LineString";
    json["coordinates"] = {
        { start_point.getX(), start_point.getY() },
        { end_point.getX(), end_point.getY() } };
    OGRGeometryUniquePtr geometry(OGRGeometryFactory::createFromGeoJson(
        json.dump().c_str()));
    LineSampleInfo info;
    if (border->Contains(geometry.get())) {
      SampleBetweenPoints(
          start_point, end_point, dem_mat, geotrans, hori_angle, vert_shift,
          info);
    } else {
      geometry.reset(geometry->Intersection(border));
      switch (geometry->getGeometryType()) {
        case wkbLineString: {
          json = nlohmann::json::parse(geometry->exportToJson());
          json.at("coordinates").get_to(coors);
          start_point = { coors.front().first, coors.front().second };
          end_point = { coors.back().first, coors.back().second };
          SampleBetweenPoints(
              start_point, end_point, dem_mat, geotrans, hori_angle,
              vert_shift, info);
          break;
        }
        case wkbMultiLineString: {
          for (const auto& line_string : geometry->toMultiLineString()) {
            json = nlohmann::json::parse(line_string->exportToJson());
            json.at("coordinates").get_to(coors);
            start_point = { coors.front().first, coors.front().second };
            end_point = { coors.back().first, coors.back().second };
            SampleBetweenPoints(
                start_point, end_point, dem_mat, geotrans, hori_angle,
                vert_shift, info);
          }
          break;
        }
      }
    }
    for (int j = 0; j < info.points.size(); j++) {
      OGRFeatureUniquePtr feature(OGRFeature::CreateFeature(
          layer->GetLayerDefn()));
      feature->SetGeometryDirectly(info.points[j]);
      feature->SetField("nbr", info.normal_bias_ratios[j]);
      layer->CreateFeature(feature.get());
    }
  }
}

void RoadDeformationDetection::SampleBetweenPoints(
    const OGRPoint& point1,
    const OGRPoint& point2,
    cv::Mat dem_mat,
    double* geotrans,
    double hori_angle,
    double vert_shift,
    LineSampleInfo& info) {
  double total_dist(point1.Distance(&point2));
  int spaces_count(static_cast<int>(ceil(
      total_dist * ideal_sample_space_in_pixels_ / dem_geotrans_[1])));
  double
      x_space((point2.getX() - point1.getX()) / spaces_count),
      y_space((point2.getY() - point1.getY()) / spaces_count),
      final_shift(
          vert_shift * sin(atan2(y_space, x_space) - hori_angle) /
          total_dist * spaces_count);
  std::vector<double> xs{ point1.getX() }, ys{ point1.getY() };
  xs.reserve(spaces_count + 1);
  ys.reserve(spaces_count + 1);
  for (int i = 0; i < spaces_count; i++) {
    xs.push_back(xs.back() + x_space);
    ys.push_back(ys.back() + y_space);
  }
  std::vector<float> elevations(spaces_count + 1);
  utils::CalcInterpolations(
      dem_mat, geotrans, spaces_count + 1, xs.data(), ys.data(),
      elevations.data(), false, dem_nodata_value_);
  LineSampleInfo cur_info;
  for (int i = 0; i < spaces_count;
      i++, cur_info.normal_bias_ratios.clear(), cur_info.points.clear()) {
    if (double normal_bias_ratio(abs(
            (elevations[i + 1] - elevations[i]) * final_shift));
        normal_bias_ratio > threshold_) {
      cur_info.points.push_back(new OGRPoint(xs[i], ys[i]));
      cur_info.normal_bias_ratios.push_back(normal_bias_ratio);
      if (i != spaces_count - 1) continue;
    }
    if (cur_info.points.size() >= min_continual_count_) {
      info.points.insert(
          info.points.end(), cur_info.points.begin(), cur_info.points.end());
      info.normal_bias_ratios.insert(
          info.normal_bias_ratios.end(), cur_info.normal_bias_ratios.begin(),
          cur_info.normal_bias_ratios.end());
    } else {
      for (auto& point : cur_info.points)
        delete point;
    }
  }

}

} // namespace road_deformation_detection
} // namespace rs_toolset