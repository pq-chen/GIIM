#include <rs-toolset/utils.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <list>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gdalwarper.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace fs = std::filesystem;

namespace rs_toolset {
namespace utils {

OGRGeometryUniquePtr ApplyGeotransOnGeometry(
    OGRGeometry* geometry,
    double* geotrans) {
  OGRGeometryUniquePtr new_geometry(geometry->clone());
  switch (double coor[2]; new_geometry->getGeometryType()) {
    case wkbLineString: {
      for (auto& point : new_geometry->toLineString()) {
        GDALApplyGeoTransform(
            geotrans, point.getX(), point.getY(), coor, coor + 1);
        point.setX(coor[0]);
        point.setY(coor[1]);
      }
      break;
    }
    case wkbPolygon: {
      for (auto& linear_ring : new_geometry->toPolygon()) {
        for (auto& point : linear_ring) {
          GDALApplyGeoTransform(
              geotrans, point.getX(), point.getY(), coor, coor + 1);
          point.setX(coor[0]);
          point.setY(coor[1]);
        }
      }
      break;
    }
    case wkbMultiLineString: {
      for (auto& line_string : new_geometry->toMultiLineString()) {
        for (auto& point : line_string) {
          GDALApplyGeoTransform(
              geotrans, point.getX(), point.getY(), coor, coor + 1);
          point.setX(coor[0]);
          point.setY(coor[1]);
        }
      }
      break;
    }
    case wkbMultiPolygon: {
      for (auto& polygon : new_geometry->toMultiPolygon()) {
        for (auto& linear_ring : polygon) {
          for (auto& point : linear_ring) {
            GDALApplyGeoTransform(
                geotrans, point.getX(), point.getY(), coor, coor + 1);
            point.setX(coor[0]);
            point.setY(coor[1]);
          }
        }
      }
      break;
    }
    default: {
      return nullptr;
    }
  }
  return new_geometry;
}

OGRGeometryUniquePtr CreateBorder(
    GDALDataset* dataset,
    OGRSpatialReference* spatial_ref,
    int downsample_factor,
    double buffer,
    double tol) {
  // Polygonize the mask raster dataset to the mask vector dataset
  auto mask_raster_dataset(CreateMaskDataset(
      dataset->GetRasterBand(1), downsample_factor));
  auto driver(GetGDALDriverManager()->GetDriverByName("Memory"));
  GDALDatasetUniquePtr mask_vector_dataset(driver->Create(
      "", 0, 0, 0, GDT_Unknown, nullptr));
  auto mask_layer(mask_vector_dataset->CreateLayer(
      "", const_cast<OGRSpatialReference*>(dataset->GetSpatialRef()),
      wkbPolygon));
  OGRFieldDefn value_field("value", OFTInteger);
  mask_layer->CreateField(&value_field);
  GDALPolygonize(
      mask_raster_dataset->GetRasterBand(1), nullptr, mask_layer, 0, nullptr,
      nullptr, nullptr);

  // Find the biggeest geometry among all features with label 255
  std::pair max_idx(-1, -1.0);
  for (const auto& feature : mask_layer) {
    if (feature->GetFieldAsInteger(0) == 255) {
      if (auto area(feature->GetGeometryRef()->toPolygon()->get_Area());
          area > max_idx.second) {
        max_idx = {static_cast<int>(feature->GetFID()), area};
      }
    }
  }
  OGRFeatureUniquePtr feature(mask_layer->GetFeature(max_idx.first));
  OGRGeometryUniquePtr border(feature->StealGeometry());

  // Buffer and simplify the border
  double geotrans[6];
  dataset->GetGeoTransform(geotrans);
  border.reset(border->Simplify(downsample_factor * tol * geotrans[1]));
  border.reset(border->Buffer(downsample_factor * buffer * geotrans[1]));
  border.reset(border->Simplify(downsample_factor * tol * geotrans[1]));

  if (spatial_ref) border->transformTo(spatial_ref);
  return border;
}

GDALDatasetUniquePtr CreateDatasetFromMat(
    const cv::Mat& mat,
    const std::string& path,
    double* geotrans,
    OGRSpatialReference* spatial_ref,
    int bands_count,
    int* bands_map) {
  auto driver(GetRasterDriverByPath(path));
  if (!driver) return nullptr;

  GDALDataType data_type;
  switch(mat.depth()) {
    case CV_8U: {
      data_type = GDT_Byte;
      break;
    }
    case CV_16U: {
      data_type = GDT_UInt16;
      break;
    }
    case CV_16S: {
      data_type = GDT_Int16;
      break;
    }
    case CV_32S: {
      data_type = GDT_Int32;
      break;
    }
    case CV_32F: {
      data_type = GDT_Float32;
      break;
    }
    case CV_64F: {
      data_type = GDT_Float64;
      break;
    }
    default: {
      return nullptr;
    }
  }
  auto x_size(mat.cols),
      y_size(mat.rows),
      bytes_count(GDALGetDataTypeSizeBytes(data_type));
  if (!bands_count) bands_count = mat.channels();
  GDALDatasetUniquePtr dataset(driver->Create(
      path.c_str(), x_size, y_size, bands_count, data_type, nullptr));
  if (geotrans) dataset->SetGeoTransform(geotrans);
  if (spatial_ref) dataset->SetSpatialRef(spatial_ref);
  dataset->RasterIO(
      GF_Write, 0, 0, x_size, y_size, mat.data, x_size, y_size, data_type,
      bands_count, bands_map, bytes_count * mat.channels(),
      bytes_count * mat.channels() * x_size, bytes_count);
  return dataset;
}

OGRGeometryUniquePtr CreateGeometryFromEnve(const OGREnvelope& enve) {
  auto linear_ring(new OGRLinearRing);
  linear_ring->addPoint(enve.MinX, enve.MaxY);
  linear_ring->addPoint(enve.MinX, enve.MinY);
  linear_ring->addPoint(enve.MaxX, enve.MinY);
  linear_ring->addPoint(enve.MaxX, enve.MaxY);
  linear_ring->addPoint(enve.MinX, enve.MaxY);
  auto polygon(new OGRPolygon);
  polygon->addRingDirectly(linear_ring);
  return OGRGeometryUniquePtr(polygon);
}

std::vector<cv::Mat> CreateHists(
    const cv::Mat& source_mat,
    const cv::Mat& mask_mat) {
  if (source_mat.depth() != CV_8U && source_mat.depth() != CV_16U) return {};

  int bands_count(source_mat.channels()),
      hist_size[1]{source_mat.depth() == CV_8U ? 256 : 65536};
  float _ranges[2]{0.0, static_cast<float>(hist_size[0])};
  const float* ranges[1]{_ranges};
  std::vector<cv::Mat> output_mats(bands_count);
#pragma omp parallel for schedule(static, bands_count)
  for (int b(0); b < bands_count; ++b) {
    const int bands[1]{b};
    cv::calcHist(
        &source_mat, 1, bands, mask_mat, output_mats[b], 1, hist_size, ranges);
  }
  return output_mats;
}

cv::Mat CreateHistMatchingLut(
    const std::vector<cv::Mat>& source_hist_mats,
    const std::vector<cv::Mat>& target_hist_mats) {
  if (source_hist_mats.size() != target_hist_mats.size()) return cv::Mat();

  auto bands_count(static_cast<int>(source_hist_mats.size())),
      depth(target_hist_mats[0].rows <= 256 ? CV_8U : CV_16U);
  cv::Mat lut_mat(1, source_hist_mats[0].rows, CV_MAKETYPE(depth, bands_count));
#pragma omp parallel for schedule(static, bands_count)
  for (int b(0); b < bands_count; ++b) {
    // Normalize histogram mats to create cumulative distribution function(CDF) mats
    cv::Mat 
        _source_hist_mat(source_hist_mats[b].clone()),
        _target_hist_mat(target_hist_mats[b].clone());
    _source_hist_mat.at<float>(0) = 0.0;
    _target_hist_mat.at<float>(0) = 0.0;
    cv::normalize(_source_hist_mat, _source_hist_mat, 1.0, 0.0, cv::NORM_L1);
    cv::normalize(_target_hist_mat, _target_hist_mat, 1.0, 0.0, cv::NORM_L1);
    auto source_cdf(std::make_unique<float[]>(_source_hist_mat.rows)),
        target_cdf(std::make_unique<float[]>(_target_hist_mat.rows));
    source_cdf[0] = 0.0;
    target_cdf[0] = 0.0;
    for (int i(1); i < _source_hist_mat.rows; ++i)
      source_cdf[i] = source_cdf[i - 1] + _source_hist_mat.at<float>(i);
    for (int i(1); i < _target_hist_mat.rows; ++i)
      target_cdf[i] = target_cdf[i - 1] + _target_hist_mat.at<float>(i);

    // For each source CDF value, find the closest target CDF value and its index
    for (int i(0), idx(0); i < _source_hist_mat.rows; ++i) {
      if (source_cdf[i] > target_cdf[idx]) {
        while (idx != _target_hist_mat.rows - 1 &&
            (source_cdf[i] > target_cdf[idx + 1])) {
          ++idx;
        }
        if (idx != _target_hist_mat.rows - 1 &&
            target_cdf[idx + 1] - source_cdf[i] <
                source_cdf[i] - target_cdf[idx]) {
          ++idx;
        }
      }
      if (depth == CV_8U) {
        *(lut_mat.ptr<uint8_t>(0) + bands_count * i + b) =
            static_cast<uint8_t>(idx);
      } else {
        *(lut_mat.ptr<uint16_t>(0) + bands_count * i + b) =
            static_cast<uint16_t>(idx);
      }
    }
  }
  return lut_mat;
}

template <typename T>
std::unique_ptr<uint8_t[]> CreateMaskDatasetImpl(
    GDALRasterBand* raster_band,
    int downsample_factor) {
  auto data_type(raster_band->GetRasterDataType());
  auto x_size(raster_band->GetXSize()),
      y_size(raster_band->GetYSize()),
      new_x_size(raster_band->GetXSize() / downsample_factor),
      new_y_size(raster_band->GetYSize() / downsample_factor),
      bytes_count(GDALGetDataTypeSizeBytes(data_type));
  auto size(static_cast<uint64_t>(new_x_size) * new_y_size);
  auto nodata_value(static_cast<T>(raster_band->GetNoDataValue()));
  auto source_data(std::make_unique<T[]>(size));
  raster_band->RasterIO(
      GF_Read, 0, 0, x_size, y_size, source_data.get(), new_x_size, new_y_size,
      data_type, bytes_count, bytes_count * new_x_size);
  auto mask_data(std::make_unique<uint8_t[]>(size));
  memset(mask_data.get(), 0, size * sizeof(uint8_t));
#pragma omp parallel for schedule(dynamic)
  for (int row(0); row < new_y_size; ++row) {
    auto source_ptr(source_data.get() + row * new_x_size);
    auto mask_ptr(mask_data.get() + row * new_x_size);
    for (int col(0); col < new_x_size; ++col, ++source_ptr, ++mask_ptr)
      if (*source_ptr != nodata_value) *mask_ptr = 255;
  }
  return mask_data;
}

GDALDatasetUniquePtr CreateMaskDataset(
    GDALRasterBand* raster_band,
    int downsample_factor) {
  std::unique_ptr<uint8_t[]> mask_data(nullptr);
  switch(raster_band->GetRasterDataType()) {
    case GDT_Byte: {
      mask_data = CreateMaskDatasetImpl<uint8_t>(
          raster_band, downsample_factor);
      break;
    }
    case GDT_UInt16: {
      mask_data = CreateMaskDatasetImpl<uint16_t>(
          raster_band, downsample_factor);
      break;
    }
    case GDT_Int16: {
      mask_data = CreateMaskDatasetImpl<int16_t>(
          raster_band, downsample_factor);
      break;
    }
    case GDT_UInt32: {
      mask_data = CreateMaskDatasetImpl<uint32_t>(
          raster_band, downsample_factor);
      break;
    }
    case GDT_Int32: {
      mask_data = CreateMaskDatasetImpl<int32_t>(
          raster_band, downsample_factor);
      break;
    }
    case GDT_Float32: {
      mask_data = CreateMaskDatasetImpl<float>(raster_band, downsample_factor);
      break;
    }
    case GDT_Float64: {
      mask_data = CreateMaskDatasetImpl<double>(raster_band, downsample_factor);
      break;
    }
    default: {
      return nullptr;
    }
  }
  auto new_x_size(raster_band->GetXSize() / downsample_factor),
      new_y_size(raster_band->GetYSize() / downsample_factor);
  auto driver(GetGDALDriverManager()->GetDriverByName("MEM"));
  GDALDatasetUniquePtr dataset(driver->Create(
      "", new_x_size, new_y_size, 1, GDT_Byte, nullptr));
  double geotrans[6];
  raster_band->GetDataset()->GetGeoTransform(geotrans);
  geotrans[1] *= downsample_factor;
  geotrans[5] *= downsample_factor;
  dataset->SetGeoTransform(geotrans);
  dataset->SetSpatialRef(raster_band->GetDataset()->GetSpatialRef());
  dataset->GetRasterBand(1)->RasterIO(
      GF_Write, 0, 0, new_x_size, new_y_size, mask_data.get(), new_x_size,
      new_y_size, GDT_Byte, 1, new_x_size);
  return dataset;
}

cv::Mat CreateMatFromDataset(
    GDALDataset* dataset,
    int* range,
    int bands_count,
    int* bands_map) {
  auto data_type(dataset->GetRasterBand(1)->GetRasterDataType());
  int bytes_count, depth;
  switch(data_type) {
    case GDT_Byte: {
      bytes_count = GDALGetDataTypeSizeBytes(GDT_Byte);
      depth = CV_8U;
      break;
    }
    case GDT_UInt16: {
      bytes_count = GDALGetDataTypeSizeBytes(GDT_UInt16);
      depth = CV_16U;
      break;
    }
    case GDT_Int16: {
      bytes_count = GDALGetDataTypeSizeBytes(GDT_Int16);
      depth = CV_16S;
      break;
    }
    case GDT_Int32: {
      bytes_count = GDALGetDataTypeSizeBytes(GDT_Int32);
      depth = CV_32S;
      break;
    }
    case GDT_Float32: {
      bytes_count = GDALGetDataTypeSizeBytes(GDT_Float32);
      depth = CV_32F;
      break;
    }
    case GDT_Float64: {
      bytes_count = GDALGetDataTypeSizeBytes(GDT_Float64);
      depth = CV_64F;
      break;
    }
    default: {
      return cv::Mat();
    }
  }
  int full_range[4]{
      0, 0, dataset->GetRasterXSize() , dataset->GetRasterYSize()};
  if (!range) range = full_range;
  if (!bands_count) bands_count = dataset->GetRasterCount();
  cv::Mat mat(range[3], range[2], CV_MAKETYPE(depth, bands_count));
  dataset->RasterIO(
      GF_Read, range[0], range[1], range[2], range[3], mat.data, range[2],
      range[3], data_type, bands_count, bands_map,
      bytes_count * dataset->GetRasterCount(), 
      bytes_count * dataset->GetRasterCount() * range[2], bytes_count);
  return mat;
}

bool JointMultiLineString(OGRGeometryUniquePtr& geometry) {
  if (geometry->getGeometryType() != wkbMultiLineString) return false;

  auto _geometry(geometry->toMultiLineString());
  auto middle_it(_geometry->begin());
  OGRPoint point1, point2;
  do {
    (*middle_it++)->EndPoint(&point1);
    if (middle_it == _geometry->end()) {
      middle_it = _geometry->begin();
      break;
    }
    (*middle_it)->StartPoint(&point2);
  } while (point1 == point2);
  auto line_string(new OGRLineString);
  line_string->addPoint(&*(*middle_it)->begin());
  for (auto it(middle_it); it != _geometry->end(); ++it)
    line_string->addSubLineString(*it, 1, -1);
  for (auto it(_geometry->begin()); it != middle_it; ++it)
    line_string->addSubLineString(*it, 1, -1);
  geometry.reset(line_string);
  return true;
}

void CreateRasterPyra(
    GDALDataset* dataset,
    const std::string& compress_method,
    bool clean) {
  CPLSetConfigOption("COMPRESS_OVERVIEW", compress_method.c_str());
  CPLSetConfigOption("GDAL_OVR_CHUNKYSIZE", "4096");
  if (clean)
    dataset->BuildOverviews("", 0, nullptr, 0, nullptr, nullptr, nullptr);
  if (dataset->GetRasterBand(1)->GetOverviewCount() != 0) return;

  spdlog::debug("Creating raster pyramids");
  int x_size(dataset->GetRasterXSize()),
      y_size(dataset->GetRasterYSize()),
      downsample_factor(1),
      min_pixels_count(256),  // The width and the height of the last pyramid are not less than 256
      overviews_count(0),
      overviews[16]{0};
  while (DIV_ROUND_UP(x_size, downsample_factor) > min_pixels_count ||
      DIV_ROUND_UP(y_size, downsample_factor) > min_pixels_count) {
    downsample_factor *= 2;
    overviews[overviews_count++] = downsample_factor;
  }
  dataset->BuildOverviews(
      "NEAREST", overviews_count, overviews, 0, nullptr, nullptr, nullptr);
  spdlog::info("Creating raster pyramids - done");
}

void CreateRange(
    int block_idx,
    int block_cols_count,
    int block_rows_count,
    int block_x_size,
    int block_y_size,
    int last_block_x_size,
    int last_block_y_size,
    int(&range)[4]) {
  auto block_row(block_idx / block_cols_count),
      block_col(block_idx % block_cols_count);
  range[0] = block_col * block_x_size;
  range[1] = block_row * block_y_size;
  range[2] =
      block_col != (block_cols_count - 1) ? block_x_size : last_block_x_size;
  range[3] =
      block_row != (block_rows_count - 1) ? block_y_size : last_block_y_size;
}

RPCTransPtr CreateRPCTrans(
    const std::string& path,
    RPCSource rpc_source) {
  GDALDatasetUniquePtr dataset(GDALDataset::Open(
      path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (!dataset) return RPCTransPtr(nullptr, nullptr);

  char** rpc_info_text(nullptr);
  switch (rpc_source) {
    case RPCSource::kInternal: {
      rpc_info_text = dataset->GetMetadata("RPC");
      break;
    }
    case RPCSource::kRPBFile: {
      rpc_info_text = LoadRPBFile(path.substr(0, path.rfind('.')) + ".rpb");
    }
  }
  if (GDALRPCInfo rpc_info; !GDALExtractRPCInfo(rpc_info_text, &rpc_info)) {
    CSLDestroy(rpc_info_text);
    return RPCTransPtr(nullptr, nullptr);
  } else {
    return RPCTransPtr(
        GDALCreateRPCTransformer(&rpc_info, false, 0.0, nullptr),
        [](void* p) { GDALDestroyRPCTransformer(p); });
  }
}

int DateMatching(const std::string& string) {
  std::regex reg(
      "20[0-9]{2}((0[1-9]|1[0-2])"
      "(0[1-9]|1[0-9]|2[0-8])|(0[13-9]|1[0-2])(29|30)|(0[13578]|1[02])31)");
  std::smatch result;
  if (std::regex_search(string, result, reg)) {
    return std::stoi(result.str());
  } else {
    return 0;
  }
}

std::string GetDate() {
  auto now(std::chrono::system_clock::to_time_t(
      std::chrono::system_clock::now()));
  tm _tm;
  localtime_s(&_tm, &now);
  char buf[64];
  std::strftime(buf, sizeof(buf), "%A_%d_%B_%Y_%Hh_%Mm_%Ss", &_tm);
  return buf;
}

GDALDriver* GetRasterDriverByPath(const std::string& path) {
  if (path.empty()) return GetGDALDriverManager()->GetDriverByName("MEM");
  std::string driver_name;
  if (std::regex_match(path, std::regex(".*\\.(tif|tiff|TIF|TIFF)"))) {
    driver_name = "GTiff";
  } else if (std::regex_match(path, std::regex(".*\\.img"))) {
    driver_name = "HFA";
  } else {
    return nullptr;
  }
  return GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
}

GDALDriver* GetVectorDriverByPath(const std::string& path) {
  if (path.empty()) return GetGDALDriverManager()->GetDriverByName("Memory");
  std::string driver_name;
  if (std::regex_match(path, std::regex(".*\\.shp"))) {
    driver_name = "Esri Shapefile";
  } else {
    return nullptr;
  }
  return GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
}

void GraftSeamline(
    OGRGeometryUniquePtr& seamline,
    OGRGeometryUniquePtr& geometry1,
    OGRGeometryUniquePtr& geometry2,
    double tol) {
  if (tol) seamline.reset(seamline->Simplify(tol));

  std::list<OGRPoint> seamline_points;
  for (const auto& point : seamline->toLineString())
    seamline_points.push_back(point);
  std::vector<OGRPoint> old_points;
  std::vector<OGRPoint>::iterator it1, it2;
  OGRLinearRing* linear_ring(nullptr);

  // Graft the seamline to the geometry1
  spdlog::debug("Grafting the seamline to the geometry1");
  auto _geometry1(geometry1->toPolygon());
  for (const auto& point : _geometry1->getExteriorRing())
    old_points.push_back(point);
  it1 = std::find(
      old_points.begin(), old_points.end(), seamline_points.front());
  while (it1 == old_points.end()) {
    seamline_points.pop_front();
    it1 = std::find(
        old_points.begin(), old_points.end(), seamline_points.front());
  }
  it2 = std::find(old_points.begin(), old_points.end(), seamline_points.back());
  while (it2 == old_points.end()) {
    seamline_points.pop_back();
    it2 = std::find(
        old_points.begin(), old_points.end(), seamline_points.back());
  }
  linear_ring = new OGRLinearRing;
  if (it2 > it1) {
    for (auto it(old_points.begin()); it != it1; ++it)
      linear_ring->addPoint(&*it);
    for (auto it(seamline_points.begin()); it != seamline_points.end(); ++it)
      linear_ring->addPoint(&*it);
    for (auto it(it2 + 1); it != old_points.end(); ++it)
      linear_ring->addPoint(&*it);
  } else {
    for (auto it(old_points.begin()); it != it2; ++it)
      linear_ring->addPoint(&*it);
    for (auto it(seamline_points.rbegin()); it != seamline_points.rend(); ++it)
      linear_ring->addPoint(&*it);
    for (auto it(it1 + 1); it != old_points.end(); ++it)
      linear_ring->addPoint(&*it);
  }
  _geometry1->removeRing(0);
  _geometry1->addRingDirectly(linear_ring);
  spdlog::debug("Grafting the seamline to the geometry1 - done");

  // Graft the seamline to the geometry2
  spdlog::debug("Grafting the seamline to the geometry2");
  auto _geometry2(geometry2->toPolygon());
  old_points.clear();
  for (const auto& point : _geometry2->getExteriorRing())
    old_points.push_back(point);
  it1 = std::find(
      old_points.begin(), old_points.end(), seamline_points.front());
  while (it1 == old_points.end()) {
    seamline_points.pop_front();
    it1 = std::find(
        old_points.begin(), old_points.end(), seamline_points.front());
  }
  it2 = std::find(old_points.begin(), old_points.end(), seamline_points.back());
  while (it2 == old_points.end()) {
    seamline_points.pop_back();
    it2 = std::find(
        old_points.begin(), old_points.end(), seamline_points.back());
  }
  linear_ring = new OGRLinearRing;
  if (it2 > it1) {
    for (auto it(old_points.begin()); it != it1; ++it)
      linear_ring->addPoint(&*it);
    for (auto it(seamline_points.begin()); it != seamline_points.end(); ++it)
      linear_ring->addPoint(&*it);
    for (auto it(it2 + 1); it != old_points.end(); ++it)
      linear_ring->addPoint(&*it);
  } else {
    for (auto it(old_points.begin()); it != it2; ++it)
      linear_ring->addPoint(&*it);
    for (auto it(seamline_points.rbegin()); it != seamline_points.rend(); ++it)
      linear_ring->addPoint(&*it);
    for (auto it(it1 + 1); it != old_points.end(); ++it)
      linear_ring->addPoint(&*it);
  }
  _geometry2->removeRing(0);
  _geometry2->addRingDirectly(linear_ring);
  spdlog::debug("Grafting the seamline to the geometry2 - done");
}

void InitGdal(const std::string& app_path) {
  _putenv_s("PROJ_LIB", fs::path(app_path).parent_path().string().c_str());
  GDALAllRegister();
  CPLSetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS");
  CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
  GDALSetCacheMax64(8LL * 1024 * 1024 * 1024);
}

void InitSpdlog(const std::string& name, spdlog::level::level_enum level) {
  auto console_sink(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
  console_sink->set_level(spdlog::level::info);
  auto file_sink(std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "logs/" + name + "-" + GetDate() + ".txt", true));
  file_sink->set_level(spdlog::level::debug);
  spdlog::set_default_logger(std::make_shared<spdlog::logger>(
      name, spdlog::sinks_init_list({console_sink, file_sink})));
  spdlog::default_logger()->flush_on(spdlog::level::debug);
  spdlog::set_level(level);
}

std::map<std::string, std::string> rpb_names{
    {"errBias", "ERR_BIAS"}, {"errRand", "ERR_RAND"},
    {"lineOffset", "LINE_OFF"}, {"sampOffset", "SAMP_OFF"},
    {"latOffset", "LAT_OFF"}, {"longOffset", "LONG_OFF"},
    {"heightOffset", "HEIGHT_OFF"}, {"lineScale", "LINE_SCALE"},
    {"sampScale", "SAMP_SCALE"}, {"latScale", "LAT_SCALE"}, 
    {"longScale", "LONG_SCALE"}, {"heightScale", "HEIGHT_SCALE"},
    {"lineNumCoef", "LINE_NUM_COEFF"}, {"lineDenCoef", "LINE_DEN_COEFF"},
    {"sampNumCoef", "SAMP_NUM_COEFF"}, {"sampDenCoef", "SAMP_DEN_COEFF"}};

char** LoadRPBFile(const std::string& rpb_path) {
  std::ifstream rpb(rpb_path);
  std::string line, word, temp, name, value;
  char** rpc_info(nullptr);
  int idx(0);
  while (std::getline(rpb, line)) {
    std::istringstream ss1(line);
    ss1 >> word;
    if (rpb_names.find(word) != rpb_names.end()) {
      name = rpb_names.at(word);
      if (name.find("COEFF") == std::string::npos) {
        ss1 >> value;
        ss1 >> value;
      } else {
        value.clear();
        for (int i(0); i < 20; ++i) {
          std::getline(rpb, line);
          std::istringstream ss2(line);
          ss2 >> temp;
          value.append(temp);
          value.back() = ' ';
        }
        value.pop_back();
      }
      value.pop_back();
      rpc_info = CSLSetNameValue(rpc_info, name.c_str(), value.c_str());
    }
  }
  return rpc_info;
}

cv::Mat MapMatWithLut(const cv::Mat& source_mat, const cv::Mat& lut_mat) {
  if (lut_mat.depth() != CV_8U && lut_mat.depth() != CV_16U) return cv::Mat();

  cv::Mat output_mat(source_mat.size(), source_mat.type());
  if (lut_mat.depth() == CV_8U) {
    cv::LUT(source_mat, lut_mat, output_mat);
  } else {
    auto bands_count(source_mat.channels());
    auto lut_ptr(lut_mat.ptr<uint16_t>(0));
#pragma omp parallel for schedule(dynamic)
    for (int row(0); row < source_mat.rows; ++row) {
      auto source_ptr(source_mat.ptr<uint16_t>(row));
      auto output_ptr(output_mat.ptr<uint16_t>(row));
      for (int col(0); col < source_mat.cols; ++col) {
        for (int b(0); b < bands_count; ++b, ++source_ptr, ++output_ptr)
          *output_ptr = *(lut_ptr + bands_count * *source_ptr + b);
      }
    }
  }
  return output_mat;
}

void RearrangeNeighborGeometries(
    OGRGeometry* seamline,
    OGRGeometryUniquePtr& geometry1,
    OGRGeometryUniquePtr& geometry2) {
  std::vector<OGRPoint> seamline_points;
  for (const auto& point : seamline->toLineString())
    seamline_points.push_back(point);
  std::vector<OGRPoint> old_points;
  std::vector<OGRPoint>::iterator it1, it2;
  OGRLinearRing* linear_ring(nullptr);

  // Rearrange points in the geometry1
  if (auto _geometry1(geometry1->toPolygon()); 
      std::find(
          seamline_points.begin(), seamline_points.end(),
          *(_geometry1->getExteriorRing()->begin())) != seamline_points.end()) {
    for (const auto& point : _geometry1->getExteriorRing())
      old_points.push_back(point);
    int trunc1(-1), trunc2(-1);
    do {
      it1 = std::find(
          old_points.begin(), old_points.end(), 
          *(seamline_points.begin() + (++trunc1)));
    } while (it1 == old_points.end());
    do {
      it2 = std::find(
          old_points.begin(), old_points.end(),
          *(seamline_points.rbegin() + (++trunc2)));
    } while (it2 == old_points.end());
    if (auto it(std::find(it1 + 1, old_points.end(), *it1));
        it != old_points.end() &&
        (it - it2 + 1 + trunc1 + trunc2) != seamline_points.size()) {
      it1 = it;
    } else if (auto it(std::find(it2 + 1, old_points.end(), *it2));
        it != old_points.end() &&
        (it - it1 + 1 + trunc1 + trunc2) != seamline_points.size()) {
      it2 = it;
    }
    linear_ring = new OGRLinearRing;
    if (it2 > it1) {
      for (auto it(it1 + (it2 - it1) / 2); it != it2; ++it)
        linear_ring->addPoint(&*it);
      for (auto it(seamline_points.rbegin()); it != seamline_points.rend();
          ++it)
        linear_ring->addPoint(&*it);
      for (auto it(it1 + 1); it != it1 + (it2 - it1) / 2 + 1; ++it)
        linear_ring->addPoint(&*it);
    } else {
      for (auto it(it2 + (it1 - it2) / 2); it != it1; ++it)
        linear_ring->addPoint(&*it);
      for (auto it(seamline_points.begin()); it != seamline_points.end(); ++it)
        linear_ring->addPoint(&*it);
      for (auto it(it2 + 1); it != it2 + (it1 - it2) / 2 + 1; ++it)
        linear_ring->addPoint(&*it);
    }
    _geometry1->removeRing(0);
    _geometry1->addRingDirectly(linear_ring);
  }

  // Rearrange points in the geometry2
  if (auto _geometry2(geometry2->toPolygon());
      std::find(
          seamline_points.begin(), seamline_points.end(),
          *(_geometry2->getExteriorRing()->begin())) != seamline_points.end()) {
    old_points.resize(0);
    for (const auto& point : _geometry2->getExteriorRing())
      old_points.push_back(point);
    int trunc1(-1), trunc2(-1);
    do {
      it1 = std::find(
          old_points.begin(), old_points.end(),
          *(seamline_points.begin() + (++trunc1)));
    } while (it1 == old_points.end());
    do {
      it2 = std::find(
          old_points.begin(), old_points.end(),
          *(seamline_points.rbegin() + (++trunc2)));
    } while (it2 == old_points.end());
    if (auto it(std::find(it1 + 1, old_points.end(), *it1));
        it != old_points.end() &&
        (it - it2 + 1 + trunc1 + trunc2) != seamline_points.size()) {
      it1 = it;
    } else if (auto it(std::find(it2 + 1, old_points.end(), *it2));
        it != old_points.end() &&
        (it - it1 + 1 + trunc1 + trunc2) != seamline_points.size()) {
      it2 = it;
    }
    linear_ring = new OGRLinearRing;
    if (it2 > it1) {
      for (auto it(it1 + (it2 - it1) / 2); it != it2; ++it)
        linear_ring->addPoint(&*it);
      for (auto it(seamline_points.rbegin()); it != seamline_points.rend();
          ++it)
        linear_ring->addPoint(&*it);
      for (auto it(it1 + 1); it != it1 + (it2 - it1) / 2 + 1; ++it)
        linear_ring->addPoint(&*it);
    } else {
      for (auto it(it2 + (it1 - it2) / 2); it != it1; ++it)
        linear_ring->addPoint(&*it);
      for (auto it(seamline_points.begin()); it != seamline_points.end(); ++it)
        linear_ring->addPoint(&*it);
      for (auto it(it2 + 1); it != it2 + (it1 - it2) / 2 + 1; ++it)
        linear_ring->addPoint(&*it);
    }
    _geometry2->removeRing(0);
    _geometry2->addRingDirectly(linear_ring);
  }
}

bool WarpByGeometry(
    const std::vector<GDALDataset*>& source_datasets,
    const std::vector<OGRGeometry*>& geometries,
    GDALDatasetUniquePtr& output_dataset,
    const std::vector<int>& bands_map,
    GDALResampleAlg resample_arg,
    double blend_dist,
    double nodata_value) {
  if (source_datasets.size() != geometries.size()) return false;

  auto bands_count(
      bands_map.empty() ? output_dataset->GetRasterCount()
      : static_cast<int>(bands_map.size()));
  auto actual_bands_map(std::make_unique<int[]>(bands_count));
  if (bands_map.empty()) {
    for (int b(0); b < bands_count; ++b)
      actual_bands_map[b] = b + 1;
  } else {
    for (int b(0); b < bands_count; ++b)
      actual_bands_map[b] = bands_map[b];
  }
  auto nodata_values(std::make_unique<double[]>(bands_count));
  for (int b(0); b < bands_count; ++b)
    nodata_values[b] = nodata_value;
  std::unique_ptr<void, void(*)(void*)> trans_arg(
      nullptr, [](void* p) { GDALDestroyGenImgProjTransformer(p); });
  std::unique_ptr<GDALWarpOptions> warp_options(GDALCreateWarpOptions());
  warp_options->dfCutlineBlendDist = blend_dist;
  warp_options->eResampleAlg = resample_arg;
  warp_options->eWorkingDataType = 
      output_dataset->GetRasterBand(1)->GetRasterDataType();
  warp_options->hDstDS = output_dataset.get();
  warp_options->nBandCount = bands_count;
  warp_options->padfDstNoDataReal = nodata_values.get();
  warp_options->padfSrcNoDataReal = nodata_values.get();
  warp_options->panDstBands = actual_bands_map.get();
  warp_options->panSrcBands = actual_bands_map.get();
  warp_options->pfnTransformer = GDALGenImgProjTransform;
  warp_options->dfWarpMemoryLimit = 8.0 * 1024 * 1024 * 1024;
  warp_options->papszWarpOptions = CSLSetNameValue(
      warp_options->papszWarpOptions, "SKIP_NOSOURCE", "YES");
  warp_options->papszWarpOptions = CSLSetNameValue(
      warp_options->papszWarpOptions, "CUTLINE_ALL_TOUCHED", "TRUE");
  GDALWarpOperation warp_operation;
  for (int i(0); i < source_datasets.size(); ++i) {
    trans_arg.reset(GDALCreateGenImgProjTransformer2(
        source_datasets[i], output_dataset.get(), nullptr));
    OGRGeometryUniquePtr cutline(nullptr);
    if (geometries[i]) {
      double geotrans[6], inv_geotrans[6];
      source_datasets[i]->GetGeoTransform(geotrans);
      GDALInvGeoTransform(geotrans, inv_geotrans);
      cutline.reset(geometries[i]->clone());
      if (cutline->getSpatialReference()) {
        cutline->transformTo(const_cast<OGRSpatialReference*>(
            source_datasets[i]->GetSpatialRef()));
      }
      cutline = ApplyGeotransOnGeometry(cutline.get(), inv_geotrans);
      warp_options->hCutline = cutline.get();
    }
    warp_options->hSrcDS = source_datasets[i];
    warp_options->pTransformerArg = trans_arg.get();
    warp_operation.Initialize(warp_options.get());
    warp_operation.ChunkAndWarpMulti(
        0, 0, output_dataset->GetRasterXSize(),
        output_dataset->GetRasterYSize());
  }
  return true;
}

}  // namespace utils
}  // namespace rs_toolset