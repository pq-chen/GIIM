#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_UTILS_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_UTILS_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else // LCMAKE_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif // LCMAKE_EXPORTS

#else // _WIN32
#define RS_TOOLSET_API
#endif // _WIN32

#include <cmath>
#include <cstdint>
#include <cstring>

#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <gdalwarper.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>


namespace rs_toolset {
namespace utils {

/// <summary>
/// Transform the given geometry's coordinates between image space and object space
/// </summary>
/// <param name="geometry">The given geometry. Available type includes OGRPolygon, OGRMultiPolygon, OGRLineString and OGRMultiLineString</param>
/// <param name="geotrans">The geotransform</param>
/// <returns>The transformed geometry</returns>
RS_TOOLSET_API OGRGeometryUniquePtr ApplyGeoTransToPolygon(
    OGRGeometry* geometry,
    double* geotrans);

/// <summary>
/// Calculate histogram of the given mat
/// </summary>
/// <param name="mat">The given mat</param>
/// <param name="mask_mat">The mask mat</param>
/// <returns>Histogram mats</returns>
RS_TOOLSET_API std::vector<cv::Mat> CalcHist(
    const cv::Mat& mat,
    const cv::Mat& mask_mat);

/// <summary>
/// Create border geometry for the given source raster dataset
/// </summary>
/// <param name="source_raster_dataset">The given source raster dataset</param>
/// <param name="overview_idx">The overview index operated on, default is 0</param>
/// <returns>The border geometry</returns>
RS_TOOLSET_API OGRGeometryUniquePtr CreateBorderGeometry(
    GDALDataset* source_raster_dataset,
    int overview_idx = 0);

/// <summary>
/// Creaet raster dataset from the given mat
/// </summary>
/// <param name="mat">The given mat</param>
/// <param name="path">The raster dataset path</param>
/// <param name="geotrans">The geotransform, default is nullptr</param>
/// <param name="spatial_ref">The spatial reference, default is nullptr</param>
/// <param name="bands_count">The bands' count, default is 0 means all bands</param>
/// <param name="bands_map">The bands map, default is nullptr means all bands</param>
/// <param name="driver_name">The dataset driver name, default is MEM</param>
/// <returns>The output raster dataset</returns>
RS_TOOLSET_API GDALDatasetUniquePtr CreateDatasetFromMat(
    const cv::Mat& mat,
    const std::string& path,
    double* geotrans = nullptr,
    OGRSpatialReference* spatial_ref = nullptr,
    int bands_count = 0,
    int* bands_map = nullptr,
    const std::string& driver_name = "MEM");

/// <summary>
/// Create the LUT mat from the the source histogram to the target histogram
/// </summary>
/// <param name="source_hist_mats">Source histogram mats</param>
/// <param name="target_hist_mats">Target histogram mats</param>
/// <returns>The output LUT mat</returns>
RS_TOOLSET_API cv::Mat CreateHistLUT(
    const std::vector<cv::Mat>& source_hist_mats,
    const std::vector<cv::Mat>& target_hist_mats);

/// <summary>
/// Create the mask raster dataset for the given source raster dataset
/// </summary>
/// <param name="source_raster_band">The given source raster band</param>
/// <param name="overview_idx">The overview index operated on, default is 0</param>
/// <returns>The mask raster dataset</returns>
RS_TOOLSET_API GDALDatasetUniquePtr CreateMaskRasterDataset(
    GDALRasterBand* source_raster_band,
    int overview_idx = 0);

/// <summary>
/// Creaet mat from the given raster dataset
/// </summary>
/// <param name="dataset">The given raster dataset</param>
/// <param name="range">The range, default is nullptr means full range</param>
/// <param name="bands_count">The bands' count, default is 0 means all bands</param>
/// <param name="bands_map">The bands' map, default is nullptr means all bands</param>
/// <returns>The output mat</returns>
RS_TOOLSET_API cv::Mat CreateMatFromDataset(
    GDALDataset* dataset,
    int* range = nullptr,
    int bands_count = 0,
    int* bands_map = nullptr);

/// <summary>
/// Creat raster pyramids for the given raster dataset
/// </summary>
/// <param name="dataset">The given raster dataset</param>
/// <param name="compress_method">Compress method. Available method includes DEFLATE, LZW, PACKBITS and etc, default is DEFLATE</param>
/// <param name="clean">Whether needs to clean the existed overviews, default is false</param>
/// <returns></returns>
RS_TOOLSET_API void CreateRasterPyra(
    GDALDataset* dataset,
    const std::string& compress_method = "DEFLATE",
    bool clean = false);

/// <summary>
/// Create the range for the given block index
/// </summary>
/// <param name="block_idx">The given block index</param>
/// <param name="block_cols_count">Block columns' count</param>
/// <param name="block_rows_count"> Block rows' count</param>
/// <param name="block_x_size">block x size</param>
/// <param name="block_y_size">block y size</param>
/// <param name="last_block_x_size">Last block x size</param>
/// <param name="last_block_y_size">Last block y size</param>
/// <param name="range">The range</param>
/// <returns></returns>
RS_TOOLSET_API void CreateRange(
    int block_idx,
    int block_cols_count,
    int block_rows_count,
    int block_x_size,
    int block_y_size,
    int last_block_x_size,
    int last_block_y_size,
    int(&range)[4]);

enum class RS_TOOLSET_API RPCSource {
  kInternal = 0,
  kRPBFile = 1, 
  kRPCFile = 2,
};
typedef std::unique_ptr<void, void(*)(void*)> RPCTransPtr;

/// <summary>
/// Create RPC transform for the given raster path
/// </summary>
/// <param name="path">The given raster path</param>
/// <param name="reserved">Whether is reserved</param>
/// <param name="rpc_source">The RPC source, default is RPCSource::kInternal</param>
/// <returns>The output RPC transform pointer</returns>
RS_TOOLSET_API RPCTransPtr CreateRPCTrans(
    const std::string& path,
    bool reserved,
    RPCSource rpc_source = RPCSource::kInternal);

/// <summary>
/// Find date in the given string
/// </summary>
/// <param name="string">The given string</param>
/// <returns>The date number</returns>
RS_TOOLSET_API int DateMatching(const std::string& string);

/// <summary>
/// Find the biggest polygon from the given multipolygon
/// </summary>
/// <param name="geometry">The given multipolygon</param>
/// <returns>The biggest polygon</returns>
RS_TOOLSET_API OGRGeometryUniquePtr FindBiggestPolygon(OGRGeometry* geometry);

/// <summary>
/// Transform the current date into string
/// </summary>
/// <returns>The date string</returns>
RS_TOOLSET_API std::string GetDate();

/// <summary>
/// Initialize GDAL configuration
/// </summary>
/// <param name="app_path">The application path</param>
/// <returns></returns>
RS_TOOLSET_API void InitGdal(const std::string& app_path);

/// <summary>
/// Initialize Spdlog configuration
/// </summary>
/// <param name="logger_name">The logger name</param>
/// <param name="log_level">The logger level, default is spdlog::level::info</param>
/// <returns></returns>
RS_TOOLSET_API void InitSpdlog(
    const std::string& logger_name,
    spdlog::level::level_enum log_level = spdlog::level::info);

/// <summary>
/// Load RPB file
/// </summary>
/// <param name="rpb_path">RPB file path</param>
/// <returns>RPB string list</returns>
RS_TOOLSET_API char** LoadRPBFile(const std::string& rpb_path);

/// <summary>
/// Transform the source mat with the LUT mat
/// </summary>
/// <param name="source_mat">The source mat</param>
/// <param name="lut_mat">The LUT mat</param>
/// <returns>The output transformed mat</returns>
RS_TOOLSET_API cv::Mat TransformMat(
    const cv::Mat& source_mat,
    const cv::Mat& lut_mat);

/// <summary>
/// Warp source raster datasets to the result raster dataset by the given geometry as a cutline
/// </summary>
/// <param name="source_rasters_dataset">The source raster datasets</param>
/// <param name="result_raster_dataset">The result raster dataset</param>
/// <param name="geometry">The given geometry as a cutline and can be nullptr</param>
/// <param name="resample_arg">Resample argument, default is GRA_Bilinear</param>
/// <param name="nodata_value">Nodata value, default is 0.</param>
/// <returns></returns>
RS_TOOLSET_API void WarpByGeometry(
    const std::vector<GDALDataset*>& source_rasters_dataset,
    GDALDataset* result_raster_dataset,
    OGRGeometry* geometry,
    GDALResampleAlg resample_arg = GRA_Bilinear,
    double nodata_value = 0.0);

/// <summary>
/// Get the interpolation from the given mat
/// </summary>
/// <param name="mat">The given mat</param>
/// <param name="points_count">The points' count</param>
/// <param name="xs">An array of x coordinates</param>
/// <param name="ys">An array of y coordinates</param>
/// <param name="interpolations">An array of interpolation results</param>
/// <param name="nodata_value">Nodata value, default is 0.</param>
/// <returns></returns>
template <typename T>
inline void GetInterpolation(
    const cv::Mat& mat,
    int points_count,
    double* xs,
    double* ys,
    T* interpolations,
    double nodata_value = 0.0) {
  int bands_count = mat.channels();
  for (int i = 0; i < points_count; i++) {
    int x(static_cast<int>(floor(xs[i]))), 
        y(static_cast<int>(floor(ys[i])));
    if (x < -1 || x >= mat.cols || y < -1 || y >= mat.rows) {
      for (int b = 0; b < bands_count; b++)
        interpolations[bands_count * i + b] = static_cast<T>(nodata_value);
      continue;
    }
    const T* data_pos[4] {
      (x == -1 || y == -1) ? nullptr : mat.ptr<T>(y, x),
      (x == mat.cols - 1 || y == -1) ? nullptr : mat.ptr<T>(y, x + 1),
      (x == -1 || y == mat.rows - 1) ? nullptr : mat.ptr<T>(y + 1, x),
      (x == mat.cols - 1 || y == mat.rows - 1) ? nullptr
          : mat.ptr<T>(y + 1, x + 1) };
    for (int b = 0; b < bands_count; b++) {
      T values[4]{ 
        data_pos[0] ? *(data_pos[0] + b) : static_cast<T>(nodata_value),
        data_pos[1] ? *(data_pos[1] + b) : static_cast<T>(nodata_value),
        data_pos[2] ? *(data_pos[2] + b) : static_cast<T>(nodata_value),
        data_pos[3] ? *(data_pos[3] + b) : static_cast<T>(nodata_value) };
      double weights[4]{ 
          (x + 1 - xs[i]) * (y + 1 - ys[i]), (xs[i] - x)* (y + 1 - ys[i]),
          (x + 1 - xs[i]) * (ys[i] - y), (xs[i] - x) * (ys[i] - y) };
      double sum_value(0), sum_weight(0);
      for (int j = 0; j < 4; j++) {
        if (values[j] != nodata_value) {
          sum_value += weights[j] * values[j];
          sum_weight += weights[j];
        }
      }
      interpolations[bands_count * i + b] = static_cast<T>(
          sum_weight != 0. ? round(sum_value / sum_weight) : nodata_value);
    }
  }
}

template <typename T>
inline std::unique_ptr<uint8_t[]> CreateMaskRasterDatasetImpl(
    GDALRasterBand* source_raster_band) {
  auto dataset_type(source_raster_band->GetRasterDataType());
  int x_size(source_raster_band->GetXSize()),
      y_size(source_raster_band->GetYSize()),
      bytes_count(GDALGetDataTypeSizeBytes(dataset_type));
  uint64_t size(static_cast<uint64_t>(x_size) * y_size);
  double nodata_value(source_raster_band->GetNoDataValue());
  auto data(std::make_unique<T[]>(size));
  source_raster_band->RasterIO(
      GF_Read, 0, 0, x_size, y_size, data.get(), x_size, y_size, dataset_type,
      bytes_count, bytes_count * x_size);
  auto mask_data(std::make_unique<uint8_t[]>(size));
  memset(mask_data.get(), 0, size * sizeof(uint8_t));
#pragma omp parallel for schedule(dynamic)
  for (int row = 0; row < y_size; row++) {
    T* source_ptr(data.get() + row * x_size);
    uint8_t* mask_ptr(mask_data.get() + row * x_size);
    for (int col = 0; col < x_size; col++) {
      if (*source_ptr != static_cast<T>(nodata_value))
        *mask_ptr = 255;
      source_ptr++;
      mask_ptr++;
    }
  }
  return mask_data;
}

} // utils
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_UTILS_H_