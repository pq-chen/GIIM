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
/// Transform the given geometry's coordinates between the image space and the object space
/// </summary>
/// <param name="geometry">The given geometry. Available types include OGRPolygon, OGRMultiPolygon, OGRLineString and OGRMultiLineString</param>
/// <param name="geotrans">The geotransform</param>
/// <returns>The transformed geometry</returns>
RS_TOOLSET_API OGRGeometryUniquePtr ApplyGeoTransToPolygon(
    OGRGeometry* geometry,
    double* geotrans);

/// <summary>
/// Calculate histogram mats of the given mat
/// </summary>
/// <param name="mat">The given mat</param>
/// <param name="mask_mat">The mask mat, default is empty mat</param>
/// <returns>Histogram mats</returns>
RS_TOOLSET_API std::vector<cv::Mat> CalcHist(
    const cv::Mat& mat,
    const cv::Mat& mask_mat = cv::Mat());

/// <summary>
/// Create the border geometry for the given source raster dataset
/// </summary>
/// <param name="source_raster_dataset">The given source raster dataset</param>
/// <param name="overview_idx">The overview index operated on, default is 0</param>
/// <returns>The border geometry</returns>
RS_TOOLSET_API OGRGeometryUniquePtr CreateBorderGeometry(
    GDALDataset* source_raster_dataset,
    int overview_idx = 0);

/// <summary>
/// Creaet the raster dataset from the given mat
/// </summary>
/// <param name="mat">The given mat</param>
/// <param name="path">The raster dataset path</param>
/// <param name="geotrans">The geotransform, default is nullptr</param>
/// <param name="spatial_ref">The spatial reference, default is nullptr</param>
/// <param name="bands_count">The bands' count, default is 0 means all bands</param>
/// <param name="bands_map">The bands' map, default is nullptr means all bands</param>
/// <param name="driver_name">The raster dataset driver name, default is MEM</param>
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
/// Create the mask raster dataset for the given source raster band
/// </summary>
/// <param name="source_raster_band">The given source raster band</param>
/// <param name="overview_idx">The overview index operated on, default is 0</param>
/// <returns>The output mask raster dataset</returns>
RS_TOOLSET_API GDALDatasetUniquePtr CreateMaskRasterDataset(
    GDALRasterBand* source_raster_band,
    int overview_idx = 0);

/// <summary>
/// Creaet the mat from the given raster dataset
/// </summary>
/// <param name="dataset">The given raster dataset</param>
/// <param name="range">The range, default is nullptr means the full range</param>
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
/// <param name="compress_method">The compress method. Available methods include DEFLATE, LZW, PACKBITS and etc, default is DEFLATE</param>
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
/// <param name="block_cols_count">The block columns' count</param>
/// <param name="block_rows_count">The block rows' count</param>
/// <param name="block_x_size">The block x size</param>
/// <param name="block_y_size">The block y size</param>
/// <param name="last_block_x_size">The last block x size</param>
/// <param name="last_block_y_size">The last block y size</param>
/// <param name="range">The output range</param>
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
  //TODO: kRPCFile = 2,
};
typedef std::unique_ptr<void, void(*)(void*)> RPCTransPtr;

/// <summary>
/// Create the RPC transform for the given raster path
/// </summary>
/// <param name="path">The given raster path</param>
/// <param name="reserved">Whether the RPC transform is reserved</param>
/// <param name="rpc_source">The RPC source. Available RPC sources include kInternal and kRPBFile, default is kInternal</param>
/// <returns>The output RPC transform</returns>
RS_TOOLSET_API RPCTransPtr CreateRPCTrans(
    const std::string& path,
    bool reserved,
    RPCSource rpc_source = RPCSource::kInternal);

/// <summary>
/// Find the date number in the given string
/// </summary>
/// <param name="string">The given string</param>
/// <returns>The output date number</returns>
RS_TOOLSET_API int DateMatching(const std::string& string);

/// <summary>
/// Find the biggest polygon from the given multipolygon
/// </summary>
/// <param name="geometry">The given multipolygon</param>
/// <returns>The output biggest polygon</returns>
RS_TOOLSET_API OGRGeometryUniquePtr FindBiggestPolygon(OGRGeometry* geometry);

/// <summary>
/// Transform the current date to the date string
/// </summary>
/// <returns>The output date string</returns>
RS_TOOLSET_API std::string GetDate();

/// <summary>
/// Initialize GDAL configuration and the proj.db file should be under the same directory with the application
/// </summary>
/// <param name="app_path">The application path</param>
/// <returns></returns>
RS_TOOLSET_API void InitGdal(const std::string& app_path);

/// <summary>
/// Initialize Spdlog configuration
/// <param name="name">The logger name</param>
/// <param name="level">The logger level, default is spdlog::level::info</param>
/// <returns></returns>
RS_TOOLSET_API void InitSpdlog(
    const std::string& name,
    spdlog::level::level_enum level = spdlog::level::info);

/// <summary>
/// Load RPB file to the RPB string list
/// </summary>
/// <param name="rpb_path">The RPB file path</param>
/// <returns>The output RPB string list</returns>
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
/// Warp source raster datasets to the output raster dataset by the given geometry as a cutline
/// </summary>
/// <param name="source_rasters_dataset">Source raster datasets</param>
/// <param name="output_raster_dataset">The output raster dataset</param>
/// <param name="geometry">The given geometry which can be nullptr</param>
/// <param name="resample_arg">The resample argument, default is GRA_Bilinear</param>
/// <param name="nodata_value">The nodata value, default is 0.0</param>
/// <returns></returns>
RS_TOOLSET_API void WarpByGeometry(
    const std::vector<GDALDataset*>& source_rasters_dataset,
    GDALDataset* output_raster_dataset,
    OGRGeometry* geometry,
    GDALResampleAlg resample_arg = GRA_Bilinear,
    double nodata_value = 0.0);

/// <summary>
/// Calculate interpolations from the given mat
/// </summary>
/// <param name="mat">The given mat</param>
/// <param name="points_count">The points' count</param>
/// <param name="xs">The x coordinate array</param>
/// <param name="ys">The y coordinate array</param>
/// <param name="interpolations">The output interpolation array</param>
/// <param name="nodata_value">The nodata value, default is 0.0</param>
/// <returns></returns>
template <typename T>
inline void CalculateInterpolations(
    const cv::Mat& mat,
    int points_count,
    double* xs,
    double* ys,
    T* interpolations,
    double nodata_value = 0.0) {
  int bands_count = mat.channels();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < points_count; i++) {
    // Calculate the left-top point position and check whether out of range
    int x(static_cast<int>(floor(xs[i]))), 
        y(static_cast<int>(floor(ys[i])));
    if (x < -1 || x >= mat.cols || y < -1 || y >= mat.rows) {
      for (int b = 0; b < bands_count; b++)
        interpolations[bands_count * i + b] = static_cast<T>(nodata_value);
      continue;
    }

    // Calculate four neighbor point positions
    const T* data_pos[4] {
      (x == -1 || y == -1) ? nullptr : mat.ptr<T>(y, x),
      (x == mat.cols - 1 || y == -1) ? nullptr : mat.ptr<T>(y, x + 1),
      (x == -1 || y == mat.rows - 1) ? nullptr : mat.ptr<T>(y + 1, x),
      (x == mat.cols - 1 || y == mat.rows - 1) ? nullptr
          : mat.ptr<T>(y + 1, x + 1) };

    // Traverse all bands and calculate the interpolations
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
          sum_weight != 0.0 ? round(sum_value / sum_weight) : nodata_value);
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
  double nodata_value(source_raster_band->GetNoDataValue());
  uint64_t size(static_cast<uint64_t>(x_size) * y_size);
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