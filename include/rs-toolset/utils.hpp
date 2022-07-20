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

/**
 * @brief Transform the given geometry's coordinates between the image space and the object space
 * @param[in] geometry The given geometry. Available types include OGRPolygon, OGRMultiPolygon, OGRLineString and OGRMultiLineString
 * @param[in] geotrans The geotransform
 * @return The output transformed geometry
*/
RS_TOOLSET_API OGRGeometryUniquePtr ApplyGeoTransToPolygon(
    OGRGeometry* geometry,
    double* geotrans);

/**
 * @brief Create histogram mats of the given mat
 * @param[in] mat The given mat
 * @param[in] mask_mat The mask mat, default is empty
 * @return The output histogram mats
*/
RS_TOOLSET_API std::vector<cv::Mat> CreateHist(
    const cv::Mat& mat,
    const cv::Mat& mask_mat = cv::Mat());

/**
 * @brief Create the border geometry from the given source raster
 * @param[in] source_raster_dataset The given source raster dataset
 * @param[in] downsample_factor The downsample factor on the source raster, default is 1
 * @return The output border geometry
*/
RS_TOOLSET_API OGRGeometryUniquePtr CreateBorderGeometry(
    GDALDataset* source_raster_dataset,
    int downsample_factor = 1);

/**
 * @brief Create the raster dataset from the given mat
 * @param[in] mat The given mat
 * @param[in] path The raster path
 * @param[in] geotrans The output geotransform, default is nullptr
 * @param[in] spatial_ref The output spatial reference, default is nullptr
 * @param[in] bands_count The output bands' count, default is 0 means all bands
 * @param[in] bands_map The output bands' map, default is nullptr means all bands
 * @param[in] driver_name The driver nam of the output raster dataset, default is MEM
 * @return The output raster dataset
*/
RS_TOOLSET_API GDALDatasetUniquePtr CreateDatasetFromMat(
    const cv::Mat& mat,
    const std::string& path,
    double* geotrans = nullptr,
    OGRSpatialReference* spatial_ref = nullptr,
    int bands_count = 0,
    int* bands_map = nullptr,
    const std::string& driver_name = "MEM");

/**
 * @brief Create the histogram matching LUT mat from the the source histogram to the target histogram
 * @param[in] source_hist_mats The source histogram mats
 * @param[in] target_hist_mats The Target histogram mats
 * @return The output histogram matching LUT mat
*/
RS_TOOLSET_API cv::Mat CreateHistMatchingLut(
    const std::vector<cv::Mat>& source_hist_mat,
    const std::vector<cv::Mat>& target_hist_mat);

/**
 * @brief Create the mask raster dataset from the given source raster band
 * @param[in] source_raster_band The given source raster band
 * @param[in] downsample_factor The downsample factor on the source raster band, default is 1
 * @return The output mask raster dataset
*/
RS_TOOLSET_API GDALDatasetUniquePtr CreateMaskRasterDataset(
    GDALRasterBand* source_raster_band,
    int downsample_factor = 1);

/**
 * @brief Create a mat from the given raster
 * @param[in] dataset The given raster dataset
 * @param[in] range The output raster range, default is nullptr means the full range
 * @param[in] bands_count The output bands' count, default is 0 means all bands
 * @param[in] bands_map The output bands' map, default is nullptr means all bands
 * @return The output mat
*/
RS_TOOLSET_API cv::Mat CreateMatFromDataset(
    GDALDataset* dataset,
    int* range = nullptr,
    int bands_count = 0,
    int* bands_map = nullptr);

/**
 * @brief Creat raster pyramids for the given raster
 * @param[in] dataset The given raster dataset
 * @param[in] compress_method The output compress method. Available methods include DEFLATE, LZW, PACKBITS and etc, default is DEFLATE
 * @param[in] clean Whether needs to clean the existed overviews, default is false
*/
RS_TOOLSET_API void CreateRasterPyra(
    GDALDataset* dataset,
    const std::string& compress_method = "DEFLATE",
    bool clean = false);

/**
 * @brief Create a range for the given block index
 * @param[in] block_idx The given block index
 * @param[in] block_cols_count The block columns' count
 * @param[in] block_rows_count The block rows' count
 * @param[in] block_x_size The block x size
 * @param[in] block_y_size The block y size
 * @param[in] last_block_x_size The last block x size
 * @param[in] last_block_y_size The last block y size
 * @param[out] range The output range
*/
RS_TOOLSET_API void CreateRange(
    int block_idx,
    int block_cols_count,
    int block_rows_count,
    int block_x_size,
    int block_y_size,
    int last_block_x_size,
    int last_block_y_size,
    int(&range)[4]);

/** @brief RPC source class*/
enum class RS_TOOLSET_API RPCSource {
  kInternal = 0,
  kRPBFile = 1, 
  //TODO kRPCFile = 2,
};
typedef std::unique_ptr<void, void(*)(void*)> RPCTransPtr;

/**
 * @brief Create a RPC transform unique pointer from the given raster
 * @param[in] path The given raster path
 * @param[in] reserved Whether the RPC transform is reserved
 * @param[in] rpc_source The RPC source, default is kInternal
 * @return The output RPC transform unique pointer
*/
RS_TOOLSET_API RPCTransPtr CreateRPCTrans(
    const std::string& path,
    bool reserved,
    RPCSource rpc_source = RPCSource::kInternal);

/**
 * @brief Find the date number in the given string
 * @param[in] string The given string
 * @return The output date number
*/
RS_TOOLSET_API int DateMatching(const std::string& string);

/**
 * @brief Find the biggest polygon from the given geometry
 * @param[in,out] geometry The given geometry, available types is OGRMultiPolygon
*/
RS_TOOLSET_API void FindBiggestPolygon(OGRGeometryUniquePtr& geometry);

/**
 * @brief Transform the current date to the date string
 * @return The output date string
*/
RS_TOOLSET_API std::string GetDate();

/**
 * @brief Initialize the GDAL configuration
 * @param[in] app_path The application path
*/
RS_TOOLSET_API void InitGdal(const std::string& app_path);

/**
 * @brief Initialize the Spdlog configuration
 * @param[in] name The logger name
 * @param[in] level The logger level, default is spdlog::level::info
*/
RS_TOOLSET_API void InitSpdlog(
    const std::string& name,
    spdlog::level::level_enum level = spdlog::level::info);

/**
 * @brief Load the RPB file into the RPC string list
 * @param[in] rpb_path The RPB file path
 * @return The output RPC string list
*/
RS_TOOLSET_API char** LoadRPBFile(const std::string& rpb_path);

/**
 * @brief Transform the source mat with the LUT mat
 * @param[in] source_mat The source mat
 * @param[in] lut_mat The LUT mat
 * @return The output transformed mat
*/
RS_TOOLSET_API cv::Mat TransformMatWithLut(
    const cv::Mat& source_mat,
    const cv::Mat& lut_mat);

/**
 * @brief Warp source datasets to the output dataset by the given geometries as cutlines
 * @param[in] source_datasets Source datasets
 * @param[in] geometries The given geometries
 * @param[out] output_dataset The output dataset
 * @param[in] resample_arg The resample argument, default is GRA_Bilinear
 * @param[in] blend_dist The blend distance, default is 0.0
 * @param[in] nodata_value The nodata value, default is 0.0
 * @return Running state
*/
RS_TOOLSET_API bool WarpByGeometry(
    const std::vector<GDALDataset*>& source_datasets,
    const std::vector<OGRGeometry*>& geometries,
    GDALDatasetUniquePtr& output_dataset,
    GDALResampleAlg resample_arg = GRA_Bilinear,
    double blend_dist = 0.0,
    double nodata_value = 0.0);

/**
 * @brief Calculate interpolations from the given mat
 * @tparam T The interpolation type
 * @param mat The given mat
 * @param points_count The points' count
 * @param xs The x coordinate array
 * @param ys The y coordinate array
 * @param interpolations The output interpolation array
 * @param nodata_value The nodata value, default is 0.0
*/
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
      double sum_value(0.0), sum_weight(0.0);
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
    GDALRasterBand* source_raster_band,
    int downsample_factor) {
  auto dataset_type(source_raster_band->GetRasterDataType());
  int x_size(source_raster_band->GetXSize()),
      y_size(source_raster_band->GetYSize()),
      new_x_size(source_raster_band->GetXSize() / downsample_factor),
      new_y_size(source_raster_band->GetYSize() / downsample_factor),
      bytes_count(GDALGetDataTypeSizeBytes(dataset_type));
  double nodata_value(source_raster_band->GetNoDataValue());
  uint64_t size(static_cast<uint64_t>(new_x_size) * new_y_size);
  auto data(std::make_unique<T[]>(size));
  source_raster_band->RasterIO(
      GF_Read, 0, 0, x_size, y_size, data.get(), new_x_size, new_y_size,
      dataset_type, bytes_count, bytes_count * new_x_size);
  auto mask_data(std::make_unique<uint8_t[]>(size));
  memset(mask_data.get(), 0, size * sizeof(uint8_t));
#pragma omp parallel for schedule(dynamic)
  for (int row = 0; row < new_y_size; row++) {
    T* source_ptr(data.get() + row * new_x_size);
    uint8_t* mask_ptr(mask_data.get() + row * new_x_size);
    for (int col = 0; col < new_x_size; col++) {
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