#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_UTILS_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_UTILS_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else  // RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif  // RS_TOOLSET_EXPORTS

#else  // _WIN32
#define RS_TOOLSET_API
#endif  // _WIN32

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gdalwarper.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

namespace rs_toolset {
namespace utils {

/**
 * @brief Transform the given geometry's coordinates between the image space and the object space
 * @param[in] geometry The given geometry, available types include OGRPolygon, OGRMultiPolygon, OGRLineString and OGRMultiLineString
 * @param[in] geotrans The geotransform
 * @return The output geometry
*/
RS_TOOLSET_API OGRGeometryUniquePtr ApplyGeotransOnGeometry(
    OGRGeometry* geometry,
    double* geotrans);

template <typename T1, typename T2>
inline void CalcInterpolationsImpl(
    const cv::Mat& mat,
    int points_count,
    double* xs,
    double* ys,
    T1* interpolations,
    bool round_result,
    const T2& nodata_value) {
  int bands_count(mat.channels());
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < points_count; i++) {
    // Calculate the left-top point position and check whether is out of range or not
    int x(static_cast<int>(floor(xs[i]))), y(static_cast<int>(floor(ys[i])));
    if (x < -1 || x >= mat.cols || y < -1 || y >= mat.rows) {
      for (int b = 0; b < bands_count; b++)
        interpolations[bands_count * i + b] = static_cast<T1>(nodata_value);
      continue;
    }

    // Calculate four neighbor point values and weights
    const T2* values[4]{
        (x != -1 && y != -1 && *mat.ptr<T2>(y, x) != nodata_value)
            ? mat.ptr<T2>(y, x) : nullptr,
        (x != mat.cols - 1 && y != -1 &&
        *mat.ptr<T2>(y, x + 1) != nodata_value)
            ? mat.ptr<T2>(y, x + 1) : nullptr,
        (x != -1 && y != mat.rows - 1 &&
        *mat.ptr<T2>(y + 1, x) != nodata_value)
            ? mat.ptr<T2>(y + 1, x) : nullptr,
        (x != mat.cols - 1 && y != mat.rows - 1 &&
        *mat.ptr<T2>(y + 1, x + 1) != nodata_value)
            ? mat.ptr<T2>(y + 1, x + 1) : nullptr};
    double weights[4]{
        values[0] ? (x + 1 - xs[i]) * (y + 1 - ys[i]) : 0.0,
        values[1] ? (xs[i] - x) * (y + 1 - ys[i]) : 0.0,
        values[2] ? (x + 1 - xs[i]) * (ys[i] - y) : 0.0,
        values[3] ? (xs[i] - x) * (ys[i] - y) : 0.0};

    // Traverse all bands and calculate the interpolations
    for (int b = 0; b < bands_count; b++) {
      double sum_value(0.0), sum_weight(0.0);
      for (int j = 0; j < 4; j++) {
        if (values[j]) {
          sum_value += weights[j] * *(values[j] + b);
          sum_weight += weights[j];
        }
      }
      if (round_result) {
        interpolations[bands_count * i + b] = static_cast<T1>(
            sum_weight ? round(sum_value / sum_weight) : nodata_value);
      } else {
        interpolations[bands_count * i + b] = static_cast<T1>(
            sum_weight ? sum_value / sum_weight : nodata_value);
      }
    }
  }
}

/**
 * @brief Calculate interpolations from the given mat
 * @tparam T The interpolation type
 * @param[in] mat The given mat
 * @param[in] points_count The points' count
 * @param[in] xs The x coordinate array
 * @param[in] ys The y coordinate array
 * @param[out] interpolations The output interpolation array
 * @param[in] round_result whether rounds the result or not, default is false
 * @param[in] nodata_value The output nodata value, default is 0.0
*/
template <typename T>
inline void CalcInterpolations(
    const cv::Mat& mat,
    int points_count,
    double* xs,
    double* ys,
    T* interpolations,
    bool round_result = false,
    double nodata_value = 0.0) {
  switch(mat.depth()) {
    case CV_8U: {
      CalcInterpolationsImpl<T, uint8_t>(
          mat, points_count, xs, ys, interpolations, round_result,
          static_cast<uint8_t>(nodata_value));
      break;
    }
    case CV_16U: {
      CalcInterpolationsImpl<T, uint16_t>(
          mat, points_count, xs, ys, interpolations, round_result,
          static_cast<uint16_t>(nodata_value));
      break;
    }
    case CV_16S: {
      CalcInterpolationsImpl<T, int16_t>(
          mat, points_count, xs, ys, interpolations, round_result,
          static_cast<int16_t>(nodata_value));
      break;
    }
    case CV_32S: {
      CalcInterpolationsImpl<T, int32_t>(
          mat, points_count, xs, ys, interpolations, round_result,
          static_cast<int32_t>(nodata_value));
      break;
    }
    case CV_32F: {
      CalcInterpolationsImpl<T, float>(
          mat, points_count, xs, ys, interpolations, round_result,
          static_cast<float>(nodata_value));
      break;
    }
    case CV_64F: {
      CalcInterpolationsImpl<T, double>(
          mat, points_count, xs, ys, interpolations, round_result,
          static_cast<double>(nodata_value));
      break;
    }
    default: {
      return;
    }
  }
}

/**
 * @brief Calculate interpolations from the given mat
 * @tparam T The interpolation type
 * @param[in] mat The given mat
 * @param[in] points_count The points' count
 * @param[in] xs The x coordinate array
 * @param[in] ys The y coordinate array
 * @param[out] interpolations The output interpolation array
 * @param[in] round_result whether rounds the result or not, default is false
 * @param[in] nodata_value The output nodata value, default is 0.0
*/
template <typename T>
inline void CalcInterpolations(
    const cv::Mat& mat,
    double* geotrans,
    int points_count,
    double* xs,
    double* ys,
    T* interpolations,
    bool round_result = false,
    double nodata_value = 0.0) {
  auto _xs(std::make_unique<double[]>(points_count)),
      _ys(std::make_unique<double[]>(points_count));
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < points_count; i++) {
    _xs[i] = (xs[i] - geotrans[0]) / geotrans[1] - 0.5;
    _ys[i] = (ys[i] - geotrans[3]) / geotrans[5] - 0.5;
  }
  CalcInterpolations(
      mat, points_count, _xs.get(), _ys.get(), interpolations, round_result,
      nodata_value);
}

/**
 * @brief Create a border from the given raster dataset
 * @param[in] dataset The given raster dataset
 * @param[in] spatial_ref Transform the border into the given spatial reference
 * @param[in] overview_idx The overview index on the raster dataset, default is 1. Negative means the original size and out of range means the last overview
 * @param[in] buffer The buffer distance in pixels, default is -2.0
 * @param[in] tol The tolerance in pixels for simplifying the border. default is 1.5
 * @return The output border
*/
RS_TOOLSET_API OGRGeometryUniquePtr CreateBorder(
    GDALDataset* dataset,
    OGRSpatialReference* spatial_ref = nullptr,
    int overview_idx = 1,
    double buffer = -2.0,
    double tol = 1.5);

/**
 * @brief Create an unique pointer for the given CPL string
 * @param string_list The given CPL string
 * @return The output unique pointer 
*/
RS_TOOLSET_API std::unique_ptr<char, void(*)(char*)> CreateCplString(
    char* string);

/**
 * @brief Create an unique pointer for the given CSL string list
 * @param string_list The given CSL string list
 * @return The output unique pointer 
*/
RS_TOOLSET_API std::unique_ptr<char*, void(*)(char**)> CreateCslStringList(
    char** string_list);

/**
 * @brief Create a raster dataset from the given mat
 * @param[in] mat The given mat
 * @param[in] path The output raster path
 * @param[in] geotrans The output geotransform, default is nullptr
 * @param[in] spatial_ref The output spatial reference, default is nullptr
 * @param[in] bands_count The output bands' count, default is 0 means all bands
 * @param[in] bands_map The output bands' map, default is nullptr means all bands
 * @return The output raster dataset
*/
RS_TOOLSET_API GDALDatasetUniquePtr CreateDatasetFromMat(
    const cv::Mat& mat,
    const std::string& path,
    double* geotrans = nullptr,
    OGRSpatialReference* spatial_ref = nullptr,
    int bands_count = 0,
    int* bands_map = nullptr);

/**
 * @brief Create a geometry from the given envelope
 * @param[in] enve The given envelope
 * @return The output geometry
*/
RS_TOOLSET_API OGRGeometryUniquePtr CreateGeometryFromEnve(
    const OGREnvelope& enve);

/**
 * @brief Create histogram mats for the given source mat
 * @param[in] source_mat The given source mat
 * @param[in] mask_mat The mask mat, default is empty mat
 * @return The output histogram mats
*/
RS_TOOLSET_API std::vector<cv::Mat> CreateHists(
    const cv::Mat& source_mat,
    const cv::Mat& mask_mat = cv::Mat());

/**
 * @brief Create a histogram matching LUT mat from the the source histogram mats to the target histogram mats
 * @param[in] source_hist_mats The source histogram mats
 * @param[in] target_hist_mats The target histogram mats
 * @return The output histogram matching LUT mat
*/
RS_TOOLSET_API cv::Mat CreateHistMatchingLut(
    const std::vector<cv::Mat>& source_hist_mats,
    const std::vector<cv::Mat>& target_hist_mats);

/**
 * @brief Create a mask raster dataset from the given raster band
 * @param[in] raster_band The given raster band
 * @param[in] overview_idx The overview index on the raster dataset, default is 1. Negative means the original size and out of range means the last overview
 * @return The output mask raster dataset
*/
RS_TOOLSET_API GDALDatasetUniquePtr CreateMaskDataset(
    GDALRasterBand* raster_band,
    int overview_idx = 1);

/**
 * @brief Create a mat from the given raster dataset
 * @param[in] dataset The given raster dataset
 * @param[in] range The output range, default is nullptr means the full range
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
 * @brief Create a mat from the given raster band
 * @param[in] raster_band The given raster bands
 * @param[in] range The output range, default is nullptr means the full range
 * @return The output mat
*/
RS_TOOLSET_API cv::Mat CreateMatFromRasterBand(
    GDALRasterBand* raster_band,
    int* range = nullptr);

/**
 * @brief Joint the given OGRMultiLineString into a OGRLineString
 * @param[in,out] geometry given OGRMultiLineString
 * @return Running state
*/
RS_TOOLSET_API bool JointMultiLineString(OGRGeometryUniquePtr& geometry);

/**
 * @brief Create raster pyramids for the given raster dataset
 * @param[in] dataset The given raster dataset
 * @param[in] compress_method The compress method, available methods include DEFLATE, LZW, PACKBITS and etc, default is DEFLATE
 * @param[in] clean Whether needs to clean the existing pyramids or not, default is false
*/
RS_TOOLSET_API void CreateRasterPyra(
    GDALDataset* dataset,
    const std::string& compress_method = "DEFLATE",
    bool clean = false);

/**
 * @brief Create range for the given block index
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
    int (&range)[4]);

/** @brief RPC source class */
enum class RS_TOOLSET_API RPCSource {
  kInternal = 0,
  kRPBFile = 1, 
  //TODO kRPCFile = 2,
};
using RPCTransPtr = std::unique_ptr<void, void (*)(void*)>;

/**
 * @brief Create a RPC transform unique pointer from the given raster
 * @param[in] path The given raster path
 * @param[in] rpc_source The RPC source, default is kInternal
 * @return The output RPC transform unique pointer
*/
RS_TOOLSET_API RPCTransPtr CreateRPCTrans(
    const std::string& path,
    RPCSource rpc_source = RPCSource::kInternal);

/**
 * @brief Find the date in the given string
 * @param[in] string The given string
 * @return The acquired date
*/
RS_TOOLSET_API int DateMatching(const std::string& string);

/**
 * @brief Extract the biggest polygon in the given geometry
 * @param[in,out] geometry The given geometry which must be OGRMultiPolygon
 * @return Running state
*/
RS_TOOLSET_API bool ExtractBiggestPolygon(OGRGeometryUniquePtr& geometry);

/**
 * @brief Create a string of the current date
 * @return The output string
*/
RS_TOOLSET_API std::string GetDate();

/**
 * @brief Get the corresonding raster driver pointer from the given path
 * @param[in] path The given path
 * @param[in] Compression Check the compression method if not empty
 * @return The corresonding raster driver pointer
*/
RS_TOOLSET_API GDALDriver* GetRasterDriverByPath(
    const std::string& path,
    const std::string& compression = "");

/**
 * @brief Get the corresonding vector driver pointer from the given path
 * @param[in] path The given path
 * @return The corresonding vector driver pointer
*/
RS_TOOLSET_API GDALDriver* GetVectorDriverByPath(const std::string& path);

/**
 * @brief Graft the given seamline to neighbor geometries
 * @param[in,out] seamline The given seamline
 * @param[in,out] geometry1 The geometry1, type is OGRPolygon
 * @param[in,out] geometry2 The geometry2, type is OGRPolygon
 * @param[in] tol The tolerance in pixels for simplifying the seamline, default is 0.0
*/
RS_TOOLSET_API void GraftSeamline(
    OGRGeometryUniquePtr& seamline,
    OGRGeometryUniquePtr& geometry1,
    OGRGeometryUniquePtr& geometry2,
    double tol = 0.0);

/**
 * @brief Initialize the GDAL configuration
 * @param[in] app_path The application path, which always acquired by argv[0]
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
 * @brief Load RPC informaiton from the given RPB file into a string list unique pointer
 * @param[in] rpb_path The given RPB file path
 * @return The output string list unique pointer
*/
RS_TOOLSET_API std::unique_ptr<char*, void(*)(char**)> LoadRPBFile(
    const std::string& rpb_path);

/**
 * @brief Map the given source mat with the LUT mat
 * @param[in] source_mat The given source mat
 * @param[in] lut_mat The LUT mat
 * @return The output mat
*/
RS_TOOLSET_API cv::Mat MapMatWithLut(
    const cv::Mat& source_mat,
    const cv::Mat& lut_mat);

/**
 * @brief Rearrange points in neighbor geometries to expose the given seamline
 * @param[in] seamline The given seamline
 * @param[in,out] geometry1 The geometry1
 * @param[in,out] geometry2 The geometry2
*/
RS_TOOLSET_API void RearrangeNeighborGeometries(
    OGRGeometry* seamline,
    OGRGeometryUniquePtr& geometry1,
    OGRGeometryUniquePtr& geometry2);

/**
 * @brief Warp source datasets to the output dataset using the given geometries as cutlines. This method does not use overviews
 * @param[in] source_datasets The source datasets
 * @param[in] geometries The given geometries
 * @param[in,out] output_dataset The output dataset
 * @param[in] bands_map The bands' map, default is empty means all bands
 * @param[in] resample_arg The output resample argument, default is GRA_Bilinear
 * @param[in] blend_dist The output blend distance, default is 0.0
 * @param[in] nodata_value The output nodata value, default is 0.0
 * @return Running state
*/
RS_TOOLSET_API bool WarpByGeometry(
    const std::vector<GDALDataset*>& source_datasets,
    const std::vector<OGRGeometry*>& geometries,
    GDALDatasetUniquePtr& output_dataset,
    const std::vector<int>& bands_map = {},
    GDALResampleAlg resample_arg = GRA_Bilinear,
    double blend_dist = 0.0,
    double nodata_value = 0.0);

/**
 * @brief Warp source datasets to the output dataset using the given geometries as cutlines. This method uses overviews as possible
 * @param[in] source_paths The source paths
 * @param[in] geometries The given geometries
 * @param[in,out] output_dataset The output dataset
 * @param[in] bands_map The bands' map, default is empty means all bands
 * @param[in] resample_arg The output resample argument, default is GRA_Bilinear
 * @param[in] blend_dist The output blend distance, default is 0.0
 * @param[in] nodata_value The output nodata value, default is 0.0
 * @return Running state
*/
RS_TOOLSET_API bool WarpByGeometry(
    const std::vector<std::string>& source_paths,
    const std::vector<OGRGeometry*>& geometries,
    GDALDatasetUniquePtr& output_dataset,
    const std::vector<int>& bands_map = {},
    GDALResampleAlg resample_arg = GRA_Bilinear,
    double blend_dist = 0.0,
    double nodata_value = 0.0);

}  // namespace utils
}  // namespace rs_toolset

#endif  // RS_TOOLSET_INCLUDE_RS_TOOLSET_UTILS_H_