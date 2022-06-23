#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_
#define RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_

#include <memory>
#include <string>
#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/pansharpening.h>


namespace rs_toolset {
namespace pansharpening {

class PansharpeningBase : virtual public PansharpeningInterface {
 public:
  PansharpeningBase(int block_size) : block_size_(block_size) {
    spdlog::info(
        "Initializing the basic pansharpening with\n- Block size: {}",
        block_size);
  }
  PansharpeningBase(const PansharpeningBase&) = delete;
  PansharpeningBase& operator=(const PansharpeningBase&) = delete;
  virtual ~PansharpeningBase() = default;

 protected:
  typedef std::unique_ptr<void, void(*)(void*)> RPCTransPtr;
  struct Data {
    cv::Mat mat;
    std::vector<cv::Mat> mats;
  };

  /// <summary>
  /// Create the target range and the source range with the RPC information
  /// </summary>
  /// <param name="target_block_idx">The target block idx</param>
  /// <param name="target_block_cols_count">The target block cols' count</param>
  /// <param name="target_block_rows_count">The target block rows' count</param>
  /// <param name="target_block_x_size">The target block x size</param>
  /// <param name="target_block_y_size">The target block y size</param>
  /// <param name="last_target_block_x_size">The last target block x size</param>
  /// <param name="last_target_block_y_size">The last target block y size</param>
  /// <param name="source_x_size">The source block x size</param>
  /// <param name="source_y_size">The source block y size</param>
  /// <param name="forward_trans_arg">The RPCTransPtr from the target RPC raster pixel(row/column) to the coordinate(longitude/latitude)</param>
  /// <param name="backward_trans_arg">The RPCTransPtr from the coordinate(longitude/latitude) to the source RPC raster pixel(row/column)</param>
  /// <param name="target_range">The output target range</param>
  /// <param name="source_range">The output source range</param>
  void CreateRangesForRPC(
      int target_block_idx,
      int target_block_cols_count,
      int target_block_rows_count,
      int target_block_x_size,
      int target_block_y_size,
      int last_target_block_x_size,
      int last_target_block_y_size,
      int source_x_size,
      int source_y_size,
      const RPCTransPtr& forward_trans_arg,
      const RPCTransPtr& backward_trans_arg,
      int(&target_range)[4],
      int(&source_range)[4]);

  /// <summary>
  /// Create the PAN mat and the upsampled MS mat with the RPC information
  /// </summary>
  /// <param name="pan_path">The PAN raster path</param>
  /// <param name="ms_path">The MS raster path</param>
  /// <param name="pan_range">The PAN raster range</param>
  /// <param name="ms_range">The MS raster range</param>
  /// <param name="pan_trans_arg">The RPCTransPtr from the PAN raster pixel(row/column) to the coordinate(longitude/latitude)</param>
  /// <param name="ms_trans_arg">The RPCTransPtr from the coordinate(longitude/latitude) to the MS raster pixel(row/column)</param>
  /// <returns>The output PAN mat and the upsampled MS mat</returns>
  Data CreateUpsampledDataWithRPC(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range,
      int* ms_range,
      const RPCTransPtr& pan_trans_arg,
      const RPCTransPtr& ms_trans_arg);

  /// <summary>
  /// Create the downsampled PAN mat and the MS mat with the RPC information
  /// </summary>
  /// <param name="pan_path">The PAN raster path</param>
  /// <param name="ms_path">The MS raster path</param>
  /// <param name="pan_range">The PAN raster range</param>
  /// <param name="ms_range">The MS raster range</param>
  /// <param name="pan_trans_arg">The RPCTransPtr from the coordinate(longitude/latitude) to the PAN raster pixel(row/column)</param>
  /// <param name="ms_trans_arg">The RPCTransPtr from the MS raster pixel(row/column) to the coordinate(longitude/latitude)</param>
  /// <returns>The output downsampled PAN mat and the MS mat</returns>
  Data CreateDownsampledDataWithRPC(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range,
      int* ms_range,
      const RPCTransPtr& pan_trans_arg,
      const RPCTransPtr& ms_trans_arg);

  /// <summary>
  /// Create the PAN mat and the upsampled MS mat with the geotransform
  /// </summary>
  /// <param name="pan_path">The PAN raster path</param>
  /// <param name="ms_path">The MS raster path</param>
  /// <param name="pan_range">The PAN raster range</param>
  /// <returns>The output PAN mat and the upsampled MS mat</returns>
  Data CreateUpsampledDataWithGeotrans(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range);

  /// <summary>
  /// Create the downsampled PAN mat and the MS mat with the geotransform
  /// </summary>
  /// <param name="pan_path">The PAN raster path</param>
  /// <param name="ms_path">The MS raster path</param>
  /// <param name="ms_range">The MS raster range</param>
  /// <returns>The output downsampled PAN mat and the MS mat</returns>
  Data CreateDownsampledDataWithGeotrans(
      const std::string& pan_path,
      const std::string& ms_path,
      int* ms_range);

  int block_size_;

 private:
  /// <summary>
  /// Create the resampled source mat from the source mat with the target mat size in the ranges
  /// </summary>
  /// <param name="source_mat">The source mat</param>
  /// <param name="target_mat">The target mat</param>
  /// <param name="source_range">The source raster range</param>
  /// <param name="target_range">The target raster range</param>
  /// <param name="forward_trans_arg">The RPCTransPtr from the target RPC raster pixel(row/column) to the coordinate(longitude/latitude)</param>
  /// <param name="backward_trans_arg">The RPCTransPtr from the coordinate(longitude/latitude) to the source RPC raster pixel(row/column)</param>
  /// <returns>The output resampled source mat from</returns>
  std::vector<cv::Mat> CreateResampledMatsWithRpc(
      const cv::Mat& source_mat,
      const cv::Mat& target_mat,
      int* source_range,
      int* target_range,
      const RPCTransPtr& forward_trans_arg,
      const RPCTransPtr& backward_trans_arg);

  /// <summary>
  /// Make all mats in the given data with same mask range
  /// </summary>
  /// <param name="data">The given data</param>
  void MaskData(Data& data);
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_