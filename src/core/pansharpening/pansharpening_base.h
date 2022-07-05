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
        "Creating a pansharpening base with\n - Block size: {}", block_size);
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

  /**
   * @brief Create a target range and a source range with the RPC information
   * @param[in] target_block_idx The target block index
   * @param[in] target_block_cols_count The target block columns' count
   * @param[in] target_block_rows_count The target block rows' count
   * @param[in] target_block_x_size The target block x size
   * @param[in] target_block_y_size The target block y size
   * @param[in] last_target_block_x_size The last target block x size
   * @param[in] last_target_block_y_size The last target block y size
   * @param[in] source_x_size The source block x size
   * @param[in] source_y_size The source block y size
   * @param[in] forward_trans_arg The RPCTransPtr from the target RPC raster pixel(row/column) to the coordinate(longitude/latitude)
   * @param[in] backward_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the source RPC raster pixel(row/column)
   * @param[out] target_range The output target range
   * @param[out] source_range The output source range
  */
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

  /**
   * @brief Create the PAN mat and the upsampled MS mat with the RPC information
   * @param[in] pan_path The PAN raster path
   * @param[in] ms_path The MS raster path
   * @param[in] pan_range The PAN raster range
   * @param[in] ms_range The MS raster range
   * @param[in] pan_trans_arg The RPCTransPtr from the PAN raster pixel(row/column) to the coordinate(longitude/latitude)
   * @param[in] ms_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the MS raster pixel(row/column)
   * @return The output PAN mat and the upsampled MS mat
  */
  Data CreateUpsampledDataWithRPC(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range,
      int* ms_range,
      const RPCTransPtr& pan_trans_arg,
      const RPCTransPtr& ms_trans_arg);

  /**
   * @brief Create the downsampled PAN mat and the MS mat with the RPC information
   * @param[in] pan_path The PAN raster path
   * @param[in] ms_path The MS raster path
   * @param[in] pan_range The PAN raster range
   * @param[in] ms_range The MS raster range
   * @param[in] pan_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the PAN raster pixel(row/column)
   * @param[in] ms_trans_arg The RPCTransPtr from the MS raster pixel(row/column) to the coordinate(longitude/latitude)
   * @return The output downsampled PAN mat and the MS mat
  */
  Data CreateDownsampledDataWithRPC(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range,
      int* ms_range,
      const RPCTransPtr& pan_trans_arg,
      const RPCTransPtr& ms_trans_arg);

  /**
   * @brief Create the PAN mat and the upsampled MS mat with the geotransform
   * @param[in] pan_path The PAN raster path
   * @param[in] ms_path The MS raster path
   * @param[in] pan_range The PAN raster range
   * @return The output PAN mat and the upsampled MS mat
  */
  Data CreateUpsampledDataWithGeotrans(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range);

  /**
   * @brief Create the downsampled PAN mat and the MS mat with the geotransform
   * @param[in] pan_path The PAN raster path
   * @param[in] ms_path The MS raster path
   * @param[in] ms_range The MS raster range
   * @return The output downsampled PAN mat and the MS mat
  */
  Data CreateDownsampledDataWithGeotrans(
      const std::string& pan_path,
      const std::string& ms_path,
      int* ms_range);

  int block_size_;

 private:
  /**
   * @brief Create the resampled source mat with the target mat size
   * @param[in] source_mat The source mat
   * @param[in] target_mat The target mat
   * @param[in] source_range The source raster range
   * @param[in] target_range The target raster range
   * @param[in] forward_trans_arg The RPCTransPtr from the target RPC raster pixel(row/column) to the coordinate(longitude/latitude)
   * @param[in] backward_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the source RPC raster pixel(row/column)
   * @return The output resampled source mat
  */
  std::vector<cv::Mat> CreateResampledMatsWithRpc(
      const cv::Mat& source_mat,
      const cv::Mat& target_mat,
      int* source_range,
      int* target_range,
      const RPCTransPtr& forward_trans_arg,
      const RPCTransPtr& backward_trans_arg);

  /**
   * @brief Make all mats in the given data with the same mask range
   * @param[in,out] data The given data
  */
  void MaskData(Data& data);
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_