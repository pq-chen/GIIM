#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_
#define RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_

#include <string>
#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/pansharpening.h>
#include <rs-toolset/utils.hpp>

namespace rs_toolset {
namespace pansharpening {

class PansharpeningBase : virtual public PansharpeningInterface {
 public:
  explicit PansharpeningBase(int block_size) : block_size_(block_size) {
    spdlog::info(
        "Creating a pansharpening base with\n - Block size: {}", block_size);
  }
  PansharpeningBase(const PansharpeningBase&) = delete;
  PansharpeningBase& operator=(const PansharpeningBase&) = delete;
  virtual ~PansharpeningBase() = default;

 protected:
  struct Data {
    cv::Mat pan_mat;
    std::vector<cv::Mat> ms_mats;
  };

  /**
   * @brief Create the source range and the target range with the given RPCTransPtrs
   * @param[in] source_x_size The source x size
   * @param[in] source_y_size The source y size
   * @param[in] target_block_idx The target block index
   * @param[in] target_block_cols_count The target block columns' count
   * @param[in] target_block_rows_count The target block rows' count
   * @param[in] target_block_x_size The target block x size
   * @param[in] target_block_y_size The target block y size
   * @param[in] last_target_block_x_size The last target block x size
   * @param[in] last_target_block_y_size The last target block y size
   * @param[in] source_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the source RPC raster pixel(row/column)
   * @param[in] target_trans_arg The RPCTransPtr from the target RPC raster pixel(row/column) to the coordinate(longitude/latitude)
   * @param[out] source_range The output source range
   * @param[out] target_range The output target range
  */
  void CreateRangesForRPC(
      int source_x_size,
      int source_y_size,
      int target_block_idx,
      int target_block_cols_count,
      int target_block_rows_count,
      int target_block_x_size,
      int target_block_y_size,
      int last_target_block_x_size,
      int last_target_block_y_size,
      const utils::RPCTransPtr& source_trans_arg,
      const utils::RPCTransPtr& target_trans_arg,
      int(&source_range)[4],
      int(&target_range)[4]);

  /**
   * @brief Create the downsampled data with the RPC information
   * @param[in] pan_dataset The PAN raster dataset
   * @param[in] ms_dataset The MS raster dataset
   * @param[in] pan_range The PAN raster range
   * @param[in] ms_range The MS raster range
   * @param[in] pan_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the PAN raster pixel(row/column)
   * @param[in] ms_trans_arg The RPCTransPtr from the MS raster pixel(row/column) to the coordinate(longitude/latitude)
   * @return The output downsampled data
  */
  Data CreateDownsampledDataWithRPC(
      GDALDataset* pan_dataset,
      GDALDataset* ms_dataset,
      int* pan_range,
      int* ms_range,
      const utils::RPCTransPtr& pan_trans_arg,
      const utils::RPCTransPtr& ms_trans_arg);

  /**
   * @brief Create the downsampled data with the geotrans
   * @param[in] pan_dataset The PAN raster dataset
   * @param[in] ms_dataset The MS raster dataset
   * @param[in] ms_range The MS raster range
   * @return The output downsampled data
  */
  Data CreateDownsampledDataWithGeotrans(
      GDALDataset* pan_dataset,
      GDALDataset* ms_dataset,
      int* ms_range);

  /**
   * @brief Create the upsampled data with the RPC information
   * @param[in] pan_dataset The PAN raster dataset
   * @param[in] ms_dataset The MS raster dataset
   * @param[in] pan_range The PAN raster range
   * @param[in] ms_range The MS raster range
   * @param[in] pan_trans_arg The RPCTransPtr from the PAN raster pixel(row/column) to the coordinate(longitude/latitude)
   * @param[in] ms_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the MS raster pixel(row/column)
   * @return The output upsampled data
  */
  Data CreateUpsampledDataWithRPC(
      GDALDataset* pan_dataset,
      GDALDataset* ms_dataset,
      int* pan_range,
      int* ms_range,
      const utils::RPCTransPtr& pan_trans_arg,
      const utils::RPCTransPtr& ms_trans_arg);

  /**
   * @brief Create the upsampled data with the geotrans
   * @param[in] pan_dataset The PAN raster dataset
   * @param[in] ms_dataset The MS raster dataset
   * @param[in] pan_range The PAN raster range
   * @return The output upsampled data
  */
  Data CreateUpsampledDataWithGeotrans(
      GDALDataset* pan_dataset,
      GDALDataset* ms_dataset,
      int* pan_range);

  int block_size_;

 private:
  /**
   * @brief Create the resampled source mats with the target mat size
   * @param[in] source_mat The source mat
   * @param[in] source_range The source raster range
   * @param[in] target_range The target raster range
   * @param[in] source_trans_arg The RPCTransPtr from the coordinate(longitude/latitude) to the source RPC raster pixel(row/column)
   * @param[in] target_trans_arg The RPCTransPtr from the target RPC raster pixel(row/column) to the coordinate(longitude/latitude)
   * @return The output resampled source mats
  */
  std::vector<cv::Mat> CreateResampledMatsWithRpc(
      const cv::Mat& source_mat,
      int* source_range,
      int* target_range,
      const utils::RPCTransPtr& source_trans_arg,
      const utils::RPCTransPtr& target_trans_arg);

  /**
   * @brief Make all mats in the given data with the same mask range
   * @param[in,out] data The given data
  */
  void MaskData(Data& data);
};

}  // namespace pansharpening
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_