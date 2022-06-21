#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_
#define RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_

#include <memory>
#include <string>
#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/pansharpening.h>


namespace rs_toolset {
namespace pansharpening {

class PansharpeningBase : virtual public PansharpeningInterface {
 public:
  PansharpeningBase(int block_size) : block_size_(block_size) {}
  PansharpeningBase(const PansharpeningBase&) = delete;
  PansharpeningBase& operator=(const PansharpeningBase&) = delete;
  virtual ~PansharpeningBase() = default;

 protected:
  typedef std::unique_ptr<void, void(*)(void*)> RPCTransPtr;
  struct Data {
    cv::Mat mat;
    std::vector<cv::Mat> mats;
  };

  void CreateRangesForRPC(
      int target_block_idx,
      int target_block_cols_count,
      int target_block_rows_count,
      int last_target_block_x_size,
      int last_target_block_y_size,
      int target_block_x_size,
      int target_block_y_size,
      int source_x_size,
      int source_y_size,
      const RPCTransPtr& forward_trans_arg,
      const RPCTransPtr& backward_trans_arg,
      int(&target_range)[4],
      int(&source_range)[4]);

  Data CreateUpsampledDataWithRPC(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range,
      int* ms_range,
      const RPCTransPtr& pan_trans_arg,
      const RPCTransPtr& ms_trans_arg);

  Data CreateDownsampledDataWithRPC(
      const std::string& pan_path,
      const std::string& ms_path,
      int* pan_range,
      int* ms_range,
      const RPCTransPtr& pan_trans_arg,
      const RPCTransPtr& ms_trans_arg);

  Data CreateUpsampledDataWithGeotrans(
      const std::string& pan_path,
      const std::string& ms_path,
      int* target_range);

  Data CreateDownsampledDataWithGeotrans(
      const std::string& pan_path,
      const std::string& ms_path,
      int* target_range);

  int block_size_;

 private:
  std::vector<cv::Mat> CreateResampledMats(
      const cv::Mat& source_mat,
      const cv::Mat& target_mat,
      int* source_range,
      int* target_range,
      const RPCTransPtr& forward_trans_arg,
      const RPCTransPtr& backward_trans_arg);

  void MaskData(Data& data);
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_PANSHARPENING_BASE_H_