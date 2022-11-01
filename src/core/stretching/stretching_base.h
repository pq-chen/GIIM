#ifndef RS_TOOLSET_SRC_CORE_STRETCHING_STRETCHING_BASE_H_
#define RS_TOOLSET_SRC_CORE_STRETCHING_STRETCHING_BASE_H_

#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretching.h>

namespace rs_toolset {
namespace stretching {

class StretchingBase : virtual public StretchingInterface {
 public:
  bool Run(std::vector<cv::Mat>& mats) override; 

  bool Run(cv::Mat& mat) override;

  void Clear() override;

 protected:
  int depth_ = -1;
  int bands_count_ = -1;
  std::vector<cv::Mat> lut_mats_;
};

}  // namespace stretching
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_STRETCHING_STRETCHING_BASE_H_