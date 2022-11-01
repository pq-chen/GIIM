#ifndef RS_TOOLSET_SRC_CORE_STRETCHING_PERCENT_CLIP_STRETCHING_H_
#define RS_TOOLSET_SRC_CORE_STRETCHING_PERCENT_CLIP_STRETCHING_H_

#pragma warning(disable:4250)

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretching.h>
#include "stretching_base.h"

namespace rs_toolset {
namespace stretching {

class PercentClipImpl final : public StretchingBase, public PercentClip {
 public:
  explicit PercentClipImpl(double low_percent, double high_percent) 
      : low_percent_(low_percent), high_percent_(high_percent) {
    spdlog::info(
        "Creating a percent clip stretching with\n"
        " - Low percent: {}\n"
        " - High percent: {}", low_percent_, high_percent_);
  }
  PercentClipImpl(const PercentClipImpl&) = delete;
  PercentClipImpl& operator=(const PercentClipImpl&) = delete;
  ~PercentClipImpl() = default;

  bool AccumulateStat(const std::vector<cv::Mat>& mats) override;

  bool AccumulateStat(const cv::Mat& mat) override;
  
  bool CreateLutMats() override;

  void Clear() override;

  void SetPercent(double low_percent, double high_percent) override;

  void GetPercent(double& low_percent, double& high_percent) override;

 private:
  double low_percent_;
  double high_percent_;
  std::vector<cv::Mat> hist_mats_;
};

}  // namespace stretching
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_STRETCHING_PERCENT_CLIP_STRETCHING_H_