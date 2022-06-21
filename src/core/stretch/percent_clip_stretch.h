#ifndef RS_TOOLSET_SRC_CORE_STRETCH_PERCENT_CLIP_STRETCH_H_
#define RS_TOOLSET_SRC_CORE_STRETCH_PERCENT_CLIP_STRETCH_H_

#pragma warning(disable:4250)

#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>
#include "stretch_base.h"


namespace rs_toolset {
namespace stretch {

class PercentClipImpl final : public StretchBase, public PercentClip {
 public:
  PercentClipImpl(
      double low_percent,
      double high_percent) 
      : low_percent_(low_percent),
        high_percent_(high_percent) {}
  PercentClipImpl(const PercentClipImpl&) = delete;
  PercentClipImpl& operator=(const PercentClipImpl&) = delete;
  ~PercentClipImpl() = default;

  void AddBlock(
      const cv::Mat& mat,
      int idx) override;

  void SetPercent(double low_percent, double high_percent) override {
    low_percent_ = low_percent;
    high_percent_ = high_percent;
  }
  void GetPercent(double& low_percent, double& high_percent) override {
    low_percent = low_percent;
    high_percent = high_percent_;
  }

 protected:
  void CreateThreshold(
      std::vector<int>& low_thres,
      std::vector<int>& high_thres) override;

 private:
  std::vector<cv::Mat> hist_mats_;

  double low_percent_;
  double high_percent_;
};

} // stretch 
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_STRETCH_PERCENT_CLIP_STRETCH_H_