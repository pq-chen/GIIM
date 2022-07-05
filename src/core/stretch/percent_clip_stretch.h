#ifndef RS_TOOLSET_SRC_CORE_STRETCH_PERCENT_CLIP_STRETCH_H_
#define RS_TOOLSET_SRC_CORE_STRETCH_PERCENT_CLIP_STRETCH_H_

#pragma warning(disable:4250)

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretch.h>
#include "stretch_base.h"


namespace rs_toolset {
namespace stretch {

class PercentClipImpl final : public StretchBase, public PercentClip {
 public:
  PercentClipImpl(
      double low_percent,
      double high_percent) 
      : low_percent_(low_percent), high_percent_(high_percent) {
    spdlog::info(
        "Creating a percent clip stretch with\n"
        " - Low percent: {}\n"
        " - High percent: {}", low_percent_, high_percent_);
  }
  PercentClipImpl(const PercentClipImpl&) = delete;
  PercentClipImpl& operator=(const PercentClipImpl&) = delete;
  ~PercentClipImpl() = default;

  bool AddStatForSingleBlock(
      const cv::Mat& mat,
      int band) override;
  bool AddStatForMultiBlock(const cv::Mat& mat) override;

  void SetPercent(double low_percent, double high_percent) override {
    low_percent_ = low_percent;
    high_percent_ = high_percent;
  }
  void GetPercent(double& low_percent, double& high_percent) override {
    low_percent = low_percent_;
    high_percent = high_percent_;
  }

  void ClearStat() { hist_mats_.resize(0); }

 protected:
  void CreateThres(
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