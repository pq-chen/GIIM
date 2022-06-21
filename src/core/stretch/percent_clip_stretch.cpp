#include "percent_clip_stretch.h"

#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>
#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace stretch {

void PercentClipImpl::AddBlock(
    const cv::Mat& mat,
    int idx) {
  if (hist_mats_.size() <= idx) {
    hist_mats_.push_back(utils::CalcHist(mat, cv::Mat())[0]);
  } else {
    hist_mats_[idx] += utils::CalcHist(mat, cv::Mat())[0];
  }
}

void PercentClipImpl::CreateThreshold(
    std::vector<int>& low_thres,
    std::vector<int>& high_thres) {
  int bands_count(static_cast<int>(hist_mats_.size()));
  low_thres.resize(bands_count);
  high_thres.resize(bands_count);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++) {
    hist_mats_[b].at<float>(0) = 0;
    cv::normalize(hist_mats_[b], hist_mats_[b], 1.0, 0.0, cv::NORM_L1);
    double sum(hist_mats_[b].at<float>(0));
    for (int i = 1; i < hist_mats_[b].rows; i++) {
      sum += hist_mats_[b].at<float>(i);
      if (sum > low_percent_) {
        low_thres[b] = i;
        break;
      }
    }
    sum = hist_mats_[b].at<float>(hist_mats_[b].rows - 1);
    for (int i = hist_mats_[b].rows - 2; i >= 0; i--) {
      sum += hist_mats_[b].at<float>(i);
      if (sum > high_percent_) {
        high_thres[b] = i;
        break;
      }
    }
  }
  hist_mats_.resize(0);
}

std::shared_ptr<PercentClip> PercentClip::Create(
    double low_percent,
    double high_percent) {
  return std::make_shared<PercentClipImpl>(low_percent, high_percent);
}

} // stretch 
} // rs_toolset