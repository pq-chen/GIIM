#include "percent_clip_stretch.h"

#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>
#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace stretch {

bool PercentClipImpl::AddSingleBlock(
    const cv::Mat& mat,
    int band) {
  if (mat.channels() != 1) return false;
  if (hist_mats_.size() <= band) {
    hist_mats_.push_back(utils::CalcHist(mat)[0]);
  } else {
    hist_mats_[band] += utils::CalcHist(mat)[0];
  }
  return true;
}

bool PercentClipImpl::AddMultiBlock(const cv::Mat& mat) {
  if (hist_mats_.size() != 0 && hist_mats_.size() != mat.channels())
    return false;
  if (!hist_mats_.size()) {
    hist_mats_ = utils::CalcHist(mat);
  } else {
    std::vector<cv::Mat> hist_mats(utils::CalcHist(mat));
    for (int b = 0; b < mat.channels(); b++)
      hist_mats_[b] += hist_mats[b];
  }
  return true;
}

void PercentClipImpl::CreateThresholds(
    std::vector<int>& low_thres,
    std::vector<int>& high_thres) {
  int bands_count(static_cast<int>(hist_mats_.size()));
  low_thres.resize(bands_count);
  high_thres.resize(bands_count);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++) {
    hist_mats_[b].at<float>(0) = 0.0;
    cv::normalize(hist_mats_[b], hist_mats_[b], 1.0, 0.0, cv::NORM_L1);
    double sum(0.0);
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

  // Clear histogram mats after creating thresholds
  hist_mats_.resize(0);
}

std::shared_ptr<PercentClip> PercentClip::Create(
    double low_percent,
    double high_percent) {
  return std::make_shared<PercentClipImpl>(low_percent, high_percent);
}

} // stretch 
} // rs_toolset