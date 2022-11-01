#include "percent_clip_stretching.h"

#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretching.h>
#include <rs-toolset/utils.hpp>

namespace rs_toolset {
namespace stretching {

bool PercentClipImpl::AccumulateStat(const std::vector<cv::Mat>& mats) {
  auto depth(mats[0].depth());
  for (const auto& mat : mats) {
    if (mat.depth() != depth || mat.channels() != 1) {
      spdlog::error("The input mats are inconsistent");
      return false;
    }
  }
  if (depth_ == -1 && (depth == CV_8U || depth == CV_16U)) {
    depth_ = depth;
    bands_count_ = static_cast<int>(mats.size());
    hist_mats_.resize(bands_count_);
#pragma omp parallel for schedule(static, bands_count_)
    for (int b(0); b < bands_count_; ++b)
      hist_mats_[b] = utils::CreateHists(mats[b])[0];
    return true;
  } else if (depth == depth_ && mats.size() == bands_count_) {
#pragma omp parallel for schedule(static, bands_count_)
    for (int b(0); b < bands_count_; ++b)
      hist_mats_[b] += utils::CreateHists(mats[b])[0];
    return true;
  } else {
    spdlog::error("The input mats are invalid for the existing statistics");
    return false;
  }
}

bool PercentClipImpl::AccumulateStat(const cv::Mat& mat) {
  std::vector<cv::Mat> mats;
  cv::split(mat, mats);
  return AccumulateStat(mats);
}

bool PercentClipImpl::CreateLutMats() {
  spdlog::debug("Calculating thresholds");
  if (depth_ == -1) {
    spdlog::error("No statistics in the stretching class");
    return false;
  }

  lut_mats_.clear();
  lut_mats_.resize(bands_count_);
#pragma omp parallel for schedule(static, bands_count_)
  for (int b(0); b < bands_count_; ++b) {
    hist_mats_[b].at<float>(0) = 0.0f;
    cv::normalize(hist_mats_[b], hist_mats_[b], 1.0, 0.0, cv::NORM_L1);
    int low_threshold(0), high_threshold(hist_mats_[b].rows - 1);
    float sum(0.0f);
    while (sum < low_percent_)
      sum += hist_mats_[b].at<float>(++low_threshold);
    sum = hist_mats_[b].at<float>(hist_mats_[b].rows - 1);
    while (sum < high_percent_)
      sum += hist_mats_[b].at<float>(--high_threshold);
    lut_mats_[b] = cv::Mat(1, depth_ == CV_8U ? 256 : 65536, CV_16UC1);
    for (int i(0); i < lut_mats_[b].cols; ++i)
      lut_mats_[b].at<uint16_t>(i) = i;
    cv::Mat(254.0 * (lut_mats_[b] - low_threshold) /
        (high_threshold - low_threshold)).convertTo(lut_mats_[b], CV_8UC1);
  }
  spdlog::info("Calculating thresholds - done");
  return true;
}

void PercentClipImpl::Clear() {
  StretchingBase::Clear();
  hist_mats_.clear();
}

void PercentClipImpl::SetPercent(double low_percent, double high_percent) {
  low_percent_ = low_percent;
  high_percent_ = high_percent;
}

void PercentClipImpl::GetPercent(double& low_percent, double& high_percent) {
  low_percent = low_percent_;
  high_percent = high_percent_;
}

std::shared_ptr<PercentClip> PercentClip::Create(
    double low_percent,
    double high_percent) {
  return std::make_shared<PercentClipImpl>(low_percent, high_percent);
}

}  // namespace stretching 
}  // namespace rs_toolset