#include "standard_deviations_stretching.h"

#include <cmath>
#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretching.h>

namespace rs_toolset {
namespace stretching {

bool StandardDeviationsImpl::AccumulateStat(const std::vector<cv::Mat>& mats) {
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
    pixels_counts_.resize(bands_count_);
    sums_.resize(bands_count_);
    square_sums_.resize(bands_count_);
#pragma omp parallel for schedule(static, bands_count_)
    for (int b(0); b < bands_count_; ++b) {
      pixels_counts_[b] = cv::countNonZero(mats[b]);
      sums_[b] = cv::sum(mats[b])[0];
      square_sums_[b] = mats[b].dot(mats[b]);
    }
    return true;
  } else if (depth == depth_ && mats.size() == bands_count_) {
#pragma omp parallel for schedule(static, bands_count_)
    for (int b(0); b < bands_count_; ++b) {
      pixels_counts_[b] += cv::countNonZero(mats[b]);
      sums_[b] += cv::sum(mats[b])[0];
      square_sums_[b] += mats[b].dot(mats[b]);
    }
    return true;
  } else {
    spdlog::error("The input mats are invalid for the existing statistics");
    return false;
  }
}

bool StandardDeviationsImpl::AccumulateStat(const cv::Mat& mat) {
  std::vector<cv::Mat> mats;
  cv::split(mat, mats);
  return AccumulateStat(mats);
}

bool StandardDeviationsImpl::CreateLutMats() {
  spdlog::debug("Creating thresholds");
  if (depth_ == -1) {
    spdlog::error("No statistics in the stretching class");
    return false;
  }

  lut_mats_.clear();
  lut_mats_.resize(bands_count_);
#pragma omp parallel for schedule(static, bands_count_)
  for (int b(0); b < bands_count_; ++b) {
    auto mean(sums_[b] / pixels_counts_[b]),
        stddev(sqrt(square_sums_[b] / pixels_counts_[b] - mean * mean));
    auto low_threshold(static_cast<int>(round(mean - scale_ * stddev))),
        high_threshold(static_cast<int>(round(mean + scale_ * stddev)));
    lut_mats_[b] = cv::Mat(1, depth_ == CV_8U ? 256 : 65536, CV_16UC1);
    for (int i(0); i < lut_mats_[b].cols; ++i)
      lut_mats_[b].at<uint16_t>(i) = i;
    cv::Mat(254.0 * (lut_mats_[b] - low_threshold) /
        (high_threshold - low_threshold)).convertTo(lut_mats_[b], CV_8UC1);
  }
  spdlog::info("Creating thresholds - done");
}

void StandardDeviationsImpl::Clear() {
  StretchingBase::Clear();
  pixels_counts_.clear();
  sums_.clear();
  square_sums_.clear();
}

std::shared_ptr<StandardDeviations> StandardDeviations::Create(double scale) {
  return std::make_shared<StandardDeviationsImpl>(scale);
}

}  // namespace stretching
}  // namespace rs_toolset