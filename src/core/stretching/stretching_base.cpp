#include "stretching_base.h"

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

namespace rs_toolset {
namespace stretching {

bool StretchingBase::Run(std::vector<cv::Mat>& mats) {
  spdlog::info("Running a stretching task");
  if (depth_ == -1) {
    spdlog::error("No statistics in the stretching class");
    return false;
  }
  if (mats.size() != bands_count_) {
    spdlog::error(
        "The input mats size {} must be the same with "
        "the statistics bands' count {}", mats.size(), bands_count_);
    return false;
  }
  for (const auto& mat : mats) {
    if (mat.depth() != depth_ || mat.channels() != 1) {
      spdlog::error(
          "One of the input mats is invalid for the existing statistics");
      return false;
    }
  }
  if (lut_mats_.empty()) {
    spdlog::error(
        "Please call \"CreateLutMats()\" before run the stretching algorithm");
    return false;
  }
  
  cv::Mat mask = cv::Mat::zeros(mats[0].size(), CV_8UC1);
  for (int i = 0; i < bands_count_; ++i) {
    mask.setTo(1, mats[i] != 0);
  }
  cv::threshold(mask, mask, 0, 1, cv::THRESH_BINARY);
  for (int b(0); b < bands_count_; ++b) {
    cv::Mat mat(mats[b].size(), CV_8UC1);
    if (mats[b].depth() == CV_8U) {
      cv::LUT(mats[b], lut_mats_[b], mat);
    } else {
#pragma omp parallel for schedule(dynamic)
      for (int row(0); row < mats[b].rows; ++row) {
        for (int col(0); col < mats[b].cols; ++col) {
          mat.at<uint8_t>(row, col) = lut_mats_[b].at<uint8_t>(
              0, mats[b].at<uint16_t>(row, col));
        }
      }
    }
    mats[b] = (mat + 1).mul(mask);
  }
  spdlog::info("Running a stretching task - done");
  return true;
}

bool StretchingBase::Run(cv::Mat& mat) {
  std::vector<cv::Mat> mats;
  cv::split(mat, mats);
  if (Run(mats)) {
    cv::merge(mats, mat);
    return true;
  } else {
    return false;
  }
}

void StretchingBase::Clear() {
  depth_ = -1;
  bands_count_ = -1;
  lut_mats_.clear();
}

}  // namespace stretching
}  // namespace rs_toolset