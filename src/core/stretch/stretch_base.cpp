#include "stretch_base.h"

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>


namespace rs_toolset {
namespace stretch {

bool StretchBase::Run(std::vector<cv::Mat>& mats, double nodata_value) {
  spdlog::info("Running a stretch task");

  if (mats[0].depth() != CV_8U && mats[0].depth() != CV_16U) {
    spdlog::warn(
        "Input mats' type is neither 8-bit unsigned nor 16-bit unsigned");
    return false;
  }
  
  // Create a mask mat for the given mats
  cv::Mat mask;
  cv::Mat((cv::abs(mats[0] - nodata_value))).convertTo(mask, CV_8UC1);
  cv::threshold(mask, mask, 0, 1, cv::THRESH_BINARY);

  // Run the stretch algorithm on the given mats
  std::vector<int> low_thres, high_thres;
  CreateThres(low_thres, high_thres);
  int bands_count(static_cast<int>(mats.size()));
  std::vector<double> means(bands_count);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++) {
    cv::Mat((mats[b] - low_thres[b]) / (high_thres[b] - low_thres[b]) * 254.0)
        .convertTo(mats[b], CV_8UC1);
    mats[b] = (mats[b] + 1).mul(mask);
  }
  spdlog::info("Running a stretch task - done");
  return true;
}

} // stretch
} // rs_toolset