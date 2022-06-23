#include "stretch_base.h"

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>


namespace rs_toolset {
namespace stretch {

bool StretchBase::Run(std::vector<cv::Mat>& mats) {
  spdlog::info("Running the stretch task");
  std::vector<int> low_thres, high_thres;
  CreateThresholds(low_thres, high_thres);
  int bands_count(static_cast<int>(mats.size()));
  std::vector<double> means(bands_count);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++)
    cv::Mat((mats[b] - low_thres[b]) / (high_thres[b] - low_thres[b]) * 255.)
        .convertTo(mats[b], CV_8UC1);
  return true;
  spdlog::info("Running the stretch task - done");
}

} // stretch 
} // rs_toolset