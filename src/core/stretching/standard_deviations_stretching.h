#ifndef RS_TOOLSET_SRC_CORE_STRETCHING_STANDARD_DEVIATIONS_STRETCHING_H_
#define RS_TOOLSET_SRC_CORE_STRETCHING_STANDARD_DEVIATIONS_STRETCHING_H_

#pragma warning(disable:4250)

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretching.h>
#include "stretching_base.h"

namespace rs_toolset {
namespace stretching {

class StandardDeviationsImpl final 
    : public StretchingBase, public StandardDeviations {
 public:
  explicit StandardDeviationsImpl(double scale) : scale_(scale) {
    spdlog::info(
        "Creating a standard deviations stretching with\n - Scale: {}", scale);
  }
  StandardDeviationsImpl(const StandardDeviationsImpl&) = delete;
  StandardDeviationsImpl& operator=(const StandardDeviationsImpl&) = delete;
  ~StandardDeviationsImpl() = default;

  bool AccumulateStat(const std::vector<cv::Mat>& mats) override;

  bool AccumulateStat(const cv::Mat& mat) override;

  bool CreateLutMats() override;

  void Clear() override;

  void SetScale(double scale) override { scale_ = scale; }

  void GetScale(double& scale) override { scale = scale_; }

 private:
  double scale_;
  std::vector<int> pixels_counts_;
  std::vector<double> sums_;
  std::vector<double> square_sums_;
};

}  // namespace stretching
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_STRETCHING_STANDARD_DEVIATIONS_STRETCHING_H_