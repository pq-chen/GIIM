#ifndef RS_TOOLSET_SRC_CORE_STRETCH_STANDARD_DEVIATION_STRETCH_H_
#define RS_TOOLSET_SRC_CORE_STRETCH_STANDARD_DEVIATION_STRETCH_H_

#pragma warning(disable:4250)

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretch.h>
#include "stretch_base.h"


namespace rs_toolset {
namespace stretch {

class StandardDeviationImpl final 
    : public StretchBase, public StandardDeviation {
 public:
  StandardDeviationImpl(double scale) : scale_(scale) {
    spdlog::info(
        "Creating the standard deviation stretch with\nScale: {}", scale);
  }
  StandardDeviationImpl(const StandardDeviationImpl&) = delete;
  StandardDeviationImpl& operator=(const StandardDeviationImpl&) = delete;
  ~StandardDeviationImpl() = default;

  bool AddSingleBlock(
      const cv::Mat& mat,
      int band) override;
  bool AddMultiBlock(const cv::Mat& mat) override;
  void SetScale(double scale) override { scale_ = scale; }
  void GetScale(double& scale) override { scale = scale_; }

 protected:
  void CreateThresholds(
      std::vector<int>& low_thres,
      std::vector<int>& high_thres) override;

 private:
  std::vector<int> pixels_counts_; // Valid pixels count
  std::vector<double> sums_; // Valid pixels sums
  std::vector<double> square_sums_; // Valid pixels square sums

  double scale_;
};

} // stretch 
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_STRETCH_STANDARD_DEVIATION_STRETCH_H_