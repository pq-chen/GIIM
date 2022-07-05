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
        "Creating a standard deviation stretch with\n - Scale: {}", scale);
  }
  StandardDeviationImpl(const StandardDeviationImpl&) = delete;
  StandardDeviationImpl& operator=(const StandardDeviationImpl&) = delete;
  ~StandardDeviationImpl() = default;

  bool AddStatForSingleBlock(
      const cv::Mat& mat,
      int band) override;
  bool AddStatForMultiBlock(const cv::Mat& mat) override;

  void SetScale(double scale) override { scale_ = scale; }
  void GetScale(double& scale) override { scale = scale_; }

  void ClearStat() {
    pixels_counts_.resize(0);
    sums_.resize(0);
    square_sums_.resize(0);
  }

 protected:
  void CreateThres(
      std::vector<int>& low_thres,
      std::vector<int>& high_thres) override;

 private:
  std::vector<int> pixels_counts_;
  std::vector<double> sums_;
  std::vector<double> square_sums_;

  double scale_;
};

} // stretch 
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_STRETCH_STANDARD_DEVIATION_STRETCH_H_