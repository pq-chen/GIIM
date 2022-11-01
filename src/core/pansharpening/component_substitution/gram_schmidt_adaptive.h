#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_ADAPTIVE_H_ 
#define RS_TOOLSES_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_ADAPTIVE_H_

#pragma warning(disable:4250)

#include <vector>

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include "component_substitution_base.h"
#include <rs-toolset/pansharpening.h>

namespace rs_toolset {
namespace pansharpening {

class GramSchmidtAdaptiveImpl final 
    : public ComponentSubstitutionBase, public GramSchmidtAdaptive {
 public:
  explicit GramSchmidtAdaptiveImpl(int block_size)
      : ComponentSubstitutionBase(block_size, true) {
    spdlog::info("Creating a Gram-Schmidt adaptive pansharpening");
  }
  GramSchmidtAdaptiveImpl(const GramSchmidtAdaptiveImpl&) = delete;
  GramSchmidtAdaptiveImpl& operator=(const GramSchmidtAdaptiveImpl&) = delete;
  ~GramSchmidtAdaptiveImpl() override = default;

 private:
  struct Statistic {
    Eigen::MatrixXf A;
    Eigen::VectorXf b;
    int pixels_count;
    double synthetic_low_reso_pan_sum;
    double synthetic_low_reso_pan_square_sum;
    std::vector<double> upsampled_ms_sums;
    std::vector<double> product_sums;
    cv::Mat pan_hist_mat;
    cv::Mat synthetic_low_reso_pan_hist_mat;
  };

  void* CreateStatistics(
      int bands_count,
      int ms_x_size,
      int ms_y_size) override {
    return static_cast<void*>(new Statistic{
        Eigen::MatrixXf::Ones(ms_x_size * ms_y_size, bands_count + 1),
        Eigen::VectorXf(ms_x_size * ms_y_size),
        0,
        0.0,
        0.0,
        std::vector<double>(bands_count, 0.0),
        std::vector<double>(bands_count, 0.0) });
  }

  void UpdateDownsampleInfo(const Data& data, void* statistics) override;

  std::vector<double> CreateWeights(void* statistics) override;

  void UpdateUpsampleInfo(
      const Data& data,
      const std::vector<double>& weights,
      void* statistics) override;

  std::vector<double> CreateInjectionGains(void* statistics) override;

  void DestroyStatistics(void* s) override {
    delete static_cast<Statistic*>(s);
  }
};

} // namespace pansharpening
} // namespace rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_ADAPTIVE_H_