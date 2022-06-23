#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_ADAPTIVE_H_ 
#define RS_TOOLSES_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_ADAPTIVE_H_

#pragma warning(disable:4250)

#include <vector>

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

#include "component_substitution_base.h"
#include <rs-toolset/pansharpening.h>


namespace rs_toolset {
namespace pansharpening {

class GramSchmidtAdaptiveImpl final 
    : public ComponentSubstitutionBase, public GramSchmidtAdaptive {
 public:
  GramSchmidtAdaptiveImpl(int block_size)
      : ComponentSubstitutionBase(block_size, true) {}
  GramSchmidtAdaptiveImpl(const GramSchmidtAdaptiveImpl&) = delete;
  GramSchmidtAdaptiveImpl& operator=(const GramSchmidtAdaptiveImpl&) = delete;
  ~GramSchmidtAdaptiveImpl() override = default;

 private:
  struct Statistic {
    Eigen::MatrixXf A;
    Eigen::VectorXf b;
    int upsample_pixels_count;
    double synthetic_low_reso_pan_sum_;
    double synthetic_low_reso_pan_square_sum_;
    std::vector<double> upsampled_ms_sums_;
    std::vector<double> product_sums_;
    std::vector<double> upsampled_ms_means;
    std::vector<cv::Mat> pan_hist_mat;
    std::vector<cv::Mat> synthetic_low_reso_pan_hist_mat;
  };

  void* CreateStatistic(
      int bands_count,
      int x_size,
      int y_size) override {
    return static_cast<void*>(new Statistic{
        Eigen::MatrixXf::Ones(x_size * y_size, bands_count + 1),
        Eigen::VectorXf(x_size * y_size),
        0,
        0.,
        0.,
        std::vector<double>(bands_count, 0.),
        std::vector<double>(bands_count, 0.),
        std::vector<double>(bands_count, 0.) });
  }

  void UpdateDownsampleInfo(
      const Data& data,
      void* s) override;

  std::vector<double> CreateWeights(void* s) override;

  void UpdateUpsampleInfo(
      const Data& data,
      const std::vector<double>& weights,
      void* s) override;

  std::vector<double> CreateInjectionGains(void* s) override;

  std::vector<cv::Mat> CreateDeltaMats(
      const Data& data,
      const std::vector<double>& weights,
      void* s) override;

  void DestroyStatistic(void* s) override {
    delete static_cast<Statistic*>(s);
  }
};

} // pansharpening
} // rs_toolset

#endif // #define RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_ADAPTIVE_H_