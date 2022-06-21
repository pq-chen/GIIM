#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_H_ 
#define RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_H_

#pragma warning(disable:4250)

#include <vector>

#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include "component_substitution_base.h"
#include <rs-toolset/pansharpening.h>


namespace rs_toolset {
namespace pansharpening {

class GramSchmidtImpl final 
    : public ComponentSubstitutionBase, public GramSchmidt {
 public:
  GramSchmidtImpl(int block_size)
      : ComponentSubstitutionBase(block_size, false) {}
  GramSchmidtImpl(const GramSchmidtImpl&) = delete;
  GramSchmidtImpl& operator=(const GramSchmidtImpl&) = delete;
  ~GramSchmidtImpl() override = default;

 private:
  struct Statistic {
    int pixels_count;
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
        0,
        0.,
        0.,
        std::vector<double>(bands_count, 0.),
        std::vector<double>(bands_count, 0.),
        std::vector<double>(bands_count, 0.) });
  }

  std::vector<double> CreateWeights(void* s) override {
    spdlog::debug("Creating the weights");
    auto _s(static_cast<Statistic*>(s));
    int bands_count(static_cast<int>(_s->upsampled_ms_means.size()));
    spdlog::info("Creating the weights - done");
    return std::vector<double>(bands_count, 1. / bands_count);
  }

  void UpdateStatistic(
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

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_H_