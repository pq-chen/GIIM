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
  explicit GramSchmidtImpl(int block_size)
      : ComponentSubstitutionBase(block_size, false) {
    spdlog::info("Creating a Gram-Schmidt pansharpening");
  }
  GramSchmidtImpl(const GramSchmidtImpl&) = delete;
  GramSchmidtImpl& operator=(const GramSchmidtImpl&) = delete;
  ~GramSchmidtImpl() override = default;

 private:
  struct Statistics {
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
    return static_cast<void*>(new Statistics{
        0,
        0.0,
        0.0,
        std::vector<double>(bands_count, 0.0),
        std::vector<double>(bands_count, 0.0) });
  }

  std::vector<double> CreateWeights(void* statistics) override {
    auto s(static_cast<Statistics*>(statistics));
    int bands_count(static_cast<int>(s->upsampled_ms_sums.size()));
    return std::vector<double>(bands_count, 1.0 / bands_count);
  }

  void UpdateUpsampleInfo(
      const Data& data,
      const std::vector<double>& weights,
      void* statistics) override;

  std::vector<double> CreateInjectionGains(void* statistics) override;

  void DestroyStatistics(void* statistics) override {
    delete static_cast<Statistics*>(statistics);
  }
};

}  // namespace pansharpening
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_GRAM_SCHMIDT_H_