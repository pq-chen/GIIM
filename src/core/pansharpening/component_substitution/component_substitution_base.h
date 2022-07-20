#ifndef RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_COMPONENT_SUBSTITUTION_BASE_H_ 
#define RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_COMPONENT_SUBSTITUTION_BASE_H_

#include <string>
#include <vector>

#include <opencv2/opencv.hpp>

#include "../pansharpening_base.h"


namespace rs_toolset {
namespace pansharpening {

class ComponentSubstitutionBase : public PansharpeningBase {
 public:
  ComponentSubstitutionBase(
      int block_size,
      bool need_downsample_info)
      : need_downsample_info_(need_downsample_info),
        PansharpeningBase(block_size) {}
  ComponentSubstitutionBase(const ComponentSubstitutionBase&) = delete;
  ComponentSubstitutionBase& operator=(const ComponentSubstitutionBase&) =
      delete;
  ~ComponentSubstitutionBase() override = default;

  bool Run(
      const std::string& pan_path,
      const std::string& ms_path,
      const std::string& pansharpened_path,
      bool use_rpc,
      bool use_stretch,
      const std::vector<int>& pansharpened_bands_map) override;

 protected:
  /**
   * @brief Create a statistic struct
   * @param[in] bands_count The bands' count
   * @param[in] x_size x size(any kind)
   * @param[in] y_size y size(any kind)
   * @return The output statistic struct
  */
  virtual void* CreateStatistic(
      int bands_count,
      int x_size,
      int y_size) = 0;

  /**
   * @brief Update the downsample information in the statistic struct with the given downsampled data if needed
   * @param[in] data The given downsampled data
   * @param[in,out] s The statistic struct
  */
  virtual void UpdateDownsampleInfo(
      const Data& data,
      void* s) {}

  /**
   * @brief Create weights with the statistic struct for building the synthetic low resolution PAN mat
   * @param[in] s The statistic struct
   * @return The output weights
  */
  virtual std::vector<double> CreateWeights(void* s) = 0;

  /**
   * @brief Update the upsample information in the statistic struct with the given upsampled data and weights
   * @param data[in] Weights
   * @param weights[in] The given upsampled data
   * @param s[in,out] The statistic struct
  */
  virtual void UpdateUpsampleInfo(
      const Data& data,
      const std::vector<double>& weights,
      void* s) = 0;

  /**
   * @brief Create injection gains with the statistic struct for building the pansharpened mat
   * @param s[in] The statistic struct
   * @return The output injection gains
  */
  virtual std::vector<double> CreateInjectionGains(void* s) = 0;

  /**
   * @brief Create delta mats between the PAN mat and the synthetic low resolution PAN mat for all bands
   * @param data[in] The given upsampled data
   * @param weights[in] Weights
   * @param s[in] The statistic struct
   * @return The output delta mats
  */
  virtual std::vector<cv::Mat> CreateDeltaMats(
      const Data& data,
      const std::vector<double>& weights,
      void* s) = 0;

  /**
   * @brief Destroy the statistic struct
   * @param[in,out] s The statistic struct
  */
  virtual void DestroyStatistic(void* s) = 0;

 private:
  bool need_downsample_info_; // Whether needs downsample information
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_COMPONENT_SUBSTITUTION_BASE_H_