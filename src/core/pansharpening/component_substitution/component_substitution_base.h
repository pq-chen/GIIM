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
  /// <summary>
  /// Create the statistic struct
  /// </summary>
  /// <param name="bands_count">The bands' count</param>
  /// <param name="x_size">x size(any kind)</param>
  /// <param name="y_size">y size(any kind)</param>
  /// <returns>The output statistic struct</returns>
  virtual void* CreateStatistic(
      int bands_count,
      int x_size,
      int y_size) = 0;

  /// <summary>
  /// Update the downsample information in the statistic struct with the given downsampled data if needed
  /// </summary>
  /// <param name="data">The given downsampled data</param>
  /// <param name="s">The statistic struct</param>
  virtual void UpdateDownsampleInfo(
      const Data& data,
      void* s) {}

  /// <summary>
  /// Create weights with the statistic struct for building the synthetic low resolution PAN mat
  /// </summary>
  /// <param name="s">The statistic struct</param>
  /// <returns>The output weights</returns>
  virtual std::vector<double> CreateWeights(void* s) = 0;

  /// <summary>
  /// Update the upsample information in the statistic struct with the given upsampled data and weights
  /// </summary>
  /// <param name="data">Weights</param>
  /// <param name="weights">The given upsampled data</param>
  /// <param name="s">The statistic struct</param>
  virtual void UpdateUpsampleInfo(
      const Data& data,
      const std::vector<double>& weights,
      void* s) = 0;

  /// <summary>
  /// Create injection gains with the statistic struct for building the pansharpened mat
  /// </summary>
  /// <param name="s">The statistic struct</param>
  /// <returns>The output injection gains</returns>
  virtual std::vector<double> CreateInjectionGains(void* s) = 0;

  /// <summary>
  /// Create delta mats between the PAN mat and the synthetic low resolution PAN mat for all bands
  /// </summary>
  /// <param name="data">The given upsampled data</param>
  /// <param name="weights">weights</param>
  /// <param name="s">The statistic struct</param>
  /// <returns>The output delta mats</returns>
  virtual std::vector<cv::Mat> CreateDeltaMats(
      const Data& data,
      const std::vector<double>& weights,
      void* s) = 0;

  /// <summary>
  /// Destroy the statistic struct
  /// </summary>
  /// <param name="s">The statistic struct</param>
  virtual void DestroyStatistic(void* s) = 0;

 private:
  bool need_downsample_info_; // Whether needs downsample information
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_COMPONENT_SUBSTITUTION_BASE_H_