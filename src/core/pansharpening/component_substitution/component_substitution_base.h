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
  explicit ComponentSubstitutionBase(int block_size, bool need_downsample_info)
      : PansharpeningBase(block_size),
        need_downsample_info_(need_downsample_info) {}
  ComponentSubstitutionBase(const ComponentSubstitutionBase&) = delete;
  ComponentSubstitutionBase& operator=(const ComponentSubstitutionBase&) =
      delete;
  ~ComponentSubstitutionBase() override = default;

  bool Run(
      const std::string& pan_path,
      const std::string& ms_path,
      const std::string& output_path,
      bool use_rpc,
      bool use_stretching,
      const std::vector<int>& bands_map) override;

 protected:
  /**
   * @brief Create a statistics struct
   * @param[in] bands_count The bands' count
   * @param[in] x_size The MS raster x size
   * @param[in] y_size The MS raster y size
   * @return The output statistics struct
  */
  virtual void* CreateStatistics(
      int bands_count,
      int ms_x_size,
      int ms_y_size) = 0;

  /**
   * @brief Update the downsample information in the statistics struct with the given downsampled data if needed
   * @param[in] data The given downsampled data
   * @param[in,out] statistics The statistics struct
  */
  virtual void UpdateDownsampleInfo(const Data& data, void* statistics) {}

  /**
   * @brief Create weights with the statistics struct for building the synthetic low resolution PAN mat
   * @param[in] statistics The statistics struct
   * @return The output weights
  */
  virtual std::vector<double> CreateWeights(void* statistics) = 0;

  /**
   * @brief Update the upsample information in the statistics struct with the given upsampled data and weights
   * @param data[in] The weights
   * @param weights[in] The given upsampled data
   * @param statistics[in,out] The statistics struct
  */
  virtual void UpdateUpsampleInfo(
      const Data& data,
      const std::vector<double>& weights,
      void* statistics) = 0;

  /**
   * @brief Create injection gains with the statistics struct for building the pansharpened mat
   * @param statistics[in] The statistics struct
   * @return The output injection gains
  */
  virtual std::vector<double> CreateInjectionGains(void* statistics) = 0;

  /**
   * @brief Destroy the statistics struct
   * @param[in,out] statistics The statistics struct
  */
  virtual void DestroyStatistics(void* statistics) = 0;

  cv::Mat hist_matching_mat_;

 private:
  /**
   * @brief Create delta mats between the PAN mat and the synthetic low resolution PAN mat for all bands
   * @param data[in] The given upsampled data
   * @param weights[in] The weights
   * @return The output delta mats
  */
  std::vector<cv::Mat> CreateDeltaMats(
      const Data& data,
      const std::vector<double>& weights);

  bool need_downsample_info_;
};

}  // namespace pansharpening
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_COMPONENT_SUBSTITUTION_BASE_H_