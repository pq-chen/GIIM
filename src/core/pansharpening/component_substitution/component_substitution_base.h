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
      const std::vector<int>& pansharpened_bands_map)override;

 protected:
  virtual void* CreateStatistic(
      int bands_count,
      int x_size,
      int y_size) = 0;

  virtual void UpdateDownsampleInfo(
      const Data& data,
      void* s) {}

  virtual std::vector<double> CreateWeights(void* s) = 0;

  virtual void UpdateStatistic(
      const Data& data,
      const std::vector<double>& weights,
      void* s) = 0;

  virtual std::vector<double> CreateInjectionGains(void* s) = 0;

  virtual std::vector<cv::Mat> CreateDeltaMats(
      const Data& data,
      const std::vector<double>& weights,
      void* s) = 0;

  virtual void DestroyStatistic(void* s) = 0;

 private:
  bool need_downsample_info_;   
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_PANSHARPENING_COMPONENT_SUBSTITUTION_COMPONENT_SUBSTITUTION_BASE_H_