#ifndef RS_TOOLSET_SRC_CORE_STRETCH_STRETCH_BASE_H_
#define RS_TOOLSET_SRC_CORE_STRETCH_STRETCH_BASE_H_

#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>


namespace rs_toolset {
namespace stretch {

class StretchBase : virtual public StretchInterface {
 public:
  bool Run(std::vector<cv::Mat>& mats, double nodata_value) override;

 protected:
  /**
   * @brief Create low and high thresholds for all bands
   * @param[out] low_thres The output low thresholds
   * @param[out] high_thres The output high thresholds
  */
  virtual void CreateThres(
      std::vector<int>& low_thres,
      std::vector<int>& high_thres) = 0;
};

} // stretch
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_STRETCH_STRETCH_BASE_H_