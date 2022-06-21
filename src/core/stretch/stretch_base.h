#ifndef RS_TOOLSET_SRC_CORE_STRETCH_STRETCH_BASE_H_
#define RS_TOOLSET_SRC_CORE_STRETCH_STRETCH_BASE_H_

#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>


namespace rs_toolset {
namespace stretch {

class StretchBase : virtual public StretchInterface {
 public:
  bool Run(std::vector<cv::Mat>& mats) override;

 protected:
   virtual void CreateThreshold(
      std::vector<int>& low_thres,
      std::vector<int>& high_thres) = 0;
};

} // stretch 
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_STRETCH_STRETCH_BASE_H_