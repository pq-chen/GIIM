#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCH_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCH_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else // LCMAKE_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif // LCMAKE_EXPORTS

#else // _WIN32
#define RS_TOOLSET_API
#endif // _WIN32

#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>


namespace rs_toolset {
namespace stretch {

class RS_TOOLSET_API StretchInterface {
 public:
  virtual bool Run(std::vector<cv::Mat>& mats) = 0;
  virtual void AddBlock(
      const cv::Mat& mat,
      int idx) = 0;
};

class RS_TOOLSET_API PercentClip : virtual public StretchInterface {
 public:
  static std::shared_ptr<PercentClip> Create(
      double low_percent,
      double high_percent);

  virtual void SetPercent(double low_percent, double high_percent) = 0;
  virtual void GetPercent(double& low_percent, double& high_percent) = 0;
};

class RS_TOOLSET_API StandardDeviation :virtual public StretchInterface {
 public:
  static std::shared_ptr<StandardDeviation> Create(double n);

  virtual void SetScale(double scale) = 0;
  virtual void GetScale(double& scale) = 0;
};

} // stretch 
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCH_H_