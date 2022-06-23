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

/// <summary>
/// Abstract stretch class
/// </summary>
class RS_TOOLSET_API StretchInterface {
 public:
  /// <summary>
  /// Run stretch on the given mats
  /// </summary>
  /// <param name="mats">The input and output given mats</param>
  /// <returns>Running state</returns>
  virtual bool Run(std::vector<cv::Mat>& mats) = 0;

  /// <summary>
  /// Add the single block statistic with the given band
  /// </summary>
  /// <param name="mat">The single block</param>
  /// <param name="band">The given band</param>
  /// <returns>Running state</returns>
  virtual bool AddSingleBlock(
      const cv::Mat& mat,
      int band) = 0;

  /// <summary>
  /// Add the multi block statistic
  /// </summary>
  /// <param name="mat">The multi block</param>
  /// <returns>Running state</returns>
  virtual bool AddMultiBlock(const cv::Mat& mat) = 0;
};

/// <summary>
/// Percent clip stretch class
/// </summary>
class RS_TOOLSET_API PercentClip : virtual public StretchInterface {
 public:
  /// <summary>
  /// Create the percent clip stretch shared pointer
  /// </summary>
  /// <param name="low_percent">The low percent trunction</param>
  /// <param name="high_percent">The high percent trunction</param>
  /// <returns>The output percent clip shared pointer</returns>
  static std::shared_ptr<PercentClip> Create(
      double low_percent,
      double high_percent);

  virtual void SetPercent(double low_percent, double high_percent) = 0;
  virtual void GetPercent(double& low_percent, double& high_percent) = 0;
};

/// <summary>
/// Standard deviation stretch class
/// </summary>
class RS_TOOLSET_API StandardDeviation :virtual public StretchInterface {
 public:
  /// <summary>
  /// Create the standard deviation stretch shared pointer
  /// </summary>
  /// <param name="scale">The trunction scale</param>
  /// <returns>The output standard deviation shared pointer</returns>
  static std::shared_ptr<StandardDeviation> Create(double scale);

  virtual void SetScale(double scale) = 0;
  virtual void GetScale(double& scale) = 0;
};

} // stretch 
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCH_H_