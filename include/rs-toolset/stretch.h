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

/** @brief Abstract stretch class */
class RS_TOOLSET_API StretchInterface {
 public:
  /**
   * @brief Run the stretch algorithm on the given mats
   * @param[in] mats The input and output given mats
   * @param[in] nodata_value The nodata value, default is 0.0
   * @return Running state
  */
  virtual bool Run(std::vector<cv::Mat>& mats, double nodata_value = 0.0) = 0;

  /**
   * @brief Add statistic for the single-band block mat
   * @param[in] mat The single-band block mat
   * @param[in] band The given band index
   * @return Running state
  */
  virtual bool AddStatForSingleBlock(
      const cv::Mat& mat,
      int band) = 0;

  /**
   * @brief Add statistic for the multi-bands block mat
   * @param[in] mat The multi-band block mat
   * @return Running state
  */
  virtual bool AddStatForMultiBlock(const cv::Mat& mat) = 0;

  virtual void ClearStat() = 0;
};

/** @brief Percent clip stretch class implementing the percent clip stretch algorithm */
class RS_TOOLSET_API PercentClip : virtual public StretchInterface {
 public:
  /**
   * @brief Create a percent clip stretch shared pointer
   * @param[in] low_percent The low percent trunction, default is 0.0025
   * @param[in] high_percent The high percent trunction, default is 0.0025
   * @return The output percent clip shared pointer
  */
  static std::shared_ptr<PercentClip> Create(
      double low_percent = 0.0025,
      double high_percent = 0.0025);

  virtual void SetPercent(double low_percent, double high_percent) = 0;
  virtual void GetPercent(double& low_percent, double& high_percent) = 0;
};

/** @brief Standard deviation stretch class implementing the standard deviation stretch algorithm */
class RS_TOOLSET_API StandardDeviation :virtual public StretchInterface {
 public:
  /**
   * @brief Create a standard deviation stretch shared pointer
   * @param[in] scale The trunction scale, default is 2.5
   * @return The output standard deviation shared pointer
  */
  static std::shared_ptr<StandardDeviation> Create(double scale = 2.5);

  virtual void SetScale(double scale) = 0;
  virtual void GetScale(double& scale) = 0;
};

} // stretch
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCH_H_