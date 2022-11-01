#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCHING_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCHING_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else  // RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif  // RS_TOOLSET_EXPORTS

#else  // _WIN32
#define RS_TOOLSET_API
#endif  // _WIN32

#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>

namespace rs_toolset {
namespace stretching {

/** @brief Abstract stretching class */
class RS_TOOLSET_API StretchingInterface {
 public:
  /**
   * @brief Run the stretching algorithm on the splited mats
   * @param[in,out] mats The splited mats, available depths include CV_8U and CV_16U
   * @return Running state
  */
  virtual bool Run(std::vector<cv::Mat>& mats) = 0;

  /**
   * @brief Run the stretching algorithm on the multi-bands mat
   * @param[in,out] mats The multi-bands mat, available depths include CV_8U and CV_16U
   * @return Running state
  */
  virtual bool Run(cv::Mat& mat) = 0;

  /**
   * @brief Accumulate statistics for the the splited mats
   * @param[in] mats The splited mats, available depths include CV_8U and CV_16U
   * @return Running state
  */
  virtual bool AccumulateStat(const std::vector<cv::Mat>& mats) = 0;

  /**
   * @brief Accumulate statistics for the multi-bands mat
   * @param[in] mat The multi-bands mat, available depths include CV_8U and CV_16U
   * @return Running state
  */
  virtual bool AccumulateStat(const cv::Mat& mat) = 0;

  /**
   * @brief Create LUT mats for the exsiting statistics
   * @return Running state
  */
  virtual bool CreateLutMats() = 0;

  /** @brief Clear statistics and thresholds in the stretching class */
  virtual void Clear() = 0;
};

/** @brief Percent clip stretching class implementing the percent clip stretching algorithm */
class RS_TOOLSET_API PercentClip : virtual public StretchingInterface {
 public:
  /**
   * @brief Create a percent clip stretching shared pointer
   * @param[in] low_percent The low percent trunction, default is 0.005
   * @param[in] high_percent The high percent trunction, default is 0.005
   * @return The output percent clip shared pointer
  */
  static std::shared_ptr<PercentClip> Create(
      double low_percent = 0.005,
      double high_percent = 0.005);

  virtual void SetPercent(double low_percent, double high_percent) = 0;

  virtual void GetPercent(double& low_percent, double& high_percent) = 0;
};

/** @brief Standard deviations stretching class implementing the standard deviations stretching algorithm */
class RS_TOOLSET_API StandardDeviations : virtual public StretchingInterface {
 public:
  /**
   * @brief Create a standard deviations stretching shared pointer
   * @param[in] scale The trunction scale, default is 2.5
   * @return The output standard deviations shared pointer
  */
  static std::shared_ptr<StandardDeviations> Create(double scale = 2.5);

  virtual void SetScale(double scale) = 0;

  virtual void GetScale(double& scale) = 0;
};

}  // namespace stretching
}  // namespace rs_toolset

#endif  // RS_TOOLSET_INCLUDE_RS_TOOLSET_STRETCHING_H_