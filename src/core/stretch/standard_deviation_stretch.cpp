#include "standard_deviation_stretch.h"

#include <cmath>

#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>


namespace rs_toolset {
namespace stretch {

bool StandardDeviationImpl::AddSingleBlock(
    const cv::Mat& mat,
    int band) {
  if (mat.channels() != 1) return false;
  if (pixels_counts_.size() <= band) {
    pixels_counts_.push_back(cv::countNonZero(mat));
    sums_.push_back(cv::sum(mat)[0]);
    square_sums_.push_back(mat.dot(mat));
  } else {
    pixels_counts_[band] += cv::countNonZero(mat);
    sums_[band] += cv::sum(mat)[0];
    square_sums_[band] += mat.dot(mat);
  }
  return true;
}

bool StandardDeviationImpl::AddMultiBlock(const cv::Mat& mat) {
  if (pixels_counts_.size() != 0 && pixels_counts_.size() != mat.channels())
    return false;
  std::vector<cv::Mat> mats;
  cv::split(mat, mats);
  for (int b = 0; b < mats.size(); b++)
    AddSingleBlock(mats[b], b);
  return true;
}

void StandardDeviationImpl::CreateThresholds(
    std::vector<int>& low_thres,
    std::vector<int>& high_thres) {
  int bands_count(static_cast<int>(pixels_counts_.size()));
  std::vector<double> means(sums_);
  std::vector<double> stddevs(square_sums_);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++) {
    means[b] /= pixels_counts_[b];
    stddevs[b] /= pixels_counts_[b];
    stddevs[b] -= means[b] * means[b];
    stddevs[b] = sqrt(stddevs[b]);
  }
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++) {
    low_thres.push_back(static_cast<int>(round(
        means[b] - scale_ * stddevs[b])));
    high_thres.push_back(static_cast<int>(round(
        means[b] + scale_ * stddevs[b])));
  }

  // Clear statistic after creating thresholds
  pixels_counts_.resize(0);
  sums_.resize(0);
  square_sums_.resize(0);
}

std::shared_ptr<StandardDeviation> StandardDeviation::Create(
    double scale) {
  return std::make_shared<StandardDeviationImpl>(scale);
}

} // stretch 
} // rs_toolset