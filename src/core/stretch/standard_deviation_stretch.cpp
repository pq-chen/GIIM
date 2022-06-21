#include "standard_deviation_stretch.h"

#include <cmath>

#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/stretch.h>


namespace rs_toolset {
namespace stretch {

void StandardDeviationImpl::AddBlock(
    const cv::Mat& mat,
    int idx) {
  if (pixels_counts_.size() <= idx) {
    pixels_counts_.push_back(cv::countNonZero(mat));
    sums_.push_back(cv::sum(mat)[0]);
    square_sums_.push_back(mat.dot(mat));
  } else {
    pixels_counts_[idx] += cv::countNonZero(mat);
    sums_[idx] += cv::sum(mat)[0];
    square_sums_[idx] += mat.dot(mat);
  }
}

void StandardDeviationImpl::CreateThreshold(
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
  pixels_counts_.resize(0);
  sums_.resize(0);
  square_sums_.resize(0);
}

std::shared_ptr<StandardDeviation> StandardDeviation::Create(
    double n) {
  return std::make_shared<StandardDeviationImpl>(n);
}

} // stretch 
} // rs_toolset