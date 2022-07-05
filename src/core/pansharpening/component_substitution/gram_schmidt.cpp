#include "gram_schmidt.h"

#include <memory>
#include <vector>

#include <opencv2/opencv.hpp>

#include <rs-toolset/pansharpening.h>
#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace pansharpening {

void GramSchmidtImpl::UpdateUpsampleInfo(
    const Data& data,
    const std::vector<double>& weights,
    void* s) {
  spdlog::debug("Updating the statistic struct");
  auto _s(static_cast<Statistic*>(s));
  int bands_count(static_cast<int>(data.mats.size()));

  // Create a synthetic low resolution PAN mat and calculate its mean
  cv::Mat synthetic_low_reso_pan_mat;
  data.mats[0].convertTo(synthetic_low_reso_pan_mat, CV_16SC1, weights[0]);
  for (int i = 1; i < bands_count; i++)
    synthetic_low_reso_pan_mat += weights[i] * data.mats[i];
  synthetic_low_reso_pan_mat.convertTo(
      synthetic_low_reso_pan_mat, data.mat.type());

  // Update histogram mats
  std::vector<cv::Mat> 
      cur_pan_hist_mat(utils::CreateHist(data.mat)),
      cur_synthetic_low_reso_pan_hist_mat(utils::CreateHist(
          synthetic_low_reso_pan_mat));
  if (_s->pan_hist_mat.empty()) {
    _s->pan_hist_mat.push_back(cur_pan_hist_mat[0]);
    _s->synthetic_low_reso_pan_hist_mat.push_back(
        cur_synthetic_low_reso_pan_hist_mat[0]);
  } else {
    _s->pan_hist_mat[0] += cur_pan_hist_mat[0];
    _s->synthetic_low_reso_pan_hist_mat[0] +=
        cur_synthetic_low_reso_pan_hist_mat[0];
  }

  // Update the other statistics
  _s->pixels_count += cv::countNonZero(data.mat);
  _s->synthetic_low_reso_pan_sum_ += cv::sum(synthetic_low_reso_pan_mat)[0];
  _s->synthetic_low_reso_pan_square_sum_ +=
      synthetic_low_reso_pan_mat.dot(synthetic_low_reso_pan_mat);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < bands_count; b++) {
    cv::Mat temp1, temp2;
    synthetic_low_reso_pan_mat.convertTo(temp1, CV_32SC1);
    data.mats[b].convertTo(temp2, CV_32SC1);
    _s->upsampled_ms_sums_[b] += cv::sum(data.mats[b])[0];
    _s->product_sums_[b] += cv::sum(temp1.mul(temp2))[0];
  }
  spdlog::info("Updating the statistic struct - done");
}

std::vector<double> GramSchmidtImpl::CreateInjectionGains(void* s) {
  spdlog::debug("Creating injection gains");
  auto _s(static_cast<Statistic*>(s));
  int bands_count(static_cast<int>(_s->upsampled_ms_sums_.size()));
  std::vector<double> injection_gains(bands_count);
  double var(
      _s->pixels_count * _s->synthetic_low_reso_pan_square_sum_ -
      _s->synthetic_low_reso_pan_sum_ * _s->synthetic_low_reso_pan_sum_);
#pragma omp parallel for schedule(static, bands_count)
  for (int b = 0; b < injection_gains.size(); b++)
    injection_gains[b] = 
        (_s->pixels_count * _s->product_sums_[b] -
        _s->synthetic_low_reso_pan_sum_ * _s->upsampled_ms_sums_[b]) / var;
  _s->hist_matching_mat =  utils::CreateHistMatchingLut(
      _s->pan_hist_mat, _s->synthetic_low_reso_pan_hist_mat);
  spdlog::info("Creating injection gains - done");
  return injection_gains;
}

std::vector<cv::Mat> GramSchmidtImpl::CreateDeltaMats(
    const Data& data,
    const std::vector<double>& weights,
    void* s) {
  spdlog::debug("Creating delta mats");
  auto _s(static_cast<Statistic*>(s));
  int bands_count(static_cast<int>(data.mats.size()));
  cv::Mat delta_mat;

  // The delta mat minus the synthetic low resolution PAN mat
  data.mats[0].convertTo(delta_mat, CV_16SC1, -weights[0]);
  for (int i = 1; i < bands_count; i++)
    delta_mat -= weights[i] * data.mats[i];
  if (weights.size() == bands_count + 1) {
    cv::Mat constant_mat;
    cv::threshold(
        data.mat, constant_mat, 0, weights.back(), cv::THRESH_BINARY);
    delta_mat -= constant_mat;
  }

  // The delta mat plus the histogram matched PAN mat
  delta_mat += utils::TransformMatWithLut(data.mat, _s->hist_matching_mat);
  spdlog::info("Creating delta mats - done");
  return std::vector<cv::Mat>(bands_count, delta_mat);
}

std::shared_ptr<GramSchmidt> GramSchmidt::Create(int block_size) {
  return std::make_shared<GramSchmidtImpl>(block_size);
}

} // pansharpening
} // rs_toolset
