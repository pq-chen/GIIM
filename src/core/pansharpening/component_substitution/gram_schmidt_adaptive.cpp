#include "gram_schmidt_adaptive.h"

#include <cstdint>
#include <memory>
#include <vector>

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

#include <rs-toolset/pansharpening.h>
#include <rs-toolset/utils.hpp>

namespace rs_toolset {
namespace pansharpening {

void GramSchmidtAdaptiveImpl::UpdateDownsampleInfo(
    const Data& data,
    void* statistics) {
  auto s(static_cast<Statistic*>(statistics));
#pragma omp parallel for schedule(dynamic)
  for (int row(0); row < data.pan_mat.rows; ++row) {
    for (int col(0); col < data.pan_mat.cols; ++col) {
      int idx(row * data.pan_mat.cols + col);
      if (data.pan_mat.depth() == CV_8U) {
        for (int b(0); b < data.ms_mats.size(); ++b)
          s->A(idx, b) = data.ms_mats[b].at<uint8_t>(row, col);
        s->b(idx) = data.pan_mat.at<uint8_t>(row, col);
      } else {
        for (int b(0); b < data.ms_mats.size(); ++b)
          s->A(idx, b) = data.ms_mats[b].at<uint16_t>(row, col);
        s->b(idx) = data.pan_mat.at<uint16_t>(row, col);
      }
    }
  }
}

std::vector<double> GramSchmidtAdaptiveImpl::CreateWeights(void* statistics) {
  auto s(static_cast<Statistic*>(statistics));
  Eigen::VectorXf solve(
      s->A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(s->b));
  std::vector<double> weights(solve.rows());
  for (int b(0); b < solve.rows(); ++b)
    weights[b] = solve(b);
  return weights;
}

void GramSchmidtAdaptiveImpl::UpdateUpsampleInfo(
    const Data& data,
    const std::vector<double>& weights,
    void* statistics) {
  auto s(static_cast<Statistic*>(statistics));
  auto bands_count(static_cast<int>(data.ms_mats.size()));

  // Create a synthetic low resolution PAN mat and calculate its mean
  cv::Mat synthetic_low_reso_pan_mat;
  data.ms_mats[0].convertTo(synthetic_low_reso_pan_mat, CV_16SC1, weights[0]);
  for (int i(1); i < bands_count; ++i)
    synthetic_low_reso_pan_mat += weights[i] * data.ms_mats[i];
  synthetic_low_reso_pan_mat.convertTo(
      synthetic_low_reso_pan_mat, data.pan_mat.type());

  // Update histogram mats
  cv::Mat
      cur_pan_hist_mat(utils::CreateHists(data.pan_mat)[0]),
      cur_synthetic_low_reso_pan_hist_mat(utils::CreateHists(
          synthetic_low_reso_pan_mat)[0]);
  if (s->pan_hist_mat.empty()) {
    s->pan_hist_mat = cur_pan_hist_mat;
    s->synthetic_low_reso_pan_hist_mat = cur_synthetic_low_reso_pan_hist_mat;
  } else {
    s->pan_hist_mat += cur_pan_hist_mat;
    s->synthetic_low_reso_pan_hist_mat += cur_synthetic_low_reso_pan_hist_mat;
  }

  // Update the other statistics
  s->pixels_count += cv::countNonZero(data.pan_mat);
  s->synthetic_low_reso_pan_sum += cv::sum(synthetic_low_reso_pan_mat)[0];
  s->synthetic_low_reso_pan_square_sum +=
      synthetic_low_reso_pan_mat.dot(synthetic_low_reso_pan_mat);
#pragma omp parallel for schedule(static, bands_count)
  for (int b(0); b < bands_count; ++b) {
    cv::Mat temp1, temp2;
    synthetic_low_reso_pan_mat.convertTo(temp1, CV_32SC1);
    data.ms_mats[b].convertTo(temp2, CV_32SC1);
    s->upsampled_ms_sums[b] += cv::sum(data.ms_mats[b])[0];
    s->product_sums[b] += cv::sum(temp1.mul(temp2))[0];
  }
}

std::vector<double> GramSchmidtAdaptiveImpl::CreateInjectionGains(
    void* statistics) {
  auto s(static_cast<Statistic*>(statistics));
  auto bands_count(static_cast<int>(s->upsampled_ms_sums.size()));
  std::vector<double> injection_gains(bands_count);
  double var(
      s->pixels_count * s->synthetic_low_reso_pan_square_sum -
          s->synthetic_low_reso_pan_sum * s->synthetic_low_reso_pan_sum);
#pragma omp parallel for schedule(static, bands_count)
  for (int b(0); b < injection_gains.size(); ++b) {
    injection_gains[b] =
        (s->pixels_count * s->product_sums[b] -
            s->synthetic_low_reso_pan_sum * s->upsampled_ms_sums[b]) / var;
  }
  hist_matching_mat_ = utils::CreateHistMatchingLut(
      {s->pan_hist_mat}, {s->synthetic_low_reso_pan_hist_mat});
  return injection_gains;
}

std::shared_ptr<GramSchmidtAdaptive> GramSchmidtAdaptive::Create(
    int block_size) {
  return std::make_shared<GramSchmidtAdaptiveImpl>(block_size);
}

}  // namespace pansharpening
}  // namespace rs_toolset