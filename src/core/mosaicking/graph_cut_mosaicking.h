#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_

#pragma warning(disable:4250)

#include <cstdint>

#include <map>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include "mosaicking_base.h"
#include <rs-toolset/mosaicking.h>


namespace rs_toolset {
namespace mosaicking {

class GraphCutImpl final : public MosaickingBase, public GraphCut {
 public:
  GraphCutImpl(
      float grad_self_low,
      float grad_self_high,
      float grad_self_exp,
      float diff_low,
      float diff_exp,
      double tol)
      : MosaickingBase(tol),
        grad_self_low_(grad_self_low),
        grad_self_high_(grad_self_high),
        grad_self_exp_(grad_self_exp),
        diff_low_(diff_low),
        diff_exp_(diff_exp) {
    spdlog::info(
        "Creating the graph cut mosaicking with\n"
        " - Gradient-self term low trunction: {}\n"
        " - Gradient-self term low trunction: {}\n"
        " - Gradient-self term exponential: {}\n"
        " - Difference term low trunction: {}\n"
        " - Difference term exponential: {}",
        grad_self_low, grad_self_high, grad_self_exp, diff_low, diff_exp);
  }
  GraphCutImpl(const GraphCutImpl&) = delete;
  GraphCutImpl& operator=(const GraphCutImpl&) = delete;
  ~GraphCutImpl() = default;

 private:
  struct GraphCutInfo {
    int idx;
    uint16_t label;
  };
  struct Data {
    float grad_self_low;
    float grad_self_high;
    float grad_self_exp;
    float diff_low;
    float diff_exp;
    cv::Mat covered_mat;
    cv::Mat covered_x_mat;
    cv::Mat covered_y_mat;
    cv::Mat new_mat;
    cv::Mat new_x_mat;
    cv::Mat new_y_mat;
    std::vector<std::pair<int, int>>* idxes_to_coors;
  };

  void PrepareData(
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      GDALDataset* label_raster_dataset,
      cv::Mat& covered_mat,
      cv::Mat& new_mat,
      cv::Mat& label_mat) override;

  void ExecuteMosaicking(
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& label_mat,
      double* geotrans,
      OGRSpatialReference* spatial_ref,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry) override;

  cv::Mat CreateLightnessMat(GDALDataset* dataset);

  static int SmoothCost(
      int site1,
      int site2,
      int label1,
      int label2,
      void* data);

  float grad_self_low_;
  float grad_self_high_;
  float grad_self_exp_;
  float diff_low_;
  float diff_exp_;
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_