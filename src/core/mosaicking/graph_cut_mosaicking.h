#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_

#pragma warning(disable:4250)

#include <cstdint>
#include <map>
#include <vector>
#include <utility>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>

#include "mosaicking_base.h"
#include <rs-toolset/mosaicking.h>

namespace rs_toolset {
namespace mosaicking {

class GraphCutImpl final : public MosaickingBase, public GraphCut {
 public:
  explicit GraphCutImpl(
      float grad_self_low,
      float grad_self_high,
      float grad_self_exp,
      float diff_low,
      float diff_exp,
      double tol);
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
    cv::Mat covered_mask_mat;
    cv::Mat covered_x_mat;
    cv::Mat covered_y_mat;
    cv::Mat new_mat;
    cv::Mat new_mask_mat;
    cv::Mat new_x_mat;
    cv::Mat new_y_mat;
    std::vector<std::pair<int, int>>* idx_to_coor;
  };

  void PrepareData(
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      GDALDataset* covered_mask_dataset,
      GDALDataset* new_mask_dataset,
      GDALDataset* label_raster_dataset,
      cv::Mat& covered_mat,
      cv::Mat& new_mat,
      cv::Mat& covered_mask_mat,
      cv::Mat& new_mask_mat,
      cv::Mat& label_mat,
      bool connection_analysis) override;

  bool ExecuteMosaicking(
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& covered_mask_mat,
      const cv::Mat& new_mask_mat,
      const cv::Mat& label_mat,
      double* geotrans,
      OGRSpatialReference* spatial_ref,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry,
      bool swap,
      int subtasks_count) override;

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

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_