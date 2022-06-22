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
      double grad_exp,
      double min_diff,
      double max_diff)
      : grad_exp_(grad_exp), min_diff_(min_diff), max_diff_(max_diff) {
    spdlog::info(
        "Creating the graph cut mosaicking with\nGradient exponential: {}\n"
        "Minimum difference: {}\nMaximum difference: {}", grad_exp, min_diff,
        max_diff);
  }
  GraphCutImpl(const GraphCutImpl&) = delete;
  GraphCutImpl& operator=(const GraphCutImpl&) = delete;
  ~GraphCutImpl() = default;

 private:
  struct GraphCutInfo {
    int idx;
    uint16_t label;
  };
  struct SmoothExtraData {
    double grad_exp;
    double min_diff;
    double max_diff;
    cv::Mat covered_mat;
    cv::Mat covered_x_mat;
    cv::Mat covered_y_mat;
    cv::Mat new_mat;
    cv::Mat new_x_mat;
    cv::Mat new_y_mat;
    std::map<int, std::pair<int, int>> valid_idx_to_coors;
  };

  void PrepareData(
      GDALDataset* covered_raster_dataset,
      GDALDataset* new_raster_dataset,
      GDALDataset* seamline_raster_dataset,
      cv::Mat& covered_mat,
      cv::Mat& new_mat,
      cv::Mat& label_mat) override;

  void CreateSeamlineGeometries(
      GDALDataset* new_raster_dataset,
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& label_mat,
      GDALDatasetUniquePtr& seamline_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry,
      OGRLayer* seamline_layer) override;

  cv::Mat CreateLumiData(GDALDataset* dataset);

  static int SmoothCost(
      int site1,
      int site2,
      int label1,
      int label2,
      void* extra_data);

  double grad_exp_;
  double min_diff_;
  double max_diff_;
};

} // mosaicking  
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_MOSAICKING_GRAPH_CUT_MOSAICKING_H_