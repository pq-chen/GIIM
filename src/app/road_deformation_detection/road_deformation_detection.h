#ifndef RS_TOOLSET_SRC_APP_ROAD_DEFORMATION_DETECTION_ROAD_DEFORMATION_DETECTION_H_
#define RS_TOOLSET_SRC_APP_ROAD_DEFORMATION_DETECTION_ROAD_DEFORMATION_DETECTION_H_

#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>


namespace rs_toolset {
namespace road_deformation_detection {

class RoadDeformationDetection {
 public:
  static std::shared_ptr<RoadDeformationDetection> Create(
      const std::string& road_path,
      const std::string& road_layer_name,
      const std::string& dem_path,
      double ideal_sample_space_in_pixels,
      double threshold,
      int min_continual_count);
  RoadDeformationDetection(const RoadDeformationDetection&) = delete;
  RoadDeformationDetection& operator=(
      const RoadDeformationDetection&) = delete;
  ~RoadDeformationDetection();

  GDALDatasetUniquePtr Run(
      const std::string& dom_path,
      const std::string& rpb_path,
      const std::string& output_path);

 private:
  struct LineSampleInfo {
    std::vector<OGRPoint*> points;
    std::vector<double> normal_bias_ratios;
  };

  RoadDeformationDetection(
      GDALDataset* road_dataset,
      OGRLayer* road_layer,
      GDALDataset* dem_dataset,
      double* dem_geotrans,
      float dem_nodata_value,
      double ideal_sample_space_in_pixels,
      double threshold,
      int min_continual_count);

  void CalcLightAngles(
      const GDALRPCInfo& rpc_info,
      double& hori_angle,
      double& vert_shift);

  void DetectLineString(
      OGRLineString* line_string,
      OGRGeometry* border,
      cv::Mat dem_mat,
      double* geotrans,
      double hori_angle,
      double vert_shift,
      OGRLayer* layer);

  void SampleBetweenPoints(
      const OGRPoint& point1,
      const OGRPoint& point2,
      cv::Mat dem_mat,
      double* geotrans,
      double hori_angle,
      double vert_shift,
      LineSampleInfo& info);

  GDALDataset* road_dataset_;
  OGRLayer* road_layer_;
  GDALDataset* dem_dataset_;
  double dem_geotrans_[6];
  float dem_nodata_value_;
  double ideal_sample_space_in_pixels_;
  double threshold_;
  int min_continual_count_;
};

} // namespace road_deformation_detection
} // namespace rs_toolset

#endif // RS_TOOLSET_SRC_APP_ROAD_DEFORMATION_DETECTION_ROAD_DEFORMATION_DETECTION_H_