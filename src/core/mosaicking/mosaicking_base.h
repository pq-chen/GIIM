#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_

#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>

#include <rs-toolset/mosaicking.h>
#include <rs-toolset/color_balancing.h>

namespace rs_toolset {
namespace mosaicking {

class MosaickingBase : virtual public MosaickingInterface {
 public:
  explicit MosaickingBase(double tol);
  MosaickingBase(const MosaickingBase&) = delete;
  MosaickingBase& operator=(const MosaickingBase&) = delete;
  virtual ~MosaickingBase() { delete[] buffers_; }

  bool RunTaskForExisting(
      const std::string& path,
      OGRLayer* composite_table_layer,
      OGRLayer* border_layer,
      OGRGeometry* covered_border,
      std::vector<double>& borders_area,
      double rejection_ratio,
      int low_overview_trunc,
      int high_overview_trunc,
      const std::vector<int>& rgb_bands_map,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing) override;

  OGRGeometryUniquePtr RunTaskForPair(
      const std::string& path1,
      const std::string& path2,
      OGRGeometry* border1,
      OGRGeometry* border2,
      OGRGeometry* geometry1,
      OGRGeometry* geometry2,
      int low_overview_trunc,
      int high_overview_trunc,
      const std::vector<int>& rgb_bands_map,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing,
      int color_balancing_idx1,
      int color_balancing_idx2) override;
 
 protected:
  GDALDriver* memory_driver_;
  GDALDriver* mem_driver_;

 private:
  /**
   * @brief Create overlap datasets and the label dataset with the given overlap geometries for the existing
   * @param[in] factor The downsample factor
   * @param[in] geotrans The geotransform after reprojection
   * @param[in] rgb_bands_map The RGB bands' map
   * @param[in] source_raster_dataset The source raster dataset
   * @param[in] composite_table_layer The composite table layer
   * @param[in] covered_overlap_geometry The covered overlap geometry
   * @param[in] new_overlap_geometry The new overlap geometry
   * @param[in] color_balancing The color balancing shared pointer
   * @param[in] color_balancing_idx The source raster color balancing index
   * @param[out] covered_overlap_dataset The output covered overlap dataset
   * @param[out] new_overlap_dataset The output new overlap dataset
   * @param[in,out] label_raster_dataset The output label raster dataset if not nullptr
  */
  void CreateOverlapDatasets(
      double factor,
      double* geotrans,
      const std::vector<int>& rgb_bands_map,
      GDALDataset* source_raster_dataset,
      OGRLayer* composite_table_layer,
      OGRGeometry* covered_overlap_geometry,
      OGRGeometry* new_overlap_geometry,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing,
      int color_balancing_idx,
      GDALDatasetUniquePtr& covered_overlap_dataset,
      GDALDatasetUniquePtr& new_overlap_dataset,
      GDALDatasetUniquePtr& label_raster_dataset);

  /**
   * @brief Create overlap datasets and the label dataset with the given overlap geometries for the raster pair
   * @param[in] factor The downsample factor
   * @param[in] geotrans The geotransform after reprojection
   * @param[in] rgb_bands_map The RGB bands' map
   * @param[in] spatial_ref The spatial reference
   * @param[in] raster_dataset1 The first raster dataset
   * @param[in] raster_dataset2 The second raster dataset
   * @param[in] overlap_geometry1 The first overlap geometry
   * @param[in] overlap_geometry2 The second overlap geometry
   * @param[in] color_balancing The color balancing shared pointer
   * @param[in] color_balancing_idx1 The first raster color balancing index
   * @param[in] color_balancing_idx2 The second raster color balancing index
   * @param[out] overlap_dataset1 The first output overlap dataset
   * @param[out] overlap_dataset2 The second output overlap dataset
   * @param[in,out] label_raster_dataset The output label raster dataset if not nullptr
  */
  void CreateOverlapDatasets(
      double factor,
      double* geotrans,
      const std::vector<int>& rgb_bands_map,
      OGRSpatialReference* spatial_ref,
      GDALDataset* raster_dataset1,
      GDALDataset* raster_dataset2,
      OGRGeometry* overlap_geometry1,
      OGRGeometry* overlap_geometry2,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing,
      int color_balancing_idx1,
      int color_balancing_idx2,
      GDALDatasetUniquePtr& overlap_dataset1,
      GDALDatasetUniquePtr& overlap_dataset2,
      GDALDatasetUniquePtr& label_raster_dataset);

  /**
   * @brief Update the label raster and the overlap geometries
   * @param[in] idx The index of current operation
   * @param[in] covered_overlap_dataset The covered overlap dataset
   * @param[in] new_overlap_dataset The new overlap dataset
   * @param[in] valid_geometry The valid geometry
   * @param[in,out] label_raster_dataset The output label raster dataset
   * @param[in,out] covered_overlap_geometry The output covered overlap geometry except the last operation
   * @param[in,out] new_overlap_geometry The output new overlap geometry except the last operation
   * @param[in] last Whether is the last operation or not, default is false
   * @param[in] buffer_at_end The buffer distance in pixels at the end, default is 2.0
   * @return Running state
  */
  void UpdateMediums(
      int idx,
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      OGRGeometry* valid_geometry,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& covered_overlap_geometry,
      OGRGeometryUniquePtr& new_overlap_geometry,
      bool last = false,
      double buffer_at_end = 2.0);

  /**
   * @brief Prepare data for the mosaicking algorithm
   * @param[in] covered_overlap_dataset The covered overlap dataset
   * @param[in] new_overlap_dataset The new overlap dataset
   * @param[in] label_raster_dataset The label raster dataset
   * @param[out] covered_mat The output covered mat
   * @param[out] new_mat The output new mat
   * @param[out] label_mat The output label mat
  */
  virtual void PrepareData(
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      GDALDataset* label_raster_dataset,
      cv::Mat& covered_mat,
      cv::Mat& new_mat,
      cv::Mat& label_mat) = 0;

  /**
   * @brief Execute the mosaicking algorithm
   * @param[in] covered_mat The covered mat
   * @param[in] new_mat The new mat
   * @param[in] label_mat The label mat
   * @param[in] geotrans The geotransform of the above mats
   * @param[in] spatial_ref The spatial reference of the above mats
   * @param[out] label_raster_dataset The output label raster dataset
   * @param[out] label0_geometry The output label0 boundary
   * @param[out] label1_geometry The output label1 boundary
   * @note In the all-surrounded case, the label1_geometry should cover the label0_geometry
  */
  virtual void ExecuteMosaicking(
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& label_mat,
      double* geotrans,
      OGRSpatialReference* spatial_ref,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry) = 0;

  /**
   * @brief Update the source border, the composite table layer and the border layer
   * @param[in] covered_geometry The covered geometry
   * @param[in] new_geometry The new geometry
   * @param[in] covered_polygon The covered polygon
   * @param[in] rejection_ratio The regularization term to reject a new item in the output composite table
   * @param[in,out] source_border The source border
   * @param[in,out] composite_table_layer The composite table layer
   * @param[in,out] border_layer The border layer
   * @param[in,out] borders_area The borders' area
  */
  void UpdateResults(
      OGRGeometry* covered_geometry,
      OGRGeometry* new_geometry,
      OGRGeometry* covered_polygon,
      double rejection_ratio,
      OGRGeometryUniquePtr& source_border,
      OGRLayer* composite_table_layer,
      OGRLayer* border_layer,
      std::vector<double>& borders_area);

  double tol_;
  int* buffers_;
};

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_