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

  bool RunTaskForSerial(
      const std::string& path,
      OGRLayer* mosaicking_layer,
      OGRLayer* border_layer,
      OGRGeometry* covered_border,
      std::vector<double>& borders_area,
      int low_overviews_trunc,
      int high_overviews_trunc,
      double surrounded_buffer,
      double rejection_ratio,
      bool anti_surrounded,
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
      int low_overviews_trunc,
      int high_overviews_trunc,
      const std::vector<int>& rgb_bands_map,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing,
      int color_balancing_idx1,
      int color_balancing_idx2) override;
 
 protected:
  GDALDriver* memory_driver_;
  GDALDriver* mem_driver_;

 private:
  enum class Status {
    OVERLAP,
    SURROUNDED,
    ANTI_SURROUNDED,
  };

  bool InitializeOverlapGeometries(
      double factor,
      double* geotrans,
      double surrounded_buffer,
      double rejection_ratio,
      bool anti_surrounded,
      OGRPolygon* covered_polygon,
      OGRGeometryUniquePtr& source_border,
      std::vector<OGRGeometryUniquePtr>& covered_overlap_geometries,
      std::vector<OGRGeometryUniquePtr>& new_overlap_geometries,
      std::vector<OGRGeometryUniquePtr>& source_borders,
      Status& status);

  /**
   * @brief Create overlap datasets and the label dataset with the given overlap geometries for the existing
   * @param[in] factor The downsample factor
   * @param[in] geotrans The geotransform after reprojection
   * @param[in] rgb_bands_map The RGB bands' map
   * @param[in] source_raster_path The source raster path
   * @param[in] mosaicking_layer The mosaicking layer
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
      const std::string& source_raster_path,
      OGRLayer* mosaicking_layer,
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
      const std::string& raster_path1,
      const std::string& raster_path2,
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
   * @brief Create seamlines between overlap_geometries
   * @param[in] geotrans geotrans The geotrans of covered or new overlap dataset
   * @param[in] covered_overlap_dataset The covered overlap dataset
   * @param[in] new_overlap_dataset The new overlap dataset
   * @param[in,out] label_raster_dataset The output label raster dataset
   * @param[in,out] covered_overlap_geometry The output covered overlap geometry
   * @param[in,out] new_overlap_geometry The output new overlap geometry
   * @param[in] swap Whether swaps label0 geometry and label1 geometry in the mosaicking, default is false
   * @param[in] subtasks_count The substasks' count, default is 0 means no subtask dispatch
   * @return Running state
  */
  bool CreateSeamlines(
      double* geotrans,
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& covered_overlap_geometry,
      OGRGeometryUniquePtr& new_overlap_geometry,
      bool swap = false,
      int subtasks_count = 0);

  /**
   * @brief Update the label raster and the overlap geometries
   * @param[in] serial Whether is the serial case or not
   * @param[in] idx The index of current operation
   * @param[in] last Whether is the last operation or not
   * @param[in] overlap Whether is the overlap case or not
   * @param[in] geotrans The geotrans of covered or new overlap dataset
   * @param[in] valid_geometry The valid geometry
   * @param[in,out] covered_overlap_geometry The output covered overlap geometry except the last operation
   * @param[in,out] new_overlap_geometry The output new overlap geometry except the last operation'
   * @return Running state
  */
  bool UpdateMediums(
      bool serial,
      int idx,
      bool last,
      bool overlap,
      double* geotrans,
      OGRGeometry* valid_geometry,
      OGRGeometryUniquePtr& covered_overlap_geometry,
      OGRGeometryUniquePtr& new_overlap_geometry);

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
   * @param[out] label0_geometry The output label0 geometry
   * @param[out] label1_geometry The output label1 geometry
   * @param[in] swap Whether swaps label0 geometry and label1 geometry in the mosaicking, default is false
   * @param[in] subtasks_count The substasks' count, default is 0 means no subtask dispatch
   * @return Running state
   * @note In the all-surrounded case, the label1 geometry should cover the label0 geometry
  */
  virtual bool ExecuteMosaicking(
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& label_mat,
      double* geotrans,
      OGRSpatialReference* spatial_ref,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry,
      bool swap = false,
      int subtasks_count = 0) = 0;

  /**
   * @brief Update the source border, the mosaicking layer and the border layer
   * @param[in] covered_geometry The covered geometry
   * @param[in] new_geometry The new geometry
   * @param[in] covered_polygon The covered polygon
   * @param[in] rejection_ratio The regularization term to reject a new item
   * @param[in,out] source_border The source border
   * @param[in,out] mosaicking_layer The mosaicking layer
   * @param[in,out] border_layer The border layer
   * @param[in,out] borders_area The borders' area
  */
  void UpdateResults(
      OGRGeometry* covered_geometry,
      OGRGeometry* new_geometry,
      OGRGeometry* covered_polygon,
      double rejection_ratio,
      OGRGeometryUniquePtr& source_border,
      OGRLayer* mosaicking_layer,
      OGRLayer* border_layer,
      std::vector<double>& borders_area);

  double tol_;
  int* buffers_;
};

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_