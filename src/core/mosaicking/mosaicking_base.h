#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_

#include <string>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>

#include <rs-toolset/mosaicking.h>


namespace rs_toolset {
namespace mosaicking {

class MosaickingBase : virtual public MosaickingInterface {
 public:
  MosaickingBase(double tol);
  MosaickingBase(const MosaickingBase&) = delete;
  MosaickingBase& operator=(const MosaickingBase&) = delete;
  virtual ~MosaickingBase() { delete[] buffers_per_unit_; }

  bool Run(
      const std::string& raster_path,
      OGRLayer* composite_table_layer,
      OGRGeometry* covered_border,
      int last_overview_idx,
      bool use_seamline) override;
 
 protected:
  GDALDriver* memory_driver_;
  GDALDriver* mem_driver_;

 private:
  /**
   * @brief Update overlap datasets by the given overlap geometries
   * @param[in] downsample_factor The downsample factor
   * @param[in] geotrans The geotransform after reprojection
   * @param[in] composite_table_layer The composite table layer
   * @param[in] source_raster_dataset The source raster dataset
   * @param[in] covered_overlap_geometry The covered overlap geometry
   * @param[in] new_overlap_geometry The new overlap geometry
   * @param[out] covered_overlap_dataset The output covered overlap dataset
   * @param[out] new_overlap_dataset The output new overlap dataset
   * @param[out] label_raster_dataset The output label raster dataset if not nullptr
  */
  void UpdateOverlapDatasets(
      double downsample_factor,
      double* geotrans,
      OGRLayer* composite_table_layer,
      GDALDataset* source_raster_dataset,
      OGRGeometry* covered_overlap_geometry,
      OGRGeometry* new_overlap_geometry,
      GDALDatasetUniquePtr& covered_overlap_dataset,
      GDALDatasetUniquePtr& new_overlap_dataset,
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
  */
  void UpdateMediums(
      int idx,
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      OGRGeometry* valid_geometry,
      GDALDatasetUniquePtr& label_raster_dataset,
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
   * @param[in] geotrans The geotransform of the following datasets
   * @param[in] spatial_ref The spatial reference of the following datasets
   * @param[out] label_raster_dataset The output label raster dataset
   * @param[out] label0_geometry The output label0 boundary
   * @param[out] label1_geometry The output label1 boundary
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
   * @brief Create the seamline line string from the given seamline geometry
   * @param[in,out] seamline_geometry The given seamline geometry
   * @param[in,out] seamline_coors The coordinates of the seamline line sting
  */
  void CreateSeamlineLineString(
      OGRGeometryUniquePtr& seamline_geometry,
      std::vector<std::pair<double, double>>& seamline_coors);

  /**
   * @brief Rearrange points in label geometries
   * @param[in] seamline_coors The coordinates of the seamline line sting
   * @param[in,out] label0_geometry The label0 geometry
   * @param[in,out] label1_geometry The label1 geometry
  */
  void RearrangeLabelGeometries(
      const std::vector<std::pair<double, double>>& seamline_coors,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry);

  /**
   * @brief Simplify the label geometries
   * @param[in] tol The tolerance in pixels for simplifying the seamline
   * @param[in,out] seamline_line_string The seamlien line string
   * @param[in,out] label0_geometry The label0 geometry
   * @param[in,out] label1_geometry The label1 geometry
  */
  void SimplifyLabelGeometries(
      double tol,
      OGRGeometryUniquePtr& seamline_line_string,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry);

  /**
   * @brief Update the source border and the composite table
   * @param[in] covered_geometry The covered geometry
   * @param[in] new_geometry The new geometry
   * @param[in] covered_polygon The covered polygon
   * @param[in,out] source_border The source border
   * @param[in,out] composite_table_layer The composite table layer
  */
  void UpdateResults(
      OGRGeometry* covered_geometry,
      OGRGeometry* new_geometry,
      OGRGeometry* covered_polygon,
      OGRGeometryUniquePtr& source_border,
      OGRLayer* composite_table_layer);

  double tol_;
  int* buffers_per_unit_;
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_