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
  MosaickingBase();
  MosaickingBase(const MosaickingBase&) = delete;
  MosaickingBase& operator=(const MosaickingBase&) = delete;
  virtual ~MosaickingBase() = default;

  bool Run(
      const std::string& raster_path,
      OGRLayer* composite_table_layer,
      OGRGeometry* covered_border,
      bool use_seamline) override;
 
 protected:
  GDALDriver* memory_driver_;
  GDALDriver* mem_driver_;

 private:
  /**
   * @brief Create overlap datasets from the source raster dataset by the given overlap geometries
   * @param composite_table_layer: The composite table layer
   * @param source_raster_dataset: The source raster dataset
   * @param covered_overlap_geometry: The covered overlap geometry
   * @param new_overlap_geometry: The new overlap geometry
   * @param covered_overlap_dataset: The output covered overlap dataset
   * @param new_overlap_dataset: The output new overlap dataset
   * @param label_raster_dataset: The output label raster dataset if not nullptr
   * @param downsample_factor: The downsample factor
  */
  void CreateOverlapDatasets(
      OGRLayer* composite_table_layer,
      GDALDataset* source_raster_dataset,
      OGRGeometry* covered_overlap_geometry,
      OGRGeometry* new_overlap_geometry,
      GDALDatasetUniquePtr& covered_overlap_dataset,
      GDALDatasetUniquePtr& new_overlap_dataset,
      GDALDatasetUniquePtr& label_raster_dataset,
      double downsample_factor);

  /**
   * @brief Update the label raster dataset from the overlap datasets using the mosaicking method
   * @param covered_overlap_dataset: The covered overlap dataset
   * @param new_overlap_dataset: The new overlap dataset
   * @param label_raster_dataset: The output label raster dataset after resample
   * @param covered_overlap_geometry: The output covered overlap geometry except the last operation
   * @param new_overlap_geometry: The output new overlap geometry except the last operation
   * @param seamline_layer: The output seamline layer only in the last operation
   * @param overview_idx: The overview index
  */
  void UpdateLabelRasterDataset(
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_overlap_dataset,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& covered_overlap_geometry,
      OGRGeometryUniquePtr& new_overlap_geometry,
      OGRLayer* seamline_layer,
      int overview_idx);

  /// <summary>
  /// Prepare data for the mosaicking method
  /// </summary>
  /// <param name="covered_overlap_dataset">The covered overlap dataset</param>
  /// <param name="new_overlap_dataset">The new overlap dataset</param>
  /// <param name="label_raster_dataset">The label raster dataset</param>
  /// <param name="covered_mat">The output covered mat</param>
  /// <param name="new_mat">The output new mat</param>
  /// <param name="label_mat">The output label mat</param>
  virtual void PrepareData(
      GDALDataset* covered_overlap_dataset,
      GDALDataset* new_raster_dataset,
      GDALDataset* label_raster_dataset,
      cv::Mat& covered_mat,
      cv::Mat& new_mat,
      cv::Mat& label_mat) = 0;

  /// <summary>
  /// Create seamline geometries by the mosaicking method core algorithm
  /// </summary>
  /// <param name="new_overlap_dataset">The new overlap dataset</param>
  /// <param name="covered_mat">The output covered mat</param>
  /// <param name="new_mat">The output new mat</param>
  /// <param name="label_mat">The output label mat</param>
  /// <param name="label_raster_dataset">The output label raster dataset</param>
  /// <param name="label0_geometry">The output label0 geometry except the last operation</param>
  /// <param name="label1_geometry">The output label1 geometry except the last operation</param>
  /// <param name="seamline_layer">The output seamline layer only in the last operation</param>
  virtual void CreateSeamlineGeometries(
      GDALDataset* new_overlap_dataset,
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& label_mat,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry,
      OGRLayer* seamline_layer) = 0;

  /// <summary>
  /// Add the seamline to the composite table and the source border
  /// </summary>
  /// <param name="seamline_layer">The seamline layer</param>
  /// <param name="composite_table_layer">The composite table layer</param>
  /// <param name="source_border">The source border</param>
  void AddSeamline(
      OGRLayer* seamline_layer,
      OGRLayer* composite_table_layer,
      OGRGeometryUniquePtr& source_border);
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_