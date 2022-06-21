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
      OGRLayer* composite_table_layer,
      OGRGeometry* covered_border,
      const std::string& raster_path,
      bool with_seamline) override;
 
 protected:
  GDALDriver* memory_driver_;
  GDALDriver* mem_driver_;

 private:
  void CreateOverlapDatasets(
      OGRLayer* composite_table_layer,
      GDALDataset* source_raster_dataset,
      OGRGeometry* covered_overlap_geometry,
      OGRGeometry* new_overlap_geometry,
      GDALDatasetUniquePtr& covered_overlap_dataset,
      GDALDatasetUniquePtr& new_overlap_dataset,
      GDALDatasetUniquePtr& label_raster_dataset,
      double downsample_factor);

  void UpdateOverlapGeometries(
      GDALDataset* covered_raster_dataset,
      GDALDataset* new_raster_dataset,
      GDALDatasetUniquePtr& label_raster_dataset,
      OGRGeometryUniquePtr& covered_overlap_geometry,
      OGRGeometryUniquePtr& new_overlap_geometry,
      OGRLayer* seamline_layer,
      int overview_idx);

  virtual void PrepareData(
      GDALDataset* covered_raster_dataset,
      GDALDataset* new_raster_dataset,
      GDALDataset* seamline_raster_dataset,
      cv::Mat& covered_mat,
      cv::Mat& new_mat,
      cv::Mat& label_mat) = 0;

  virtual void CreateSeamlineGeometries(
      GDALDataset* new_raster_dataset,
      const cv::Mat& covered_mat,
      const cv::Mat& new_mat,
      const cv::Mat& label_mat,
      GDALDatasetUniquePtr& seamline_raster_dataset,
      OGRGeometryUniquePtr& label0_geometry,
      OGRGeometryUniquePtr& label1_geometry,
      OGRLayer* seamline_layer) = 0;

  void AddSeamline(
      OGRLayer* composite_table_layer,
      OGRLayer* seamline_layer,
      OGRGeometryUniquePtr& source_border);
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_BASE_H_