#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_CONTAINER_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_CONTAINER_H_

#pragma warning(disable:4250)

#include <cmath>

#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "mosaicking_base.h"
#include <rs-toolset/mosaicking.h>


namespace rs_toolset {
namespace mosaicking {

class MosaickingContainerImpl final
    : virtual public MosaickingContainerInterface, public MosaickingContainer {
 public:
  MosaickingContainerImpl(
      std::shared_ptr<MosaickingInterface> mosaicking,
      OGRSpatialReference* spatial_ref);
  MosaickingContainerImpl(
      std::shared_ptr<MosaickingInterface> mosaicking,
      OGRLayer* composite_table_layer,
      const std::string& rasters_dir);
  MosaickingContainerImpl(const MosaickingContainerImpl&) = delete;
  MosaickingContainerImpl& operator=(const MosaickingContainerImpl&) = delete;
  virtual ~MosaickingContainerImpl();

  void SortRasters(
      std::vector<std::string>& rasters_path,
      const SortFunc& sort_func) override;

  bool AddTask(
      const std::string& raster_path,
      int last_overview_idx,
      bool use_seamline) override;

  bool ExportCompositeTableVector(
      const std::string& composit_table_path,
      double buffer,
      double tol) override;

  std::vector<std::string> ExportAllRastersName() override;

  bool CreateMosaickingRaster(
      const std::string& mosaicking_raster_path,
      const std::string& composit_table_path,
      const std::string& rasters_dir,
      double reso,
      double blend_dist) override;
 
 private:
  std::shared_ptr<MosaickingInterface> mosaicking_;

  GDALDataset* composite_table_dataset_;
  OGRLayer* composite_table_layer_;
  OGRGeometry* covered_border_; // type is multipolygon
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_CONTAINER_H_