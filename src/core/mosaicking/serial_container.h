#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_SERIAL_CONTAINER_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_SERIAL_CONTAINER_H_

#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/mosaicking.h>

namespace rs_toolset {
namespace mosaicking {

class SerialContainerImpl final : public SerialContainer {
 public:
  explicit SerialContainerImpl(
      const std::shared_ptr<MosaickingInterface>& mosaicking,
      OGRSpatialReference* spatial_ref,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing);
  SerialContainerImpl(const SerialContainerImpl&) = delete;
  SerialContainerImpl& operator=(const SerialContainerImpl&) = delete;
  virtual ~SerialContainerImpl();

  /**
   * @brief Initialize the serial container by the external composite table
   * @param[in] composite_tabel_layer The external composite tabel layer
   * @param[in] rasters_dir The rasters' directory corresponding to the external composite table
  */
  void InitializeByExt(
      OGRLayer* composite_table_layer,
      const std::string& rasters_dir);

  bool SortRasters(
      std::vector<std::string>& paths,
      const SortFunc& sort_func) override;

  bool AddTask(
      const std::string& path,
      int low_overviews_trunc,
      int high_overviews_trunc,
      double surrounded_buffer,
      double rejection_ratio,
      bool anti_surrounded,
      const std::vector<int>& rgb_bands_map) override;

  GDALDatasetUniquePtr ExportMosaickingVector(
      const std::string& output_path,
      const std::string& query_path,
      const std::string& query_raster_name_field_name,
      bool extension) override;

  GDALDatasetUniquePtr ExportBorderVector(
      const std::string& output_path) override;

  std::vector<std::string> ExportAllRastersPath() override;

 private:
  std::shared_ptr<MosaickingInterface> mosaicking_;
  std::shared_ptr<color_balancing::ColorBalancingInterface> color_balancing_;

  GDALDataset* composite_table_dataset_;
  GDALDataset* border_dataset_;
  OGRLayer* composite_table_layer_;
  OGRLayer* border_layer_;
  OGRGeometry* covered_border_;  // type is OGRMultiPolygon

  std::vector<double> borders_area_;
};

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_MOSAICKING_SERIAL_CONTAINER_H_