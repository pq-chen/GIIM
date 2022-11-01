#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_CONTAINER_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_CONTAINER_H_

#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/mosaicking.h>

namespace rs_toolset {
namespace mosaicking {

class MosaickingContainerImpl final : public MosaickingContainer {
 public:
  explicit MosaickingContainerImpl(
      const std::shared_ptr<MosaickingInterface>& mosaicking,
      OGRSpatialReference* spatial_ref,
      double rejection_ratio,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing);
  MosaickingContainerImpl(const MosaickingContainerImpl&) = delete;
  MosaickingContainerImpl& operator=(const MosaickingContainerImpl&) = delete;
  virtual ~MosaickingContainerImpl();

  /**
   * @brief Initialize the mosaicking container by the external composite table
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
      int low_overview_trunc,
      int high_overview_trunc,
      const std::vector<int>& rgb_bands_map) override;

  GDALDatasetUniquePtr ExportCompositeTable(
      const std::string& output_path,
      const std::string& query_path,
      const std::string& query_rasters_name_field_name,
      bool with_extension) override;

  std::vector<std::string> ExportAllRastersName() override;

 private:
  std::shared_ptr<MosaickingInterface> mosaicking_;
  std::shared_ptr<color_balancing::ColorBalancingInterface> color_balancing_;

  GDALDataset* composite_table_dataset_;
  GDALDataset* border_dataset_;
  OGRLayer* composite_table_layer_;
  OGRLayer* border_layer_;
  OGRGeometry* covered_border_;  // type is OGRMultiPolygon

  std::vector<double> borders_area_;
  double rejection_ratio_;
};

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_MOSAICKING_MOSAICKING_CONTAINER_H_