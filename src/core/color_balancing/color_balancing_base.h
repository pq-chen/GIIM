#ifndef RS_TOOLSET_SRC_CORE_COLOR_BALANCING_COLOR_BALANCING_BASE_H_
#define RS_TOOLSET_SRC_CORE_COLOR_BALANCING_COLOR_BALANCING_BASE_H_

#include <string>
#include <vector>

#include <gdal_priv.h>
#include <gdalwarper.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>


namespace rs_toolset {
namespace color_balancing {

class ColorBalancingBase : virtual public ColorBalancingInterface {
public:
  ColorBalancingBase(int block_size) : block_size_(block_size) {
    spdlog::info(
        "Creating a color balancing base with\n - Block size: {}", block_size);
  }
  ColorBalancingBase(const ColorBalancingBase&) = delete;
  ColorBalancingBase& operator=(const ColorBalancingBase&) = delete;
  virtual ~ColorBalancingBase() = default;

  std::vector<std::string> ExportAllRastersName() override;

  bool CreateRasters(
      const std::vector<int>& idxes,
      const std::string& output_dir) override;

  bool WarpByGeometry(
      const std::vector<int>& idxes,
      const std::vector<OGRGeometry*>& geometries,
      GDALDatasetUniquePtr& output_dataset,
      const std::vector<int>& bands_map,
      GDALResampleAlg resample_arg,
      double blend_dist,
      double nodata_value) override;

protected:
  /**
   * @brief Write the given index raster color balancing result from the source dataset to the output dataset
   * @param[in] idx The given index
   * @param[in] source_dataset The source dataset
   * @param[in,out] output_dataset The output dataset
  */
  virtual void WriteToDataset(
      int idx,
      GDALDataset* source_dataset,
      GDALDatasetUniquePtr& output_dataset) = 0;

  int block_size_;

  std::vector<std::string> rasters_path_;
};

} // namespace color_balancing
} // namespace rs_toolset

#endif // RS_TOOLSET_SRC_CORE_COLOR_BALANCING_COLOR_BALANCING_BASE_H_