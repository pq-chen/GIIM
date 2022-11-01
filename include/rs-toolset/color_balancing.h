#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_COLOR_BALANCING_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_COLOR_BALANCING_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else // RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif // RS_TOOLSET_EXPORTS

#else // _WIN32
#define RS_TOOLSET_API
#endif // _WIN32

#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <gdalwarper.h>
#include <ogrsf_frmts.h>


namespace rs_toolset {
namespace color_balancing {

/** @brief Abstract color balancing class */
class RS_TOOLSET_API ColorBalancingInterface {
 public:
  /**
   * @brief Create a given method color balancing shared pointer
   * @param[in] method The given method
   * @param[in] block_size The block size per operation, default is 16384
   * @return The output color balancing shared pointer, nullptr if failed
  */
  static std::shared_ptr<ColorBalancingInterface> Create(
      const std::string& method,
      int block_size = 16384);

  /** 
   * @brief Create a color balancing shared pointer from the given JSON path
   * @param[in] json_path The given JSON path
   * @param[in] block_size The block size per operation, default is 16384
   * @return The output color balancing shared pointer, nullptr if failed
  */
  static std::shared_ptr<ColorBalancingInterface> CreateFromArgusJson(
      const std::string& json_path,
      int block_size = 16384);

  /**
   * @brief Calculate color balancing arguments for the given rasters' path
   * @param[in] rasters_path The given rasters' path
   * @return Running state
  */
  virtual bool CalcArgus(const std::vector<std::string>& rasters_path) = 0;

  /**
   * @brief Export the existed arguments to the given JSON path
   * @param[in] json_path The given JSON path
   * @return Running state
  */
  virtual bool ExportArgusJson(const std::string& json_path) = 0;

  /**
   * @brief Import arguments from the given JSON path
   * @param[in] json_path The given JSON path
   * @return Running state
  */
  virtual bool ImportArgusJson(const std::string& json_path) = 0;

  /**
   * @brief Export all rasters' name in the color balancing class
   * @return The output rasters' name
  */
  virtual std::vector<std::string> ExportAllRastersName() = 0;

  /**
   * @brief Create color balanced rasters with the internal argumets to the given output directory
   * @param[in] idxes The rasters' index corresponding to the internal rasters' path, empty means all rasters
   * @param[in] output_dir The given output directory
   * @return Running state
  */
  virtual bool CreateRasters(
      const std::vector<int>& idxes,
      const std::string& output_dir) = 0;

  /**
   * @brief Warp source datasets to the output dataset by the given geometries as cutlines with the internal argumets
   * @param[in] idxes The rasters' index corresponding to the internal rasters' path
   * @param[in] geometries The given geometries
   * @param[out] output_dataset The output dataset
   * @param[in] bands_map The bands map, default is empty means all bands
   * @param[in] resample_arg The resample argument, default is GRA_Bilinear
   * @param[in] blend_dist The blend distance, default is 0.0
   * @param[in] nodata_value The nodata value, default is 0.0
   * @return Running state
  */
  virtual bool WarpByGeometry(
      const std::vector<int>& idxes,
      const std::vector<OGRGeometry*>& geometries,
      GDALDatasetUniquePtr& output_dataset,
      const std::vector<int>& bands_map = {},
      GDALResampleAlg resample_arg = GRA_Bilinear,
      double blend_dist = 0.0,
      double nodata_value = 0.0) = 0;
};

/** @brief Global-to-local color balancing class implementing the global-to-local color balancing algorithm */
class RS_TOOLSET_API GlobalToLocal : virtual public ColorBalancingInterface {
 public:
  /**
   * @brief Create a global-to-local color balancing shared pointer
   * @param[in] block_size The block size per operation, default is 16384
   * @return The output global-to-local color balancing shared pointer
  */
  static std::shared_ptr<GlobalToLocal> Create(int block_size = 16384);
};

} // namespace color_balancing
} // namespace rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_COLOR_BALANCING_H_