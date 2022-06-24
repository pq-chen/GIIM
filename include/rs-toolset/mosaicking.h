#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else // LCMAKE_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif // LCMAKE_EXPORTS

#else // _WIN32
#define RS_TOOLSET_API
#endif // _WIN32

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <ogrsf_frmts.h>


namespace rs_toolset {
namespace mosaicking {

/** @brief Abstract mosaicking class for running the mosaicking algorithm in the mosaicking container */
class RS_TOOLSET_API MosaickingInterface {
 public:
  /**
   * @brief Run the mosaicking algorithm on the given raster with the composite table and the covered border
   * @param raster_path[in]: The given raster path
   * @param[in,out] composite_table_layer: The composite table layer
   * @param covered_border[in,out]: The covered border
   * @param[in] use_seamline: Whether uses the seamline, default is false
   * @return Running state
   * @warning This class should be used as the argument of the mosaicking container
  */
  virtual bool Run(
      const std::string& raster_path,
      OGRLayer* composite_table_layer,
      OGRGeometry* covered_border,
      bool use_seamline = false) = 0;
};

/** @brief Graph cut mosaicking class implementing the graph cut algorithm */
class RS_TOOLSET_API GraphCut : virtual public MosaickingInterface {
 public:
   /**
    * @brief Create the graph cut mosaicking shared pointer
    * @param[in] grad_term_exp: The gradient term exponential
    * @param[in] diff_term_low_trunc: The difference term low trunction
    * @param[in] diff_term_high_trunc: The difference term high trunction
    * @return The output graph cut mosaicking shared pointer
   */
   static std::shared_ptr<GraphCut> Create(
      double grad_term_exp,
      double diff_term_low_trunc,
      double diff_term_high_trunc);
};

/**
 * @brief Raster information class for sort
*/
struct RS_TOOLSET_API RasterInfo {
  std::string path;
  double reso;
  int date;
};

/**
 * @brief Sort two RasterInfo members by resolution
 * @param info1: The first RasterInfo member
 * @param info2: The second RasterInfo member
 * @return Sort result
*/
bool RS_TOOLSET_API SortByReso(
    const RasterInfo& info1,
    const RasterInfo& info2);

/**
 * @brief Abstract mosaicking container class
*/
class RS_TOOLSET_API MosaickingContainerInterface {
 public:
  typedef std::function<bool(const RasterInfo&, const RasterInfo&)> SortFunc;

  /**
   * @brief Sort the given rasters' path with the sort function
   * @param rasters_path: The given rasters' path
   * @param sort_func: The sort function, default is by resolution
  */
  virtual void SortRasters(
      std::vector<std::string>& rasters_path,
      const SortFunc& sort_func = rs_toolset::mosaicking::SortByReso) = 0;

  /**
   * @brief Add the given raster path task to execute mosaicking
   * @param raster_path: The given raster path
   * @param use_seamline: Whether uses the seamline
   * @return Running state
  */
  virtual bool AddTask(
      const std::string& raster_path,
      bool use_seamline = false) = 0;

  /**
   * @brief Export the composite table vector in the mosaicking container to the given path
   * @param composit_table_path: The given output composit table path
   * @param unit: The unit value for scaling buffer and tolerance
   * @param buffer: The buffer distance
   * @param tol: The tolerance for simplification
   * @return Running state
  */
  virtual bool ExportCompositeTableVector(
      const std::string& composit_table_path,
      double unit,
      double buffer = -1.,
      double tol = 1.5) = 0;

  /**
   * @brief Export all rasters' name in the composite table vector
   * @return The output rasters' name
  */
  virtual std::vector<std::string> ExportAllRastersName() = 0;

  /**
   * @brief Create the mosaicking raster from the given composit table
   * @param mosaicking_raster_path: The otuput mosaicking raster path
   * @param composit_table_path: The given composit table path, use the internal composit table if empty
   * @param rasters_dir: The raster directory corresponding to the composit table
   * @param reso: The output mosaicking raster resolution
   * @return Running state
  */
  virtual bool CreateMosaickingRaster(
      const std::string& mosaicking_raster_path,
      const std::string& composit_table_path,
      const std::string& rasters_dir,
      double reso) = 0;

  // TODO: virtual bool ExportStatusToJson() = 0;

  // TODO: virtual bool ImportStatusFromStatus() = 0;
};

/**
 * @brief Mosaicking container class
*/
class RS_TOOLSET_API MosaickingContainer 
    : virtual public MosaickingContainerInterface {
 public:
  /**
   * @brief Create the mosaicking container shared pointer
   * @param mosaicking: The mosaicking
   * @param spatial_ref: The spatial reference, which needs to SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER)
   * @return The output mosaicking container shared pointer
  */
  static std::shared_ptr<MosaickingContainer> Create(
      std::shared_ptr<MosaickingInterface> mosaicking,
      OGRSpatialReference* spatial_ref);

  /**
   * @brief Create the mosaicking container shared pointer with the given composite tabel
   * @param mosaicking: The mosaicking
   * @param composite_tabel_path: The given composite tabel path
   * @param rasters_dir: The raster directory corresponding to the composit table
   * @return The output mosaicking container shared pointer
  */
  static std::shared_ptr<MosaickingContainer> Create(
      std::shared_ptr<MosaickingInterface> mosaicking,
      const std::string& composite_tabel_path,
      const std::string& rasters_dir);
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_

