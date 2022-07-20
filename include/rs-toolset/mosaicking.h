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
   * @brief Run the mosaicking algorithm on the given raster with the former composite table and the covered border
   * @param[in] raster_path The given raster path
   * @param[in,out] composite_table_layer The former composite table layer
   * @param[in,out] composite_table_layer The border layer
   * @param[in,out] covered_border The former covered border
   * @param[in] last_overview_idx The last overview index, default is 3 and must be positive
   * @param[in] use_seamline Whether uses the seamline, default is false
   * @return Running state
   * @note This class should be used as the argument of the mosaicking container instead of running alone
  */
  virtual bool Run(
      const std::string& raster_path,
      OGRLayer* composite_table_layer,
      OGRLayer* border_layer,
      OGRGeometry* covered_border,
      int last_overview_idx = 3,
      bool use_seamline = false) = 0;
};

/** @brief Graph cut mosaicking class implementing the graph cut mosaicking algorithm */
class RS_TOOLSET_API GraphCut : virtual public MosaickingInterface {
 public:
  /**
   * @brief Create a graph cut mosaicking shared pointer
   * @param[in] grad_self_low The low trunction of the gradient-self term
   * @param[in] grad_self_high The high trunction of the gradient-self term
   * @param[in] grad_self_exp The exponential of the gradient-self term
   * @param[in] diff_low The low trunction of the difference term
   * @param[in] diff_exp The exponential of the difference term
   * @param[in] tol The tolerance in pixels for simplifying the seamline, default is 2.0. 
   * @return The output graph cut mosaicking shared pointer
  */
  static std::shared_ptr<GraphCut> Create(
      float grad_self_low,
      float grad_self_high,
      float grad_self_exp,
      float diff_low,
      float diff_exp,
      double tol = 2.0);
};

/** @brief Raster information class for sort */
struct RS_TOOLSET_API RasterInfo {
  std::string path;
  double reso;
  int date;
};

/**
 * @brief Sort two RasterInfo members by resolution
 * @param[in] info1 The first RasterInfo member
 * @param[in] info2 The second RasterInfo member
 * @return Sort result
*/
bool RS_TOOLSET_API SortByReso(
    const RasterInfo& info1,
    const RasterInfo& info2);

/** @brief Abstract mosaicking container class */
class RS_TOOLSET_API MosaickingContainerInterface {
 public:
  typedef std::function<bool(const RasterInfo&, const RasterInfo&)> SortFunc;

  /**
   * @brief Sort rasters by the given sort function
   * @param[in,out] rasters_path The rasters' path
   * @param[in] sort_func The given sort function, default is by resolution
  */
  virtual void SortRasters(
      std::vector<std::string>& rasters_path,
      const SortFunc& sort_func = rs_toolset::mosaicking::SortByReso) = 0;

  /**
   * @brief Add a raster task to execute the mosaicking algorithm
   * @param[in] raster_path The raster path
   * @param[in] last_overview_idx The last overview index, default is 3 and must be positive
   * @param[in] use_seamline Whether uses the seamline
   * @return Running state
  */
  virtual bool AddTask(
      const std::string& raster_path,
      int last_over_view_idx = 3,
      bool use_seamline = false) = 0;

  /**
   * @brief Export the internal composite table to the given path
   * @param[in] composit_table_path The exported composit table path
   * @param[in] buffer The buffer distance, only a negative number accepted
   * @param[in] tol The tolerance for simplifying the border
   * @param query_path[in] The query composite table path, default is empty
   * @param query_rasters_name_field_name[in] The query composite table field name representing the rasters' name, default is empty
   * @param with_extension[in] Whether with extension in the raster name field
   * @return Running state
   * @note The unit of "buffer" and "tol" arguments should be the same as the spatial reference's
  */
  virtual bool ExportCompositeTableVector(
      const std::string& composit_table_path,
      double buffer = 0.0,
      double tol = 0.0,
      const std::string& query_path = "",
      const std::string& query_rasters_name_field_name = "",
      bool with_extension = true) = 0;

  /**
   * @brief Export all rasters' name in the internal composite table
   * @return The output rasters' name
  */
  virtual std::vector<std::string> ExportAllRastersName() = 0;

  /**
   * @brief Create a mosaicking raster from the internal or the external composit table
   * @param[in] mosaicking_raster_path The output mosaicking raster path
   * @param[in] composit_table_path The external composit table path, use the internal composit table if empty
   * @param[in] rasters_dir The rasters' directory corresponding to the external composit table
   * @param[in] reso The output mosaicking raster resolution
   * @param[in] blend_dist The blend distance, default is 0.0
   * @return Running state
   * note The blend operation costs much more time than directly warp
  */
  virtual bool CreateMosaickingRaster(
      const std::string& mosaicking_raster_path,
      const std::string& composit_table_path,
      const std::string& rasters_dir,
      double reso,
      double blend_dist = 0.0) = 0;

  /** @todo virtual bool ExportStatusToJson() = 0; */

  /** @todo virtual bool ImportStatusFromStatus() = 0; */
};

/** @brief Mosaicking container class */
class RS_TOOLSET_API MosaickingContainer
    : virtual public MosaickingContainerInterface {
 public:
  /**
   * @brief Create a mosaicking container shared pointer
   * @param mosaicking[in] The mosaicking shared pointer
   * @param spatial_ref[in] The spatial reference, which needs to SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER)
   * @return The output mosaicking container shared pointer
  */
  static std::shared_ptr<MosaickingContainer> Create(
      std::shared_ptr<MosaickingInterface> mosaicking,
      OGRSpatialReference* spatial_ref);

  /**
   * @brief Create a mosaicking container shared pointer by the external composite tabel
   * @param mosaicking[in] The mosaicking shared pointer
   * @param composite_tabel_path[in] The external composite tabel path
   * @param rasters_dir[in] The rasters' directory corresponding to the external composit table
   * @return The output mosaicking container shared pointer
  */
  static std::shared_ptr<MosaickingContainer> Create(
      std::shared_ptr<MosaickingInterface> mosaicking,
      const std::string& external_composite_tabel_path,
      const std::string& rasters_dir);
};

} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_