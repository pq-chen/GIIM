#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else  // RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif  // RS_TOOLSET_EXPORTS

#else  // _WIN32
#define RS_TOOLSET_API
#endif  // _WIN32

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <rs-toolset/color_balancing.h>

namespace rs_toolset {
namespace mosaicking {

/** @brief Raster information class */
struct RS_TOOLSET_API RasterInfo {
  std::string path;
  double reso;
  int date;
};

/**
 * @brief Sort two RasterInfo instances by resolution
 * @param[in] info1 The first RasterInfo instance
 * @param[in] info2 The second RasterInfo instance
 * @return Sort result
*/
bool RS_TOOLSET_API SortByReso(
    const RasterInfo& info1,
    const RasterInfo& info2);

/** @brief Abstract mosaicking class for running the mosaicking algorithm */
class RS_TOOLSET_API MosaickingInterface {
 public:
  /**
   * @brief Run the mosaicking algorithm between the given raster and the former composite table
   * @param[in] path The given raster path
   * @param[in,out] composite_table_layer The former composite table layer
   * @param[in,out] border_layer The former border layer
   * @param[in,out] covered_border The former covered border, type it OGRMultiPolygon
   * @param[in,out] borders_area The former borders' area
   * @param[in] rejection_ratio The regularization term to reject a new item in the output composite table, default is 0.01
   * @param[in] low_overview_trunc The low overview trunction, default is 3 and must be non-negative
   * @param[in] high_overview_trunc The high overview trunction, default is 1 and must be non-negative
   * @param[in] rgb_bands_map The RGB bands' map, default is empty means first three bands
   * @param[in] color_balancing The color balancing shared pointer, default is nullptr
   * @return Running state
   * @note This method should be used in the mosaicking container instead of running alone
  */
  virtual bool RunTaskForExisting(
      const std::string& path,
      OGRLayer* composite_table_layer,
      OGRLayer* border_layer,
      OGRGeometry* covered_border,
      std::vector<double>& borders_area,
      double rejection_ratio = 0.01,
      int low_overview_trunc = 3,
      int high_overview_trunc = 1,
      const std::vector<int>& rgb_bands_map = {},
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing = nullptr) = 0;

  /**
   * @brief Run the mosaicking algorithm on the given raster pair with the corresponding covered geometries
   * @param[in] path1 The first raster path
   * @param[in] path2 The second raster path
   * @param[in] border1 The first border, create a new border if empty
   * @param[in] border2 The second border, create a new border if empty
   * @param[in] geometry1 The first covered geometry
   * @param[in] geometry2 The second covered geometry
   * @param[in] low_overview_trunc The low overview trunction, default is 3 and must be non-negative
   * @param[in] high_overview_trunc The high overview trunction, default is 1 and must be non-negative
   * @param[in] rgb_bands_map The RGB bands' map, default is empty means first three bands
   * @param[in] color_balancing The color balancing shared pointer, default is nullptr
   * @param[in] color_balancing_idx1 The first color balancing index, default is -1 means no color balancing
   * @param[in] color_balancing_idx2 The second color balancing index, default is -1 means no color balancing
   * @return The output seamline geometry, type is OGRLineString
   * @note The border1 and the border2 should have the same spatial reference
  */
  virtual OGRGeometryUniquePtr RunTaskForPair(
      const std::string& path1,
      const std::string& path2,
      OGRGeometry* border1,
      OGRGeometry* border2,
      OGRGeometry* geometry1,
      OGRGeometry* geometry2,
      int low_overview_trunc = 3,
      int high_overview_trunc = 1,
      const std::vector<int>& rgb_bands_map = {},
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing = nullptr,
      int color_balancing_idx1 = -1,
      int color_balancing_idx2 = -1) = 0;
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
   * @param[in] tol The tolerance in pixels for simplifying the seamline, default is 2.0
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

/** @brief Mosaicking container class */
class RS_TOOLSET_API MosaickingContainer {
 public:
  using SortFunc = std::function<bool(const RasterInfo&, const RasterInfo&)>;

  /**
   * @brief Create a mosaicking container shared pointer
   * @param[in] mosaicking The mosaicking shared pointer
   * @param[in] spatial_ref The spatial reference, which needs SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER) since GDAL 3.0
   * @param[in] rejection_ratio The regularization term to reject a new item in the output composite table, default is 0.01
   * @param[in] color_balancing The color balancing shared pointer, default is nullptr
   * @return The output mosaicking container shared pointer
  */
  static std::shared_ptr<MosaickingContainer> Create(
      const std::shared_ptr<MosaickingInterface>& mosaicking,
      OGRSpatialReference* spatial_ref,
      double rejection_ratio = 0.01,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing = nullptr);

  /**
   * @brief Create a mosaicking container shared pointer initialized by the external composite tabel
   * @param[in] mosaicking The mosaicking shared pointer
   * @param[in] composite_tabel_path The external composite tabel path
   * @param[in] rasters_dir The rasters' directory corresponding to the external composite table
   * @param[in] rejection_ratio The regularization term to reject a new item in the output composite table, default is 0.01
   * @param[in] color_balancing The color balancing shared pointer, default is nullptr
   * @return The output mosaicking container shared pointer
  */
  static std::shared_ptr<MosaickingContainer> Create(
      const std::shared_ptr<MosaickingInterface>& mosaicking,
      const std::string& composite_tabel_path,
      const std::string& rasters_dir,
      double rejection_ratio = 0.01,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing = nullptr);

  /**
   * @brief Sort rasters by the given sort function
   * @param[in,out] paths The rasters' path
   * @param[in] sort_func The given sort function, default is SortByReso
   * @return Running state
  */
  virtual bool SortRasters(
      std::vector<std::string>& paths,
      const SortFunc& sort_func = rs_toolset::mosaicking::SortByReso) = 0;

  /**
   * @brief Add a task to run the mosaicking algorithm with the internal composite table
   * @param[in] path The task raster path
   * @param[in] low_overview_trunc The low overview trunction, default is 3 and must be non-negative
   * @param[in] high_overview_trunc The high overview trunction, default is 1 and must be non-negative
   * @param[in] rgb_bands_map The RGB bands' map, default is empty means first three bands
   * @return Running state
  */
  virtual bool AddTask(
      const std::string& path,
      int low_overview_trunc = 3,
      int high_overview_trunc = 1,
      const std::vector<int>& rgb_bands_map = {}) = 0;

  /**
   * @brief Export the internal composite table to the given path
   * @param[in] output_path The exported composite table path
   * @param[in] query_path The query composite table path, default is empty
   * @param[in] query_rasters_name_field_name The query composite table field name representing the raster name, default is empty
   * @param[in] with_extension Whether remains the extension in the raster name field or not, default is true
   * @return The output composite table
  */
  virtual GDALDatasetUniquePtr ExportCompositeTable(
      const std::string& output_path,
      const std::string& query_path = "",
      const std::string& query_rasters_name_field_name = "",
      bool with_extension = true) = 0;

  /**
   * @brief Export all rasters' name in the internal composite table
   * @return The output rasters' name
  */
  virtual std::vector<std::string> ExportAllRastersName() = 0;

  /** @todo virtual bool ExportStatusToJson() = 0; */

  /** @todo virtual bool ImportStatusFromStatus() = 0; */
};

/** @brief Voronoi diagrams class */
class RS_TOOLSET_API VoronoiDiagrams {
 public:
  /**
   * @brief Create a voronoi diagrams shared pointer
   * @param[in] mosaicking The mosaicking shared pointer, default is nullptr
   * @return The output voronoi diagrams shared pointer
  */
  static std::shared_ptr<VoronoiDiagrams> Create(
      const std::shared_ptr<MosaickingInterface>& mosaicking = nullptr);

  /**
   * @brief Run the voronoi diagrams algorithm on the given rasters
   * @param[in] rasters_path The given rasters' path
   * @param[in] output_path The output composite table path
   * @param[in] spatial_ref The spatial reference, which needs SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER) since GDAL 3.0
   * @param[in] with_refinement Whether needs to refine the seamlines in the voronoi diagrams or not, default is true
   * @param[in] low_overview_trunc The low overview trunction, default is 3 and must be non-negative
   * @param[in] high_overview_trunc The high overview trunction, default is 1 and must be non-negative
   * @param[in] rgb_bands_map The RGB bands' map, default is empty means first three bands
   * @param[in] color_balancing The color balancing shared pointer, default is nullptr
   * @return The output composite table
  */
  virtual GDALDatasetUniquePtr Run(
      const std::vector<std::string>& rasters_path,
      const std::string& output_path,
      OGRSpatialReference* spatial_ref,
      bool with_refinement = true,
      int low_overview_trunc = 3,
      int high_overview_trunc = 1,
      const std::vector<int>& rgb_bands_map = {},
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing = nullptr) = 0;
};

/**
  * @brief Create a mosaicking raster from the given composite table
  * @param[in] path The output mosaicking raster path
  * @param[in] layer The given composite table layer
  * @param[in] reso The output resolution
  * @param[in] color_balancing The color balancing shared pointer, default is nullptr
  * @param[in] blend_dist The output blend distance in pixels, default is 0.0
  * @return The output mosaicking raster dataset
  * @note The blend operation costs much more time than directly warping
*/
RS_TOOLSET_API GDALDatasetUniquePtr CreateMosaickingRaster(
    const std::string& path,
    OGRLayer* layer,
    double reso,
    const std::shared_ptr<color_balancing::ColorBalancingInterface>&
        color_balancing = nullptr,
    double blend_dist = 0.0);

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_