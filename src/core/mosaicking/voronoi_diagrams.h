#ifndef RS_TOOLSET_SRC_CORE_MOSAICKING_VORONOI_DIAGRAMS_H_
#define RS_TOOLSET_SRC_CORE_MOSAICKING_VORONOI_DIAGRAMS_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/mosaicking.h>

namespace rs_toolset {
namespace mosaicking {

class VoronoiDiagramsImpl final : public VoronoiDiagrams {
 public:
  explicit VoronoiDiagramsImpl(
      const std::shared_ptr<MosaickingInterface>& mosaicking)
      : mosaicking_(mosaicking) {
    spdlog::info(
        "Creating a voronoi diagrams with\n - With mosaicking: {}",
        bool(mosaicking));
  }
  VoronoiDiagramsImpl(const VoronoiDiagramsImpl&) = delete;
  VoronoiDiagramsImpl& operator=(const VoronoiDiagramsImpl&) = delete;
  ~VoronoiDiagramsImpl() = default;

  GDALDatasetUniquePtr Run(
      const std::vector<std::string>& rasters_path,
      const std::string& output_path,
      OGRSpatialReference* spatial_ref,
      bool with_refinement,
      int low_overview_trunc,
      int high_overview_trunc,
      const std::vector<int>& rgb_bands_map,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing) override;

 private:
  /**
   * @brief Refine seamlines in the initial mosaicking network
   * @param[in] paths The rasters' path
   * @param[in] borders The rasters' border
   * @param[in] neighbor_pairs The neighbor pairs in the mosaicking network layer
   * @param[in] low_overview_trunc The low overview trunction
   * @param[in] high_overview_trunc The high overview trunction
   * @param[in] rgb_bands_map The RGB bands' map
   * @param[in] color_balancing The color balancing shared pointer
   * @param[in,out] layer The mosaicking network layer
  */
  void RefineSeamlines(
      const std::vector<std::string>& paths,
      const std::vector<OGRGeometryUniquePtr>& borders,
      const std::vector<std::pair<int, int>>& neighbor_pairs,
      int low_overview_trunc,
      int high_overview_trunc,
      const std::vector<int>& rgb_bands_map,
      const std::shared_ptr<color_balancing::ColorBalancingInterface>&
          color_balancing,
      OGRLayer* layer);

  std::shared_ptr<MosaickingInterface> mosaicking_;
};

}  // namespace mosaicking
}  // namespace rs_toolset

#endif  // RS_TOOLSET_SRC_CORE_MOSAICKING_VORONOI_DIAGRAMS_H_