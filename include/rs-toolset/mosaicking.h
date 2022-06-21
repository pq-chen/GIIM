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

class RS_TOOLSET_API MosaickingInterface {
 public:
  virtual bool Run(
      OGRLayer* composite_table_layer,
      OGRGeometry* covered_border,
      const std::string& raster_path,
      bool with_seamline = false) = 0;
};

class RS_TOOLSET_API GraphCut : virtual public MosaickingInterface {
 public:
  static std::shared_ptr<GraphCut> Create(
      double grad_exp,
      double min_diff,
      double max_diff);
};

struct RS_TOOLSET_API RasterInfo {
  std::string path;
  double reso;
  int date;
};

bool RS_TOOLSET_API SortByReso(
    const RasterInfo& info1,
    const RasterInfo& info2);

class RS_TOOLSET_API MosaickingContainerInterface {
 public:
  typedef std::function<bool(const RasterInfo&, const RasterInfo&)> SortFunc;

  virtual bool SortRasters(
      std::vector<std::string>& rasters_path,
      const SortFunc& sort_func = rs_toolset::mosaicking::SortByReso) = 0;

  virtual bool AddTask(
      const std::string& raster_path,
      bool with_seamline = false) = 0;

  virtual bool ExportCompositeTableVector(
      const std::string& composit_table_path,
      double reso,
      double buffer = -1.,
      double tol = 1.5) = 0;

  virtual bool CreateMosaickingRaster(
      const std::string& mosaicking_raster_path,
      const std::string& rasters_dir,
      double reso) = 0;

  // TODO: virtual bool ExportStatusToJson() = 0;

  // TODO: virtual bool ImportStatusFromStatus() = 0;

  static bool SortByReso(
      const RasterInfo& info1,
      const RasterInfo& info2);
};

class RS_TOOLSET_API MosaickingContainer 
    : virtual public MosaickingContainerInterface {
 public:
  static std::shared_ptr<MosaickingContainer> Create(
      std::shared_ptr<MosaickingInterface> mosaicking,
      OGRSpatialReference* spatial_ref);
  static std::shared_ptr<MosaickingContainer> Create(
      std::shared_ptr<MosaickingInterface> mosaicking,
      const std::string& seamline_vector_path,
      const std::string& rasters_dir);
};


} // mosaicking
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_MOSAICKING_H_

