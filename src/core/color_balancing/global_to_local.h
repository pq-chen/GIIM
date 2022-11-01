#ifndef RS_TOOLSET_SRC_CORE_COLOR_BALANCING_GLOBAL_TO_LOCAL_H_
#define RS_TOOLSET_SRC_CORE_COLOR_BALANCING_GLOBAL_TO_LOCAL_H_

#pragma warning(disable:4250)

#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include "color_balancing_base.h"
#include <rs-toolset/color_balancing.h>
#include <rs-toolset/utils.hpp>


namespace rs_toolset {
namespace color_balancing {

class GlobalToLocalImpl final
    : public ColorBalancingBase, public GlobalToLocal {
 public:
  GlobalToLocalImpl(int block_size) : ColorBalancingBase(block_size) {
    spdlog::info("Creating a global-to-local color balancing");
  }
  GlobalToLocalImpl(const GlobalToLocalImpl&) = delete;
  GlobalToLocalImpl& operator=(const GlobalToLocalImpl&) = delete;
  ~GlobalToLocalImpl() = default;

  bool CalcArgus(const std::vector<std::string>& rasters_path) override;

  bool ExportArgusJson(const std::string& json_path) override;

  bool ImportArgusJson(const std::string& json_path) override;

 private:
  struct Statistic {
    int pixels_count_ = 0;
    std::vector<double> sums_{ 0.0, 0.0, 0.0 };
    std::vector<double> square_sums_{ 0.0, 0.0, 0.0 };
    std::vector<float> means{ 0.0, 0.0, 0.0 };
    std::vector<float> stddevs{ 0.0, 0.0, 0.0 };
  };
  struct OverlapInfo {
    int idx1;
    int idx2;
    Statistic s1;
    Statistic s2;
    float weight;
  };

  std::vector<cv::Mat> CreateRasterData(GDALDataset* dataset, int* range) {
    std::vector<cv::Mat> mats;
    cv::split(utils::CreateMatFromDataset(dataset, range), mats);
    return mats;
  }

  std::vector<cv::Mat> CreateOverlapData(
      GDALDataset* dataset,
      int* range,
      OGRGeometry* geometry);

  void UpdateStatistic(const std::vector<cv::Mat>& mats, Statistic& s);

  void CalcStatistic(Statistic& s) {
    for (int b = 0; b < 3; b++) {
      s.means[b] = static_cast<float>(s.sums_[b] / s.pixels_count_);
      s.stddevs[b] = static_cast<float>(sqrt(
          s.square_sums_[b] / s.pixels_count_ - s.means[b] * s.means[b]));
    }
  }

  void CalcRastersGlobalArgus(
      const std::vector<Statistic>& rasters_statistic,
      const std::vector<OverlapInfo>& overlaps_info);

  void UpdateGridInfo(
      int idx,
      const std::vector<cv::Mat>& mats,
      int* raster_range,
      int* grid_range,
      double* composed_geotrans,
      cv::Mat& global_grid_pixels_count_mat);

  void WriteToDataset(
      int idx,
      GDALDataset* source_dataset,
      GDALDatasetUniquePtr& output_dataset) override;

  void ExecuteColorBalancing(
      int idx,
      int* range,
      double* composed_geotrans,
      std::vector<cv::Mat>& mats);

  int grid_size_ = 100;

  std::vector<std::vector<float>> rasters_global_argus_;
  std::vector<cv::Mat> rasters_grid_mean_mat_;
  cv::Mat global_grid_mean_mat_;
  float max_mean_;
  double inv_grid_geotrans_[6];
};

} // namespace color_balancing
} // namespace rs_toolset

#endif // RS_TOOLSET_SRC_CORE_COLOR_BALANCING_GLOBAL_TO_LOCAL_H_