#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <opencv2/opencv.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/stretching.h>
#include <rs-toolset/utils.hpp>

namespace fs = std::filesystem;
using namespace rs_toolset;

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "stretching", "A stretching tool using the available methods to create "
      "stretched unsigned 8-bit rasters from the given input unsigned 8-bit or "
      "16-bit rasters");
  options.add_options("basic")
      ("h,help", "Print usage")
      ("method", "Available methods include \"percent-clip\" and \"standard-deviations\"", cxxopts::value<std::string>()->default_value("percent-clip"))
      ("input", "The input rasters' path or directories, only supports rasters with unsigned 8-bit or 16-bit ", cxxopts::value<std::vector<std::string>>())
      ("output-dir", "The output rasters' directory", cxxopts::value<std::string>())
      ("output-suffix", "The output rasters' suffix", cxxopts::value<std::string>()->default_value(""))
      ("o,overviews", "Build overviews for the output rasters", cxxopts::value<bool>()->default_value("false"))
      ("log-level", "Available levels are trace(t), debug(d), info(i), warn(w), err(e), critical(c), off(o)", cxxopts::value<std::string>()->default_value("i"));
  options.add_options("percent-clip-feature")
      ("low-percent", "The low percent trunction", cxxopts::value<double>()->default_value("0.005"))
      ("high-percent", "The high percent trunction", cxxopts::value<double>()->default_value("0.005"));
  options.add_options("standard-deviations-feature")
      ("scale", "The trunction scale", cxxopts::value<double>()->default_value("2.5"));
  auto result(options.parse(argc, argv));

  utils::InitGdal(argv[0]);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  if (!result.count("input") || !result.count("output-dir")) {
    std::cout << "Lack of the \"input\" argument or the \"output-dir\" argument"
        << std::endl;
    return -1;
  }
  auto log_level(result["log-level"].as<std::string>());
  if (log_level == "t" || log_level == "trace") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::trace);
  } else if (log_level == "d" || log_level == "debug") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::debug);
  } else if (log_level == "i" || log_level == "info") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::info);
  } else if (log_level == "w" || log_level == "warn") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::warn);
  } else if (log_level == "e" || log_level == "err") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::err);
  } else if (log_level == "c" || log_level == "critical") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::critical);
  } else if (log_level == "o" || log_level == "off") {
    utils::InitSpdlog("Mosaicking", spdlog::level::level_enum::off);
  } else {
    std::cout << "Available levels are trace(t), debug(d), info(i), warn(w), "
        "err(e), critical(c), off(o)" << std::endl;
    return -1;
  }

  std::shared_ptr<stretching::StretchingInterface> stretching(nullptr);
  if (auto method(result["method"].as<std::string>());
      method == "percent-clip") {
    stretching = stretching::PercentClip::Create(
        result["low-percent"].as<double>(),
        result["high-percent"].as<double>());
  } else if (method == "standard-deviations") {
    stretching = stretching::StandardDeviations::Create(
        result["scale"].as<double>());
  } else {
    std::cout << "Available methods include \"percent-clip\" and "
        "\"standard-deviations\"" << std::endl;
    return -1;
  }
  std::vector<std::string> paths;
  for (const auto& path : result["input"].as<std::vector<std::string>>()) {
    if (!fs::exists(path)) {
      std::cout << "Input path " << path << " does not exist" << std::endl;
      return -1;
    }
    if (fs::is_directory(path)) {
      for (const auto& entry : fs::directory_iterator(path)) {
        if (auto subpath(entry.path().string());
            utils::GetRasterDriverByPath(subpath)) {
          paths.push_back(subpath);
        }
      }
    } else if (utils::GetRasterDriverByPath(path)) {
      paths.push_back(path);
    }
  }
  std::cout << "The following raster(s) will be operated:" << std::endl;
  for (const auto& path : paths)
    std::cout << path << std::endl;
  std::cout << paths.size() << " task(s) in total" << std::endl;
  for (int i(0); i < paths.size(); ++i) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        paths[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
    if (!dataset) {
      spdlog::error("Opening {} failed", paths[i]);
      return -1;
    }

    double geotrans[6];
    dataset->GetGeoTransform(geotrans);
    auto mat(utils::CreateMatFromDataset(dataset.get()));
    if (!stretching->AccumulateStat(mat))
      return -1;
    stretching->CreateLutMats();
    stretching->Run(mat);
    std::string output_path(
        result["output-dir"].as<std::string>() + "/" +
            fs::path(paths[i]).stem().string() +
            result["output-suffix"].as<std::string>() +
            fs::path(paths[i]).extension().string());
    utils::CreateDatasetFromMat(
        mat, output_path, geotrans,
        const_cast<OGRSpatialReference*>(dataset->GetSpatialRef()), 0, nullptr);
    stretching->Clear();
    if (result["overviews"].as<bool>()) {
      GDALDatasetUniquePtr output_dataset(GDALDataset::Open(
          output_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
      utils::CreateRasterPyra(output_dataset.get());
    }
    spdlog::info("---------- {}/{} - done ----------", i + 1, paths.size());
  }
  return 0;
}