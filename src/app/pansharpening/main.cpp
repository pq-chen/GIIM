#include <iostream>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/pansharpening.h>
#include <rs-toolset/utils.hpp>


int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "Pansharpening", 
      "A pansharpening tool using the Gram-Schmidt adaptive algorithm to fuse "
      "the panchromatic(PAN) raster and the multispectral(MS) raster over the "
      "same area.");
  options.add_options()
      ("h,help", "Print usage")
      ("pan", "The PAN raster path", cxxopts::value<std::string>())
      ("ms", "The MS raster path", cxxopts::value<std::string>())
      ("output", "The output pansharpened raster path", 
          cxxopts::value<std::string>())
      ("disable-rpc", "Disable using the RPC information "
          "if the PAN raster and the MS raster are both DOM", 
          cxxopts::value<bool>()->default_value("false"))
      ("disable-stretch", "Disable Stretching the pansharpened raster",
          cxxopts::value<bool>()->default_value("false"))
      ("o,build-overviews", "Build overviews for the pansharpened raster",
          cxxopts::value<bool>()->default_value("true"))
      ("bands-map", "The Bands' map from the MS raster to the output "
          "pansharpened raster", cxxopts::value<std::vector<int>>())
      ("block-size", "The block size per operation",
          cxxopts::value<int>()->default_value("16384"))
      ("log-level", "Availalbe levels are trace(t), debug(d), info(i), warn(w),"
          " err(e), critical(c)", cxxopts::value<std::string>());
  auto result = options.parse(argc, argv);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  std::string pan_path, ms_path, pansharpened_path;
  if (result.count("pan") && result.count("ms") && result.count("output")) {
    pan_path = result["pan"].as<std::string>();
    ms_path = result["ms"].as<std::string>();
    pansharpened_path = result["output"].as<std::string>();
  } else {
    std::cout << options.help() << std::endl;
    return 0;
  }
  bool use_rpc(!result["disable-rpc"].as<bool>()), 
      use_stretch(!result["disable-stretch"].as<bool>()),
      build_overviews(result["build-overviews"].as<bool>());
  std::vector<int> pansharpen_bands_map;
  if (result.count("bands-map")) {
    pansharpen_bands_map = result["bands-map"].as<std::vector<int>>();
  }
  int block_size(result["block-size"].as<int>());
  if (result.count("log-level")) {
    std::string log_level(result["log-level"].as<std::string>());
    if (log_level == "t" || log_level == "trace") {
      rs_toolset::utils::InitSpdlog(
          "Pansharpening", spdlog::level::level_enum::trace);
    } else if (log_level == "d" || log_level == "debug") {
      rs_toolset::utils::InitSpdlog(
          "Pansharpening", spdlog::level::level_enum::debug);
    } else if (log_level == "i" || log_level == "info") {
      rs_toolset::utils::InitSpdlog(
          "Pansharpening", spdlog::level::level_enum::info);
    } else if (log_level == "w" || log_level == "warn") {
      rs_toolset::utils::InitSpdlog(
          "Pansharpening", spdlog::level::level_enum::warn);
    } else if (log_level == "e" || log_level == "err") {
      rs_toolset::utils::InitSpdlog(
          "Pansharpening", spdlog::level::level_enum::err);
    } else if (log_level == "c" || log_level == "critical") {
      rs_toolset::utils::InitSpdlog(
          "Pansharpening", spdlog::level::level_enum::critical);
    }
  } else {
    rs_toolset::utils::InitSpdlog("Pansharpening");
  }

  rs_toolset::utils::InitGdal(argv[0]);
  auto pansharpening(rs_toolset::pansharpening::GramSchmidtAdaptive::Create(
      block_size));
  pansharpening->Run(
      pan_path, ms_path, pansharpened_path, use_rpc, use_stretch,
      pansharpen_bands_map);
  GDALDatasetUniquePtr dataset(GDALDataset::Open(
      pansharpened_path.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY));
  if (build_overviews)
    rs_toolset::utils::CreateRasterPyra(dataset.get(), "DEFLATE");
  return 0;
}