#include <iostream>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/pansharpening.h>
#include <rs-toolset/utils.hpp>


using namespace rs_toolset;

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "pansharpening", "A pansharpening tool using the Gram-Schmidt adaptive "
      "pansharpening algorithm to fuse the panchromatic(PAN) raster and the "
      "multispectral(MS) raster over the same area");
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
      ("log-level", "Availalbe levels are trace(t), debug(d), info(i), "
          "warn(w), err(e), critical(c)", cxxopts::value<std::string>());
  auto result(options.parse(argc, argv));
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return -1;
  }
  if (!result.count("pan") || !result.count("ms") || !result.count("output")) {
    std::cout << "Lack of the \"pan\", the \"ms\" or the \"output\" argument"
    << std::endl;
    return -1;
  }
  std::vector<int> pansharpen_bands_map;
  if (result.count("bands-map"))
    pansharpen_bands_map = result["bands-map"].as<std::vector<int>>();
  if (result.count("log-level")) {
    auto log_level(result["log-level"].as<std::string>());
    if (log_level == "t" || log_level == "trace") {
      utils::InitSpdlog("Pansharpening", spdlog::level::level_enum::trace);
    } else if (log_level == "d" || log_level == "debug") {
      utils::InitSpdlog("Pansharpening", spdlog::level::level_enum::debug);
    } else if (log_level == "i" || log_level == "info") {
      utils::InitSpdlog("Pansharpening", spdlog::level::level_enum::info);
    } else if (log_level == "w" || log_level == "warn") {
      utils::InitSpdlog("Pansharpening", spdlog::level::level_enum::warn);
    } else if (log_level == "e" || log_level == "err") {
      utils::InitSpdlog("Pansharpening", spdlog::level::level_enum::err);
    } else if (log_level == "c" || log_level == "critical") {
      utils::InitSpdlog("Pansharpening", spdlog::level::level_enum::critical);
    }
  } else {
    utils::InitSpdlog("Pansharpening");
  }

  utils::InitGdal(argv[0]);
  auto pansharpening(pansharpening::GramSchmidtAdaptive::Create(
      result["block-size"].as<int>()));
  pansharpening->Run(
      result["pan"].as<std::string>(), result["ms"].as<std::string>(), 
      result["output"].as<std::string>(), !result["disable-rpc"].as<bool>(),
      !result["disable-stretch"].as<bool>(), pansharpen_bands_map);
  if (result["build-overviews"].as<bool>()) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        result["output"].as<std::string>().c_str(),
        GDAL_OF_RASTER | GDAL_OF_READONLY));
    utils::CreateRasterPyra(dataset.get(), "DEFLATE");
  }
  return 0;
}