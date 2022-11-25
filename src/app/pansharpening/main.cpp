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
      "pansharpening", "A pansharpening tool using the available methods to "
      "create a pansharpened raster from the given input panchromatic(PAN) "
      "raster and the multispectral(MS) raster over the same area");
  options.add_options()
      ("h,help", "Print usage")
      ("method", "Available methods including \"GS\" and \"GSA\"", cxxopts::value<std::string>()->default_value("GSA"))
      ("pan", "The PAN raster path", cxxopts::value<std::string>())
      ("ms", "The MS raster path", cxxopts::value<std::string>())
      ("output", "The output raster path", cxxopts::value<std::string>())
      ("disable-rpc", "Disable using the RPC information if the PAN raster and the MS raster are both DOM", cxxopts::value<bool>()->default_value("false"))
      ("disable-stretching", "Disable Stretching the output raster", cxxopts::value<bool>()->default_value("false"))
      ("o,overviews", "Build overviews for the output raster", cxxopts::value<bool>()->default_value("false"))
      ("bands-map", "The bands' map from the MS raster to the output raster", cxxopts::value<std::vector<int>>())
      ("block-size", "The block size", cxxopts::value<int>()->default_value("4096"))
      ("log-level", "Available log levels including trace(t), debug(d), info(i), warn(w), err(e), critical(c), off(o)", cxxopts::value<std::string>()->default_value("i"));

  utils::InitGdal(argv[0]);
  auto result(options.parse(argc, argv));
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  if (!result.count("pan") || !result.count("ms") || !result.count("output")) {
    std::cout << "Lack of the \"pan\" argument, the \"ms\" argument or "
        "the \"output\" argument" << std::endl;
    return -1;
  }
  std::vector<int> bands_map;
  if (result.count("bands-map"))
    bands_map = result["bands-map"].as<std::vector<int>>();
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
    std::cout << "Available levels including trace(t), debug(d), info(i), "
        "warn(w), err(e), critical(c), off(o)" << std::endl;
    return -1;
  }

  std::shared_ptr<pansharpening::PansharpeningInterface> pansharpening(nullptr);
  if (auto method(result["method"].as<std::string>()); method == "GS") {
    pansharpening = pansharpening::GramSchmidt::Create(
        result["block-size"].as<int>());
  } else if (method == "GSA") {
    pansharpening = pansharpening::GramSchmidtAdaptive::Create(
        result["block-size"].as<int>());
  } else {
    std::cout << "Available methods include \"GS\" and \"GSA\"" << std::endl;
    return -1;
  }
  if (!pansharpening->Run(
          result["pan"].as<std::string>(), result["ms"].as<std::string>(), 
          result["output"].as<std::string>(), !result["disable-rpc"].as<bool>(),
          !result["disable-stretching"].as<bool>(), bands_map)) {
    return -1;
  }
  if (result["overviews"].as<bool>()) {
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        result["output"].as<std::string>().c_str(),
        GDAL_OF_RASTER | GDAL_OF_READONLY));
    utils::CreateRasterPyra(dataset.get());
  }
  return 0;
}