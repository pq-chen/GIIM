#include <filesystem>
#include <regex>
#include <string>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/utils.hpp>


namespace fs = std::filesystem;
using namespace rs_toolset;

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "color balancing", "A color balancing tool using the global-to-local "
      "algorithm to create the color balanced rasters or the arguments json");
  options.add_options()
      ("h,help", "Print usage")
      ("input", "The input rasters' path or directories, only supports rasters"
          "with 3-band, 8-bit and GTiff format",
          cxxopts::value<std::vector<std::string>>())
      ("input-argus-json", "The input arguments JSON path",
          cxxopts::value<std::string>())
      ("output", "The output color balanced rasters' directory",
          cxxopts::value<std::string>())
      ("output-argus-json", "The output arguments JSON path",
          cxxopts::value<std::string>())
      ("block-size", "The block size per operation",
          cxxopts::value<int>()->default_value("16384"))
      ("log-level", "Available levels are trace(t), debug(d), info(i), "
          "warn(w), err(e), critical(c)", cxxopts::value<std::string>());
  auto result(options.parse(argc, argv));
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  if (!result.count("input") && !result.count("input-argus-json")) {
    std::cout << "Either the \"input\" argument or the \"input-argus-json\" "
    "should be specified" << std::endl;
    return -1;
  }
  if (result.count("input") && result.count("input-argus-json")) {
    std::cout << "The \"input\" and the \"input-argus-json\" argument can "
    "not be specified at the same time" << std::endl;
    return -1;
  }
  if (!result.count("output") && !result.count("output-argus-json")) {
    std::cout << "Either the \"onput\" argument or the \"onput-argus-json\" "
    "should be specified" << std::endl;
    return -1;
  }
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
    utils::InitSpdlog("Color Balancing");
  }

  utils::InitGdal(argv[0]);
  auto global_to_local(color_balancing::GlobalToLocal::Create(
      result["block-size"].as<int>()));
  if (result.count("input")) {
    std::regex gtiff_reg(".*\\.(tif|tiff|TIF|TIFF)");
    std::vector<std::string> rasters_path;
    for (const auto& path : result["input"].as<std::vector<std::string>>()) {
      if (!fs::exists(path)) continue;
      if (fs::is_directory(path)) {
        for (const auto& entry : fs::directory_iterator(path)) {
          auto subpath(entry.path());
          if (std::regex_match(subpath.filename().string(), gtiff_reg))
            rasters_path.push_back(subpath.string());
        }
      } else if (std::regex_match(path, gtiff_reg)) {
        rasters_path.push_back(path);
      }
    }
    global_to_local->CalcArgus(rasters_path);
  } else {
    global_to_local->ImportArgusJson(
        result["input-argus-json"].as<std::string>());
  }
  if (result.count("output"))
    global_to_local->CreateRasters({}, result["output"].as<std::string>());
  if (result.count("output-argus-json"))
    global_to_local->ExportArgusJson(
        result["output-argus-json"].as<std::string>());
  return 0;
}