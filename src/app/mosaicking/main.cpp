#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>
#include <ogrsf_frmts.h>

#include <rs-toolset/mosaicking.h>
#include <rs-toolset/utils.hpp>


using namespace rs_toolset;

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "mosaicking", "A mosaicking tool using the graph cut mosaicking "
      "algorithm to create the seamline vector(polygon elements) from the "
      "given rasters");
  options.add_options()
      ("h,help", "Print usage")
      ("input", "The input rasters' path or directories, only supports DOM "
          "and tif files", cxxopts::value<std::vector<std::string>>())
      ("output-composite-table", "The output composite table path",
          cxxopts::value<std::string>())
      ("output-mosaicking-raster", "The output mosaicking raster path",
          cxxopts::value<std::string>())
      ("output-reso", "The resolution for the output mosaicking raster",
          cxxopts::value<double>())
      ("output-blend-dist", "The output blend distance for the mosaicking "
          "raster", cxxopts::value<double>()->default_value("0.0"))
      ("check-freq", "The frequency of creating the checking composite "
          "table", cxxopts::value<int>())
      ("input-composite-table", "Initializing the mosaicking container by the "
          "input composite table", cxxopts::value<std::string>())
      ("input-rasters-dir",  "The input rasters' directory corresponding to "
          "the input composit table", cxxopts::value<std::string>())
      ("last-overview-idx", "The last overview index operated on", 
          cxxopts::value<int>()->default_value("3"))
      ("epsg", "The spatial reference epsg code of the output composite table",
          cxxopts::value<int>()->default_value("4326"))
      ("grad-exp", "The gradient term exponential argument in the graph cut "
          "algorithm", cxxopts::value<double>()->default_value("0.8"))
      ("diff-low", "The difference term low trunction argument in the graph "
          "cut algorithm", cxxopts::value<double>()->default_value("0.4"))
      ("diff-high", "The difference term high trunction in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("8"))
      ("tol", "The tolerance for simplifyling the seamline",
          cxxopts::value<double>()->default_value("2.0"))
      ("log-level", "Available levels are trace(t), debug(d), info(i), "
          "warn(w), err(e), critical(c)", cxxopts::value<std::string>());
  auto result = options.parse(argc, argv);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  std::regex gtiff_reg(".*\\.(tif|tiff|TIF|TIFF)");
  std::string output_mosaicking_raster_path;
  double output_reso;
  if (result.count("output-mosaicking-raster")) {
    if (!std::regex_match(
          result["output-mosaicking-raster"].as<std::string>(), gtiff_reg)) {
      std::cout << "The output mosaicking raster path does not end with "
          "\".tif/.tiff/.TIF/.TIFF\"" << std::endl;
      return 0;
    }
    if (result.count("output-reso")) {
      output_mosaicking_raster_path =
          result["output-mosaicking-raster"].as<std::string>();
      output_reso = result["output-reso"].as<double>();
    } else {
      std::cout << "Lack of the \"output-reso\" argument "
          "if uses the \"output-mosaicking-raster\" argument" << std::endl;
      return 0;
    }
  } else if (result.count("output-reso")) {
    std::cout << "Lack of the \"output-mosaicking-raster\" argument "
        "if uses the \"output-reso\" argument" << std::endl;
    return 0;
  }
  std::regex esri_shapefile_reg(".*\\.(shp)");
  std::string input_composite_table_path, input_rasters_dir;
  if (result.count("input-composite-table")) {
    if (!std::regex_match(
          result["input-composite-table"].as<std::string>(),
          esri_shapefile_reg)) {
      std::cout << "The input composite table path does not end with \".shp\""
          << std::endl;
      return 0;
    }
    if (result.count("input-rasters-dir")) {
      input_composite_table_path =
          result["input-composite-table"].as<std::string>();
      input_rasters_dir = result["input-rasters-dir"].as<std::string>();
    } else {
      std::cout << "Lack of the \"input-rasters-dir\" argument "
          "if uses the \"input-composite-table\" argument" << std::endl;
    }
  } else if (result.count("rasters-dir")) {
    std::cout << "Lack of the \"input-composite-table\" argument "
        "if uses the \"input-rasters-dir\" argument" << std::endl;
    return 0;
  }
  auto epsg(result["epsg"].as<int>());
  auto grad_term_exp(result["grad-exp"].as<double>()),
      diff_term_low_trunc(result["diff-low"].as<double>()),
      diff_term_high_trunc(result["diff-high"].as<double>()),
      tol(result["tol"].as<double>());
  if (result.count("log-level")) {
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
    }
  } else {
    utils::InitSpdlog("Mosaicking");
  }

  utils::InitGdal(argv[0]);
  auto spatial_ref(std::make_unique<OGRSpatialReference>());
  spatial_ref->importFromEPSG(epsg);
  spatial_ref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  std::shared_ptr<mosaicking::MosaickingContainer> mosaicking(nullptr);
  if (input_composite_table_path.empty()) {
    mosaicking = mosaicking::MosaickingContainer::Create(
        mosaicking::GraphCut::Create(
            grad_term_exp, diff_term_low_trunc, diff_term_high_trunc, tol),
        spatial_ref.get());
  } else {
    mosaicking = mosaicking::MosaickingContainer::Create(
        mosaicking::GraphCut::Create(
            grad_term_exp, diff_term_low_trunc, diff_term_high_trunc, tol),
        result["input-composite-table"].as<std::string>(),
        result["input-rasters-dir"].as<std::string>());
  }
  std::vector<std::string> existing_rasters_name(
      mosaicking->ExportAllRastersName());
  std::vector<std::string> rasters_path;
  if (result.count("input")) {
    for (const auto& path : result["input"].as<std::vector<std::string>>()) {
      if (!std::filesystem::exists(path)) continue;
      if (std::filesystem::is_directory(path)) {
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
          auto subpath(entry.path());
          if (std::regex_match(subpath.filename().string(), gtiff_reg) &&
              std::find(
                existing_rasters_name.begin(), existing_rasters_name.end(),
                  subpath.filename().string()) == existing_rasters_name.end())
          rasters_path.push_back(subpath.string());
        }
      } else if (std::regex_match(path, gtiff_reg) &&
            std::find(
                existing_rasters_name.begin(), existing_rasters_name.end(),
                std::filesystem::path(path).filename().string()) ==
                existing_rasters_name.end()) {
        rasters_path.push_back(path);
      }
    }
    std::cout << "The following raster(s) will be operated:" << std::endl;
    for (const auto& raster_path : rasters_path)
      std::cout << raster_path << std::endl;
    std::cout << rasters_path.size() << " task(s) in total" << std::endl;
  }
  std::string output_composite_table_path;
  if (result.count("output-composite-table")) {
    if (!std::regex_match(
          result["output-composite-table"].as<std::string>(),
          esri_shapefile_reg)) {
      std::cout << "The output composite table path does not end with \".shp\""
          << std::endl;
      return 0;
    }
    output_composite_table_path =
        result["output-composite-table"].as<std::string>();
  } else if (!rasters_path.empty()) {
    std::cout << 
        "Lack of the \"output-composite-table\" argument" << std::endl;
    return 0;
  }
  int check_freq;
  if (result.count("check-freq")) {
    check_freq = result["check-freq"].as<int>();
  } else {
    check_freq = static_cast<int>(rasters_path.size());
  }
  mosaicking->SortRasters(rasters_path);
  for (int i = 0; i < rasters_path.size(); i++) {
    mosaicking->AddTask(
        rasters_path[i], result["last-overview-idx"].as<int>(), true);
    spdlog::info(
        "----------- {}/{} - done ----------", i + 1, rasters_path.size());
    if ((i + 1) % check_freq == 0 && (i + 1) != rasters_path.size()) {
      size_t pos(output_composite_table_path.size() - 4);
      mosaicking->ExportCompositeTableVector(
          output_composite_table_path.substr(0, pos) + '_' +
          std::to_string(i + 1) + ".shp");
    }
  }
  if (!output_composite_table_path.empty()) {
    mosaicking->ExportCompositeTableVector(output_composite_table_path);
  }
  if (!output_mosaicking_raster_path.empty()) {
    mosaicking->CreateMosaickingRaster(
        output_mosaicking_raster_path, "", "", output_reso,
        result["blend-dist"].as<double>());
  }
}