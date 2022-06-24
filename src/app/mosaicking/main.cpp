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
      "Mosaicking", 
      "A mosaicking tool using the graph cut algorithm to create the seamline "
      "vector(polygon elements) for the given rasters.");
  options.add_options()
      ("h,help", "Print usage")
      ("input", "The input rasters' path or directories", 
          cxxopts::value<std::vector<std::string>>())
      ("output-composite-table", "The output composite table vector path",
          cxxopts::value<std::string>())
      ("output-mosaicking-raster", "The output mosaicking raster path",
          cxxopts::value<std::string>())
      ("output-reso", "The output mosaicking raster's resolution",
          cxxopts::value<double>())
      ("check-freq", "The frequency of creating the checking composite "
          "table vector", cxxopts::value<int>())
      ("input-composite-table", "The input composite table path for creating "
          "the mosaicking container", cxxopts::value<std::string>())
      ("input-rasters-dir",  "The input rasters' directory corresponding to "
          "the input composit table", cxxopts::value<std::string>())
      ("epsg", "The spatial reference epsg code of the output composite table "
          "vector", cxxopts::value<int>()->default_value("4326"))
      //("c,conservative",
      //    "Use conservative strategy for the graph cut algorithm")
      ("grad-exp", "The gradient term exponential argument in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("0.8"))
      ("diff-low", "The difference term low trunction argument in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("0.4"))
      ("diff-high", "The difference term high trunction in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("8"))
      ("unit", "The exported composite table unit value for scaling buffer and "
          "tolerance", cxxopts::value<double>()->default_value("0"))
      ("buffer", "The exported composite table buffer distance", 
          cxxopts::value<double>()->default_value("-1"))
      ("tol", "The exported composite table tolerance for simplification", 
          cxxopts::value<double>()->default_value("1.5"))
      ("log-level", "Available levels are trace(t), debug(d), info(i), warn(w),"
          " err(e), critical(c)", cxxopts::value<std::string>());
  auto result = options.parse(argc, argv);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  std::string output_composite_table_path;
  std::regex esri_shapefile_reg(".*\\.(shp)");
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
  } else {
    std::cout << "Lack of the \"output-composite-table\" argument" << std::endl;
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
      unit(result["unit"].as<double>()),
      buffer(result["buffer"].as<double>()),
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
            grad_term_exp, diff_term_low_trunc, diff_term_high_trunc),
        spatial_ref.get());
  } else {
    mosaicking = mosaicking::MosaickingContainer::Create(
        mosaicking::GraphCut::Create(
              grad_term_exp, diff_term_low_trunc, diff_term_high_trunc),
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
  int check_freq;
  if (result.count("check-freq")) {
    check_freq = result["check-freq"].as<int>();
  } else {
    check_freq = static_cast<int>(rasters_path.size());
  }
  mosaicking->SortRasters(rasters_path);
  for (int i = 0; i < rasters_path.size(); i++) {
    mosaicking->AddTask(rasters_path[i], true);
    spdlog::info(
        "----------- {}/{} - done ----------", i + 1, rasters_path.size());
    if ((i + 1) == rasters_path.size()) {
      mosaicking->ExportCompositeTableVector(
          output_composite_table_path, unit, buffer, tol);
    } else if ((i + 1) % check_freq == 0) {
      size_t pos(output_composite_table_path.size() - 4);
      mosaicking->ExportCompositeTableVector(
          output_composite_table_path.substr(0, pos) + '_' + 
          std::to_string(i+ 1) + ".shp", unit, buffer, tol);
    }
  }
  if (!output_mosaicking_raster_path.empty()) {
    mosaicking->CreateMosaickingRaster(
        output_mosaicking_raster_path, "", "", output_reso);
  }
}