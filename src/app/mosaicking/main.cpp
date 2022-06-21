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


int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "Mosaicking", 
      "A mosaicking tool using the graph cut algorithm to create the seamline "
      "vector(polygon elements) for the given rasters.");
  options.add_options()
      ("h,help", "Print usage")
      ("input", "The input rasters' path or directories", 
          cxxopts::value<std::vector<std::string>>())
      ("output", "The output composite table vector path",
          cxxopts::value<std::string>())
      ("epsg", "The spatial reference epsg code of the output composite table "
          "vector", cxxopts::value<int>()->default_value("4326"))
      ("grad-exp", "The gradient exponential argument in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("0.8"))
      ("min-diff", "The minimun difference argument in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("0.4"))
      ("max-diff", "The maximum difference argument in the "
          "graph cut algorithm", cxxopts::value<double>()->default_value("8"))
      ("unit", "The exported composite table unit value scaling buffer and "
          "tolerance", cxxopts::value<double>()->default_value("1e-05"))
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
  std::vector<std::string> rasters_path;
  std::string composite_table_path;
  if (result.count("input") && result.count("output")) {
    std::regex raster_reg(".*\\.(tif|tiff|TIF|TIFF)"), vector_reg(".*\\.(shp)");
    for (const auto& path : result["input"].as<std::vector<std::string>>()) {
      if (!std::filesystem::exists(path)) continue;
      if (std::filesystem::is_directory(path)) {
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
          auto subpath(entry.path());
          if (std::regex_match(subpath.filename().string(), raster_reg))
            rasters_path.push_back(subpath.string());
        }
      } else if (std::regex_match(path, raster_reg)) {
        rasters_path.push_back(path);
      }
    }
    composite_table_path = result["output"].as<std::string>();
    if (!std::regex_match(composite_table_path, vector_reg)) {
      std::cout << options.help() << std::endl;
      return 0;
    }
  } else {
    std::cout << options.help() << std::endl;
    return 0;
  }
  int epsg(result["epsg"].as<int>());
  double grad_exp(result["grad-exp"].as<double>()),
      min_diff(result["min-diff"].as<double>()),
      max_diff(result["max-diff"].as<double>()),
      unit(result["unit"].as<double>()),
      buffer(result["buffer"].as<double>()),
      tol(result["tol"].as<double>());
  if (result.count("log-level")) {
    std::string log_level(result["log-level"].as<std::string>());
    if (log_level == "t" || log_level == "trace") {
      rs_toolset::utils::InitSpdlog(
          "Mosaicking", spdlog::level::level_enum::trace);
    } else if (log_level == "d" || log_level == "debug") {
      rs_toolset::utils::InitSpdlog(
          "Mosaicking", spdlog::level::level_enum::debug);
    } else if (log_level == "i" || log_level == "info") {
      rs_toolset::utils::InitSpdlog(
          "Mosaicking", spdlog::level::level_enum::info);
    } else if (log_level == "w" || log_level == "warn") {
      rs_toolset::utils::InitSpdlog(
          "Mosaicking", spdlog::level::level_enum::warn);
    } else if (log_level == "e" || log_level == "err") {
      rs_toolset::utils::InitSpdlog(
          "Mosaicking", spdlog::level::level_enum::err);
    } else if (log_level == "c" || log_level == "critical") {
      rs_toolset::utils::InitSpdlog(
          "Mosaicking", spdlog::level::level_enum::critical);
    }
  } else {
    rs_toolset::utils::InitSpdlog("Mosaicking");
  }

  rs_toolset::utils::InitGdal(argv[0]);
  std::unique_ptr<OGRSpatialReference> spatial_ref(new OGRSpatialReference);
  spatial_ref->importFromEPSG(epsg);
  spatial_ref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  auto mosaicking(rs_toolset::mosaicking::MosaickingContainer::Create(
      rs_toolset::mosaicking::GraphCut::Create(grad_exp, min_diff, max_diff), 
      spatial_ref.get()));
  mosaicking->SortRasters(rasters_path);
  for (int i = 0; i < rasters_path.size(); i++) {
    mosaicking->AddTask(rasters_path[i], true);
    spdlog::info(
        "----------- {}/{} - done ----------", i + 1, rasters_path.size());
  }
  mosaicking->ExportCompositeTableVector(
      composite_table_path, unit, buffer, tol);
}