#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/color_balancing.h>
#include <rs-toolset/mosaicking.h>
#include <rs-toolset/utils.hpp>

namespace fs = std::filesystem;
using namespace rs_toolset;

int main(int argc, char* argv[]) {
  cxxopts::Options options(
      "RS-Toolset_app_mosaicking", "A mosaicking tool using the available "
      "methods to create a mosaicking vector(polygon) or raster from the given "
      "input DOMs");
  options.add_options("basic")
      ("h,help", "Print usage")
      ("method", "Available mosaicking methods are \"serial\" and \"voronoi\"", cxxopts::value<std::string>()->default_value("serial"))
      ("input", "The input rasters' path and directories, only supporting 8-bit DOMs", cxxopts::value<std::vector<std::string>>())
      ("epsg", "The output spatial reference EPSG code", cxxopts::value<int>()->default_value("4326"))
      ("log-level", "Available log levels are \"trace(t)\", \"debug(d)\", \"info(i)\", \"warn(w)\", \"err(e)\", \"critical(c)\" and \"off(o)\"", cxxopts::value<std::string>()->default_value("i"));
  options.add_options("output vector")
      ("output-mosaicking-vector", "The output mosaicking vector path", cxxopts::value<std::string>())
      ("output-border-vector", "The output border vector path", cxxopts::value<std::string>());
  options.add_options("output raster")
      ("output-mosaicking-raster", "The output mosaicking raster path", cxxopts::value<std::string>())
      ("reso", "The output mosaicking raster resolution", cxxopts::value<double>())
      ("blend-dist", "The output mosaicking raster blend distance in pixels", cxxopts::value<double>()->default_value("0.0"))
      ("output-rasters-cut-dir", "The output rasters' cut directory", cxxopts::value<std::string>())
      ("compression", "The raster compression method", cxxopts::value<std::string>()->default_value(""));
  options.add_options("common")
      ("color-balancing-json", "The color balancing json file path", cxxopts::value<std::string>())
      ("low-overviews-trunc", "The low overviews trunction", cxxopts::value<int>()->default_value("3"))
      ("high-overviews-trunc", "The high overviews trunction", cxxopts::value<int>()->default_value("1"))
      ("rgb-bands-map", "The RGB bands' map", cxxopts::value<std::vector<int>>());
  options.add_options("serial")
      ("former-mosaicking-vector", "The former mosaicking vector path for initializition", cxxopts::value<std::string>())
      ("former-rasters-dir", "The former rasters' directory corresponding to the former mosaicking vector", cxxopts::value<std::string>())
      ("surrounded-buffer", "The buffer distance in pixels in the (anti-)surrounded case", cxxopts::value<double>()->default_value("100.0"))
      ("rejection-ratio", "The regularization term to reject creating a new item", cxxopts::value<double>()->default_value("0.0005"))
      ("anti-surrounded", "Whether deals with the anti-surrounded case or not", cxxopts::value<bool>()->default_value("false"))
      ("check-interval", "The interval of creating a check mosaicking vector", cxxopts::value<int>())
      ("disable-extension", "Disable using the raster name extension", cxxopts::value<bool>()->default_value("false"))
      ("disable-sort", "Disable sorting the input rasters", cxxopts::value<bool>()->default_value("false"))
      ("query-mosaicking-vector", "The query mosaicking vector path", cxxopts::value<std::string>())
      ("query-raster-name-field-name", "The query mosaicking vector field name representing the raster name", cxxopts::value<std::string>());
  options.add_options("graph cut")
      ("grad-self-low", "The low trunction of the gradient-self term", cxxopts::value<float>()->default_value("1.4"))
      ("grad-self-high", "The high trunction of the gradient-self term", cxxopts::value<float>()->default_value("40.0"))
      ("grad-self-exp", "The exponential of the gradient-self term", cxxopts::value<float>()->default_value("1.5"))
      ("diff-low", "The low trunction of the difference term", cxxopts::value<float>()->default_value("20.0"))
      ("diff-exp", "The exponential of the difference term", cxxopts::value<float>()->default_value("0.7"))
      ("tol", "The tolerance in pixels for simplifyling the seamline", cxxopts::value<double>()->default_value("2.0"));
  auto result(options.parse(argc, argv));

  utils::InitGdal(argv[0]);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  auto method(result["method"].as<std::string>());
  if (method != "serial" && method != "voronoi") {
    std::cout << "Available mosaicking methods are \"serial\" and \"voronoi\""
        << std::endl;
    return -1;
  }
  auto spatial_ref(std::make_unique<OGRSpatialReference>());
  if (auto epsg(result["epsg"].as<int>());
      spatial_ref->importFromEPSG(epsg) != OGRERR_NONE) {
    std::cout << "Importing spatial reference from EPSG code " << epsg
        << " failed" << std::endl;
    return -1;
  }
  spatial_ref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
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
    std::cout << "Available log levels are \"trace(t)\", \"debug(d)\", "
        "\"info(i)\", \"warn(w)\", \"err(e)\", \"critical(c)\" and \"off(o)\""
        << std::endl;
    return -1;
  }
  std::string output_mosaicking_vector_path;
  if (result.count("output-mosaicking-vector")) {
    if (output_mosaicking_vector_path =
            result["output-mosaicking-vector"].as<std::string>();
        output_mosaicking_vector_path.empty()) {
      std::cout << "The \"output-mosaicking-vector\" is empty" << std::endl;
      return -1;
    }
  }
  std::string output_border_vector_path;
  if (result.count("output-border-vector")) {
    if (output_border_vector_path =
            result["output-border-vector"].as<std::string>();
        output_border_vector_path.empty()) {
      std::cout << "The \"output-border-vector\" is empty" << std::endl;
      return -1;
    }
  }
  std::string output_mosaicking_raster_path;
  double reso;
  if (result.count("output-mosaicking-raster") ^ result.count("reso")) {
    std::cout << "The \"output-mosaicking-raster\" and the \"reso\" must be "
        "used together" << std::endl;
    return -1;
  } else if (result.count("output-mosaicking-raster")) {
    output_mosaicking_raster_path =
        result["output-mosaicking-raster"].as<std::string>();
    reso = result["reso"].as<double>();
  }
  std::string output_rasters_cut_dir;
  if (result.count("output-rasters-cut-dir")) {
    output_rasters_cut_dir = result["output-rasters-cut-dir"].as<std::string>();
    if (!fs::is_directory(output_rasters_cut_dir) ||
        !fs::exists(output_rasters_cut_dir)) {
      std::cout << "The \"output-rasters-cut-dir\" " << output_rasters_cut_dir
          << " is not a directory or does not exist" << std::endl;
      return -1;
    }
  }
  std::shared_ptr<color_balancing::ColorBalancingInterface> color_balancing(
      nullptr);
  if (result.count("color-balancing-json")) {
    if (color_balancing =
            color_balancing::ColorBalancingInterface::CreateFromArgusJson(
                result["color-balancing-json"].as<std::string>());
        !color_balancing) {
      return -1;
    }
  }
  std::vector<int> rgb_bands_map;
  if (result.count("rgb-bands-map"))
    rgb_bands_map = result["rgb-bands-map"].as<std::vector<int>>();
  std::string former_mosaicking_vector_path, former_rasters_dir;
  if (result.count("former-mosaicking-vector") ^
      result.count("former-rasters-dir")) {
    std::cout << "The \"former-mosaicking-vector\" and the "
        "\"former-rasters-dir\" must be used together" << std::endl;
    return -1;
  } else if (result.count("former-mosaicking-vector")) {
    former_mosaicking_vector_path =
        result["former-mosaicking-vector"].as<std::string>();
    former_rasters_dir = result["former-rasters-dir"].as<std::string>();
  }
  auto check_interval(-1);
  if (result.count("check-interval")) {
    if (check_interval = result["check-interval"].as<int>();
        check_interval <= 0) {
      std::cout << "The \"check-interval\" " << check_interval
          << " must be positive" << std::endl;
      return -1;
    }
  }
  std::string query_mosaicking_vector_path, query_raster_name_field_name;
  if (result.count("query-mosaicking-vector") ^
      result.count("query-raster-name-field-name")) {
    std::cout << "The \"query-mosaicking-vector\" and the "
        "\"query-raster-name-field-name\" must be used together" << std::endl;
    return -1;
  } else if (result.count("query-mosaicking-vector")) {
    query_mosaicking_vector_path =
        result["query-mosaicking-vector"].as<std::string>();
    query_raster_name_field_name =
        result["query-raster-name-field-name"].as<std::string>();
  }
  auto mosaicking(mosaicking::GraphCut::Create(
      result["grad-self-low"].as<float>(), result["grad-self-high"].as<float>(),
      result["grad-self-exp"].as<float>(), result["diff-low"].as<float>(),
      result["diff-exp"].as<float>(), result["tol"].as<double>()));
  if (!mosaicking) return -1;
  
  GDALDatasetUniquePtr mosaicking_vector_dataset(nullptr);
  std::vector<std::string> all_paths;
  if (method == "serial") {
    std::shared_ptr<mosaicking::SerialContainer> serial_container(nullptr);
    if (former_mosaicking_vector_path.empty()) {
      serial_container = mosaicking::SerialContainer::Create(
          mosaicking, spatial_ref.get(), color_balancing);
    } else {
      serial_container = mosaicking::SerialContainer::Create(
          mosaicking, former_mosaicking_vector_path, former_rasters_dir,
          color_balancing);
    }
    if (!serial_container) return -1;
    auto former_rasters_name(serial_container->ExportAllRastersPath());
    for (auto& name : former_rasters_name)
      name = fs::path(name).filename().string();
    auto begin(former_rasters_name.begin()), end(former_rasters_name.end());
    std::vector<std::string> input_paths;
    if (result.count("input")) {
      for (const auto& path : result["input"].as<std::vector<std::string>>()) {
        if (!fs::exists(path)) {
          std::cout << "The input path " << path << " does not exist"
              << std::endl;
          return -1;
        }
        if (fs::is_directory(path)) {
          for (const auto& entry : fs::directory_iterator(path)) {
            if (auto subpath(entry.path());
                utils::GetRasterDriverByPath(subpath.string()) &&
                std::find(begin, end, subpath.filename().string()) == end &&
                subpath.string().find("mask") == std::string::npos) {
              input_paths.push_back(subpath.string());
            }
          }
        } else if (utils::GetRasterDriverByPath(path) &&
            std::find(begin, end, fs::path(path).filename().string()) == end) {
          input_paths.push_back(path);
        }
      }
      std::cout << "The following raster(s) will be operated:" << std::endl;
      for (const auto& path : input_paths)
        std::cout << path << std::endl;
      std::cout << input_paths.size() << " task(s) in total" << std::endl;
    }
    if (!result["disable-sort"].as<bool>())
      serial_container->SortRasters(input_paths);
    for (auto i(0); i < input_paths.size(); ++i) {
      if (!serial_container->AddTask(
              input_paths[i], result["low-overviews-trunc"].as<int>(),
              result["high-overviews-trunc"].as<int>(),
              result["surrounded-buffer"].as<double>(),
              result["rejection-ratio"].as<double>(),
              result["anti-surrounded"].as<bool>(), rgb_bands_map)) {
        return -1;
      }
      spdlog::info(
          "---------- {}/{} - done ----------", i + 1, input_paths.size());
      if (check_interval != -1 && (i + 1) % check_interval == 0) {
        if (auto pos(output_mosaicking_vector_path.rfind('.'));
            !serial_container->ExportMosaickingVector(
                output_mosaicking_vector_path.substr(0, pos) + '_' +
                    std::to_string(i + 1) +
                    output_mosaicking_vector_path.substr(pos), "", "", true)) {
          return -1;
        }
      }
    }
    if (mosaicking_vector_dataset = serial_container->ExportMosaickingVector(
            output_mosaicking_vector_path, query_mosaicking_vector_path,
            query_raster_name_field_name,
            !result["disable-extension"].as<bool>());
        !mosaicking_vector_dataset) {
      return -1;
    } else if (!output_mosaicking_vector_path.empty()) {
      mosaicking_vector_dataset.reset(nullptr);
      mosaicking_vector_dataset.reset(GDALDataset::Open(
          output_mosaicking_vector_path.c_str(),
          GDAL_OF_VECTOR | GDAL_OF_READONLY));
    }
    if (!output_border_vector_path.empty() &&
        !serial_container->ExportBorderVector(output_border_vector_path)) {
      return -1;
    }
    all_paths = serial_container->ExportAllRastersPath();
  } else if (method == "voronoi") {
    auto voronoi_diagrams(mosaicking::VoronoiDiagrams::Create(mosaicking));
    if (result.count("input")) {
      for (const auto& path : result["input"].as<std::vector<std::string>>()) {
        if (!fs::exists(path)) {
          std::cout << "The input path " << path << " does not exist"
              << std::endl;
          return -1;
        }
        if (fs::is_directory(path)) {
          for (const auto& entry : fs::directory_iterator(path)) {
            if (auto subpath(entry.path().string());
                utils::GetRasterDriverByPath(subpath)) {
              all_paths.push_back(subpath);
            }
          }
        } else if (utils::GetRasterDriverByPath(path)) {
          all_paths.push_back(path);
        }
      }
    }
    if (mosaicking_vector_dataset = voronoi_diagrams->Run(
            all_paths, output_mosaicking_vector_path, spatial_ref.get(), true,
            result["low-overviews-trunc"].as<int>(),
            result["high-overviews-trunc"].as<int>(), rgb_bands_map,
            color_balancing); !mosaicking_vector_dataset)  {
      return -1;
    }
  }
  auto mosaicking_layer(mosaicking_vector_dataset->GetLayer(0));
  if (!output_mosaicking_raster_path.empty() &&
      !mosaicking::CreateMosaickingRaster(
            mosaicking_layer, all_paths, output_mosaicking_raster_path, reso,
            result["compression"].as<std::string>(),
            result["blend-dist"].as<double>(), rgb_bands_map,
            color_balancing)) {
    return -1;
  }
  if (!output_rasters_cut_dir.empty() &&
      !mosaicking::CreateRastersCut(
            mosaicking_layer, all_paths, output_rasters_cut_dir,
            query_raster_name_field_name,
            result["compression"].as<std::string>(), rgb_bands_map)) {
    return -1;
  }
  return 0;
}