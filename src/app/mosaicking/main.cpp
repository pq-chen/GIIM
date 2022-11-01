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
      "mosaicking", "A mosaicking tool using the available methods to create a "
      "mosaicking vector(polygon elements) or raster from the given input "
      "DOM rasters");
  options.add_options("basic")
      ("h,help", "Print usage")
      ("method", "Available methods include \"serial\" and \"parallel\"", cxxopts::value<std::string>()->default_value("serial"))
      ("input", "The input rasters' path or directories, only supporting DOM rasters with 8-bit", cxxopts::value<std::vector<std::string>>())
      ("output-epsg", "The output spatial reference EPSG code", cxxopts::value<int>()->default_value("4326"))
      ("output-composite-table", "The output composite table path", cxxopts::value<std::string>())
      ("log-level", "Available levels are trace(t), debug(d), info(i), warn(w), err(e), critical(c), off(o)", cxxopts::value<std::string>()->default_value("i"));
  options.add_options("mosaicking-raster-feature")
      ("output-mosaicking-raster", "The output mosaicking raster path", cxxopts::value<std::string>())
      ("output-reso", "The output mosaicking raster resolution", cxxopts::value<double>())
      ("output-blend-dist", "The output mosaicking raster blend distance in pixels", cxxopts::value<double>()->default_value("0.0"));
  options.add_options("common-feature")
      ("color-balancing-json", "The color balancing json file path", cxxopts::value<std::string>())
      ("low-overview-trunc", "The low overview trunction", cxxopts::value<int>()->default_value("3"))
      ("high-overview-trunc", "The high overview trunction", cxxopts::value<int>()->default_value("1"))
      ("rgb-bands-map", "The RGB bands' map", cxxopts::value<std::vector<int>>());
  options.add_options("serial-feature")
      ("input-composite-table", "The input composite table path initializing the mosaicking container", cxxopts::value<std::string>())
      ("input-rasters-dir", "The input rasters' directory corresponding to the input composite table", cxxopts::value<std::string>())
      ("rejection-ratio", "The regularization term to reject a new item in the output composite table", cxxopts::value<double>()->default_value("0.01"))
      ("check-interval", "The interval of creating a check composite table", cxxopts::value<int>())
      ("disable-extension", "Disable using the raster name extension in the output composite table", cxxopts::value<bool>()->default_value("false"))
      ("disable-sort", "Disable sorting the input rasters by resolution", cxxopts::value<bool>()->default_value("false"))
      ("query-composite-table", "The query composite table path", cxxopts::value<std::string>())
      ("query-rasters-name-field-name", "The query composite table field name representing the raster name", cxxopts::value<std::string>());
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
  std::string output_composite_table_path;
  if (result.count("output-composite-table")) {
    if (output_composite_table_path =
            result["output-composite-table"].as<std::string>();
        output_composite_table_path.empty()) {
      std::cout << "The \"output-composite-table\" argument is empty"
          << std::endl;
      return -1;
    }
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
  std::string output_mosaicking_raster_path;
  double output_reso;
  if (result.count("output-mosaicking-raster")) {
    if (output_mosaicking_raster_path =
            result["output-mosaicking-raster"].as<std::string>();
        output_mosaicking_raster_path.empty()) {
      std::cout << "The \"output-mosaicking-raster\" argument is empty"
          << std::endl;
      return -1;
    }
    if (result.count("output-reso")) {
      output_reso = result["output-reso"].as<double>();
    } else {
      std::cout << "Lack of the \"output-reso\" argument "
          "if uses the \"output-mosaicking-raster\" argument" << std::endl;
      return -1;
    }
  } else if (result.count("output-reso")) {
    std::cout << "Lack of the \"output-mosaicking-raster\" argument "
        "if uses the \"output-reso\" argument" << std::endl;
    return -1;
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
  std::string input_composite_table_path, input_rasters_dir;
  if (result.count("input-composite-table")) {
    if (input_composite_table_path =
            result["input-composite-table"].as<std::string>();
        input_composite_table_path.empty()) {
      std::cout << "The \"input-composite-table\" argument is empty"
          << std::endl;
      return -1;
    }
    if (result.count("input-rasters-dir")) {
      input_rasters_dir = result["input-rasters-dir"].as<std::string>();
    } else {
      std::cout << "Lack of the \"input-rasters-dir\" argument "
          "if uses the \"input-composite-table\" argument" << std::endl;
      return -1;
    }
  } else if (result.count("input-rasters-dir")) {
    std::cout << "Lack of the \"input-composite-table\" argument "
        "if uses the \"input-rasters-dir\" argument" << std::endl;
    return -1;
  }
  std::string query_composite_table_path, query_rasters_name_field_name;
  if (result.count("query-composite-table")) {
    if (query_composite_table_path =
            result["query-composite-table"].as<std::string>();
        query_composite_table_path.empty()) {
      std::cout << "The \"query-composite-table\" argument is empty"
          << std::endl;
      return -1;
    }
    if (result.count("query-rasters-name-field-name")) {
      query_rasters_name_field_name =
          result["query-rasters-name-field-name"].as<std::string>();
    } else {
      std::cout << "Lack of the \"query-rasters-name-field-name\" argument "
          "if uses the \"query-composite-table\" argument" << std::endl;
      return -1;
    }
  } else if (result.count("query-rasters-name-field-name")) {
    std::cout << "Lack of the \"query-composite-table\" argument "
        "if uses the \"query-rasters-name-field-name\" argument" << std::endl;
    return -1;
  }
  auto mosaicking(mosaicking::GraphCut::Create(
      result["grad-self-low"].as<float>(), result["grad-self-high"].as<float>(),
      result["grad-self-exp"].as<float>(), result["diff-low"].as<float>(),
      result["diff-exp"].as<float>(), result["tol"].as<double>()));
  
  GDALDatasetUniquePtr composite_table_dataset(nullptr);
  if (std::string method(result["method"].as<std::string>());
      method == "serial") {
    std::shared_ptr<mosaicking::MosaickingContainer> mosaicking_container(
        nullptr);
    if (input_composite_table_path.empty()) {
      auto spatial_ref(std::make_unique<OGRSpatialReference>());
      spatial_ref->importFromEPSG(result["output-epsg"].as<int>());
      spatial_ref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
      mosaicking_container = mosaicking::MosaickingContainer::Create(
          mosaicking, spatial_ref.get(), result["rejection-ratio"].as<double>(),
          color_balancing);
    } else {
      mosaicking_container = mosaicking::MosaickingContainer::Create(
          mosaicking, input_composite_table_path, input_rasters_dir,
          result["rejection-ratio"].as<double>(), color_balancing);
    }
    auto existing_rasters_name(mosaicking_container->ExportAllRastersName());
    auto begin(existing_rasters_name.begin()), end(existing_rasters_name.end());
    std::vector<std::string> paths;
    if (result.count("input")) {
      for (const auto& path : result["input"].as<std::vector<std::string>>()) {
        if (!fs::exists(path)) {
          std::cout << "Input path " << path << " does not exist" << std::endl;
          return -1;
        }
        if (fs::is_directory(path)) {
          for (const auto& entry : fs::directory_iterator(path)) {
            if (utils::GetRasterDriverByPath(entry.path().string()) &&
                std::find(
                    begin, end, entry.path().filename().string()) == end) {
              paths.push_back(entry.path().string());
            }
          }
        } else if (utils::GetRasterDriverByPath(path) &&
            std::find(begin, end, fs::path(path).filename().string()) == end) {
          paths.push_back(path);
        }
      }
      std::cout << "The following raster(s) will be operated:" << std::endl;
      for (const auto& path : paths)
        std::cout << path << std::endl;
      std::cout << paths.size() << " task(s) in total" << std::endl;
    }
    if (!result["disable-sort"].as<bool>())
      mosaicking_container->SortRasters(paths);
    for (int i(0); i < paths.size(); ++i) {
      if (!mosaicking_container->AddTask(
              paths[i], result["low-overview-trunc"].as<int>(),
              result["high-overview-trunc"].as<int>(), rgb_bands_map)) {
        return -1;
      }
      spdlog::info("---------- {}/{} - done ----------", i + 1, paths.size());
      if (result.count("check-interval") &&
          (i + 1) % result["check-interval"].as<int>() == 0 &&
          (i + 1) != paths.size()) {
        if (auto pos(output_composite_table_path.rfind('.'));
            !mosaicking_container->ExportCompositeTable(
                output_composite_table_path.substr(0, pos) + '_' +
                    std::to_string(i + 1) +
                    output_composite_table_path.substr(pos),
                query_composite_table_path, query_rasters_name_field_name)) {
          return -1;
        }
      }
    }
    if (!output_composite_table_path.empty()) {
      if (composite_table_dataset = mosaicking_container->ExportCompositeTable(
              output_composite_table_path, query_composite_table_path,
              query_rasters_name_field_name,
              !result["disable-extension"].as<bool>());
          !composite_table_dataset) {
        return -1;
      }
    }
  } else if (method == "parallel") {
    auto voronoi_diagrams(mosaicking::VoronoiDiagrams::Create(mosaicking));
    auto spatial_ref(std::make_unique<OGRSpatialReference>());
    spatial_ref->importFromEPSG(result["output-epsg"].as<int>());
    spatial_ref->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::vector<std::string> paths;
    if (result.count("input")) {
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
    }
    if (!output_composite_table_path.empty()) {
      if (composite_table_dataset = voronoi_diagrams->Run(
              paths, output_composite_table_path, spatial_ref.get(), true,
              result["low-overview-trunc"].as<int>(),
              result["high-overview-trunc"].as<int>(), rgb_bands_map,
              color_balancing);
          !composite_table_dataset)  {
        return -1;
      }
    }
  } else {
    std::cout << "Available methods include \"serial\" and \"parallel\""
        << std::endl;
    return -1;
  }
  if (!output_mosaicking_raster_path.empty()) {
    if (auto composite_table_layer(composite_table_dataset->GetLayer(0));
        !mosaicking::CreateMosaickingRaster(
            output_mosaicking_raster_path, composite_table_layer, output_reso,
            color_balancing, result["output-blend-dist"].as<double>())) {
      return -1;
    }
    GDALDatasetUniquePtr dataset(GDALDataset::Open(
        output_mosaicking_raster_path.c_str(),
        GDAL_OF_RASTER | GDAL_OF_READONLY));
    utils::CreateRasterPyra(dataset.get());
  }
  return 0;
}