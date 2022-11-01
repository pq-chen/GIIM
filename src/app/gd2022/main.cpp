#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <gdal_priv.h>
#include <nlohmann/json.hpp>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>


namespace fs = std::filesystem;
using namespace rs_toolset;

bool RunTask(
    const std::string& former_path,
    const std::string& later_path,
    const std::string& output_path,
    const std::map<std::string, std::string>& former_field_names,
    const std::map<std::string, std::string>& later_field_names,
    const std::string& later_raster_name_field_name,
    const std::string& masks_dir,
    const std::string& mask_reg_suffix) {
  std::string string(
      "Running a GuangDong2022 task with\n"
      " - Former path: {}\n"
      " - Later path: {}\n"
      " - Later raster name's field name: {}\n"
      " - Masks' directory: {}\n"
      " - Mask regular expression suffix: {}\n"
      " - Former field names: ");
  for (const auto& info : former_field_names)
    string.append(info.first).append(" -> ").append(info.second).append(",");
  string.pop_back();
  string.append("\n - Later field names: ");
  for (const auto& info : later_field_names)
    string.append(info.first).append(" -> ").append(info.second).append(",");
  string.pop_back();
  spdlog::info(
      string, former_path, later_path, later_raster_name_field_name, masks_dir,
      mask_reg_suffix);
  auto esri_shapefile_driver(
      GetGDALDriverManager()->GetDriverByName("Esri Shapefile"));
  GDALDatasetUniquePtr
      former_dataset(GDALDataset::Open(
          former_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY)),
      later_dataset(GDALDataset::Open(
          later_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY)),
      output_dataset(esri_shapefile_driver->Create(
          output_path.c_str(), 0, 0, 00, GDT_Unknown, nullptr));
  if (!former_dataset) {
    spdlog::warn("Opening {} failed", former_path);
    return false;
  }
  if (!later_dataset) {
    spdlog::warn("Opening {} failed", later_path);
    return false;
  }
  if (!output_dataset) {
    spdlog::warn("Creating {} failed", output_path);
    return false;
  }
  auto former_layer(former_dataset->GetLayer(0)),
      later_layer(later_dataset->GetLayer(0)),
      output_layer(output_dataset->CreateLayer(
          "", later_layer->GetSpatialRef(), wkbPolygon));

  spdlog::debug("Creating fields for the output composite table");
  // Copy fields from the later composite table
  for (const auto& info : later_field_names) {
    int idx(later_layer->GetLayerDefn()->GetFieldIndex(info.first.c_str()));
    if (idx == -1) {
      spdlog::warn(
          "The later composite table doesn't contain the \"{}\" field name",
          info.first);
      return false;
    }
    OGRFieldDefn field(later_layer->GetLayerDefn()->GetFieldDefn(idx));
    field.SetName(info.second.c_str());
    output_layer->CreateField(&field);
    spdlog::info(
        "Copying the field {} -> {} from the later composite table",
        info.first, info.second);
  }

  // Copy fields from the former composite table
  for (const auto& info : former_field_names) {
    int idx(former_layer->GetLayerDefn()->GetFieldIndex(info.first.c_str()));
    if (idx == -1) {
      spdlog::warn(
          "The former composite table doesn't contain the \"{}\" field name",
          info.first);
      return false;
    }
    OGRFieldDefn field(former_layer->GetLayerDefn()->GetFieldDefn(idx));
    field.SetName(info.second.c_str());
    output_layer->CreateField(&field);
    spdlog::info(
        "Creating the field {} -> {} from the former composite table",
        info.first, info.second);
  }
  spdlog::info("Creating fields for the output composite table - done");

  // Create features for all the intersection situation
  spdlog::info("Creating features for all the intersection situation");
  for (int i = 0; i < later_layer->GetFeatureCount(); i++) {
    auto later_feature(later_layer->GetFeature(i));
    auto later_geometry(later_feature->GetGeometryRef());
    former_layer->SetSpatialFilter(later_geometry);
    for (const auto& former_feature : former_layer) {
      OGRFeatureUniquePtr output_feature(OGRFeature::CreateFeature(
          output_layer->GetLayerDefn()));
      output_feature->SetGeometryDirectly(
          later_geometry->Intersection(former_feature->GetGeometryRef()));
      for (const auto& info : later_field_names)
        output_feature->SetField(
            info.second.c_str(),
            later_feature->GetRawFieldRef(later_feature->GetFieldIndex(
                info.first.c_str())));
      for (const auto& info : former_field_names)
        output_feature->SetField(
            info.second.c_str(),
            former_feature->GetRawFieldRef(former_feature->GetFieldIndex(
                info.first.c_str())));
      output_layer->CreateFeature(output_feature.get());
    }
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1,
        later_layer->GetFeatureCount());
  }
  former_layer->SetSpatialFilter(nullptr);
  spdlog::info("Creating features for all the intersection situation - done");

  // Remove rasters' mask from the output composite table
  spdlog::info("Removing rasters' mask from the output composite table");
  for (int i = 0; i < output_layer->GetFeatureCount(); i++) {
    auto feature(output_layer->GetFeature(i));
    OGRGeometryUniquePtr geometry(feature->GetGeometryRef()->clone());
    fs::path raster_name(feature->GetFieldAsString(
        later_raster_name_field_name.c_str()));
    std::regex reg(raster_name.stem().string() + mask_reg_suffix);
    std::string mask_path;
    for (const auto& entry : fs::directory_iterator(masks_dir))
      if (std::regex_match(entry.path().filename().string(), reg)) {
        mask_path = entry.path().string();
        break;
      }
    GDALDatasetUniquePtr mask_dataset(GDALDataset::Open(
        mask_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
    if (!mask_dataset) continue;
    auto mask_layer(mask_dataset->GetLayer(0));
    mask_layer->SetSpatialFilter(feature->GetGeometryRef());
    for (const auto& mask_feature : mask_layer)
      geometry.reset(geometry->Difference(mask_feature->GetGeometryRef()));
    feature->SetGeometryDirectly(geometry.release());
    output_layer->SetFeature(feature);
    if ((i + 1) % 10 == 0 && (i + 1) != output_layer->GetFeatureCount())
      spdlog::info(
          "----------- {}/{} - done ----------", i + 1,
          output_layer->GetFeatureCount());
  }
  spdlog::info(
      "----------- {}/{} - done ----------", output_layer->GetFeatureCount(),
      output_layer->GetFeatureCount());
  spdlog::info(
      "Removing rasters' mask from the output composite table - done");
  return true;
}

int main(int argc, char* argv[]) {
  cxxopts::Options options("GuangDong2022");
  options.add_options()
      ("h,help", "Print usage")
      ("former", "The input former composite table path",
          cxxopts::value<std::string>())
      ("later", "The input later composite table path",
          cxxopts::value<std::string>())
      ("output", "The output composite table path",
          cxxopts::value<std::string>())
      ("masks-dir", "The masks' directory",
          cxxopts::value<std::string>()->default_value(""))
      ("mask-reg-suffix", "The mask file suffix for the regular expression",
          cxxopts::value<std::string>()->default_value(".*_cloud\\.shp"))
      ("log-level", "Available levels are trace(t), debug(d), info(i), "
          "warn(w), err(e), critical(c)", cxxopts::value<std::string>());
  auto result(options.parse(argc, argv));
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  if (!result.count("former") || !result.count("later") ||
      !result.count("output")) {
    std::cout << "Lack of the \"former\", the \"later\" or the \"output\" "
      << "argument" << std::endl;
    return -1;
  }
  std::regex esri_shapefile_reg(".*\\.(shp)");
  if (!std::regex_match(
          result["former"].as<std::string>().c_str(), esri_shapefile_reg)) {
    std::cout << "The former composite table path does not end with \".shp\""
    << std::endl;
    return -1;
  }
  if (!std::regex_match(
          result["later"].as<std::string>().c_str(), esri_shapefile_reg)) {
    std::cout << "The later composite table path does not end with \".shp\""
    << std::endl;
    return -1;
  }
  if (!std::regex_match(
          result["output"].as<std::string>().c_str(), esri_shapefile_reg)) {
    std::cout << "The output composite table path does not end with \".shp\""
    << std::endl;
    return -1;
  }
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
  std::ifstream i(fs::path(argv[0]).parent_path() / "config.json");
  if (!i) {
    std::cout << "Opening the configuration file failed" << std::endl;
    return -1;
  }
  auto config(nlohmann::json::parse(i));
  std::map<std::string, std::string> former_field_names, later_field_names;
  for (const auto& info : config.at("former_field_names"))
    former_field_names[info[0]] = info[1];
  for (const auto& info : config.at("later_field_names"))
    later_field_names[info[0]] = info[1];
  RunTask(
      result["former"].as<std::string>().c_str(),
      result["later"].as<std::string>().c_str(),
      result["output"].as<std::string>().c_str(), former_field_names,
      later_field_names, config.at("later_raster_name_field_name"),
      result["masks-dir"].as<std::string>(),
      result["mask-reg-suffix"].as<std::string>());
  return 0;
}