#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include <cxxopts.hpp>
#include <gdal_priv.h>
#include <nlohmann/json.hpp>
#include <ogrsf_frmts.h>
#include <spdlog/spdlog.h>

#include <rs-toolset/utils.hpp>

namespace fs = std::filesystem;
using namespace rs_toolset;

bool Run2022GuangdongJcsTaskCheck(
    const std::vector<std::string>& later_paths,
    const std::string& former_path,
    const std::string& output_path,
    const std::vector<std::pair<std::string, std::string>>& later_to_output,
    const std::vector<std::pair<std::string, std::string>>& former_to_output,
    const std::string& later_raster_name_field_name,
    const std::string& masks_dir,
    const std::string& mask_suffix,
    const std::vector<std::string>& border_paths,
    double threshold,
    const std::string& encoding,
    std::vector<GDALDatasetUniquePtr>& later_datasets,
    GDALDatasetUniquePtr& former_dataset,
    GDALDatasetUniquePtr& output_dataset,
    std::vector<int>& later_fields_idx,
    std::vector<int>& former_fields_idx,
    int& later_raster_name_idx,
    std::vector<std::pair<OGRGeometryUniquePtr, std::string>>& border_infos) {
  for (const auto& path : later_paths) {
    if (later_datasets.emplace_back(GDALDataset::Open(
            path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
        !later_datasets.back()) {
      spdlog::error("Opening {} failed", path);
      return false;
    }
  }
  if (former_dataset.reset(GDALDataset::Open(
          former_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
      !former_dataset) {
    spdlog::error("Opening {} failed", former_path);
    return false;
  }
  if (auto driver(utils::GetVectorDriverByPath(output_path)); !driver) {
    spdlog::error("Finding the driver for {} failed", output_path);
    return false;
  } else if (output_dataset.reset(driver->Create(
          output_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
      !output_dataset) {
    spdlog::error("Creating {} failed", output_path);
    return false;
  }
  auto later_layer(later_datasets[0]->GetLayer(0));
  for (const auto& info : later_to_output) {
    if (later_fields_idx.push_back(later_layer->GetLayerDefn()->GetFieldIndex(
            utils::CreateCplString(CPLRecode(
                info.first.c_str(), "CP936", "UTF-8")).get()));
        later_fields_idx.back() == -1) {
      spdlog::error(
          "The \"later_paths\"[0] {} does not contain field name {}",
          later_paths[0], info.first);
      return false;
    } else if (later_raster_name_idx == -1 &&
        info.second == later_raster_name_field_name) {
      later_raster_name_idx = static_cast<int>(later_fields_idx.size() - 1);
    }
  }
  if (later_raster_name_idx == -1) {
    spdlog::error(
        "The \"later_paths\"[0] {} does not contain the field name {}",
        later_paths[0], later_raster_name_field_name);
    return false;
  }
  auto former_layer(former_dataset->GetLayer(0));
  for (const auto& info : former_to_output) {
    if (former_fields_idx.push_back(former_layer->GetLayerDefn()->GetFieldIndex(
            utils::CreateCplString(CPLRecode(
                info.first.c_str(), "CP936", "UTF-8")).get()));
        former_fields_idx.back() == -1) {
      spdlog::error(
          "The \"former_path\" {} does not contain the field name {}",
          former_path, info.first);
      return false;
    }
  }
  if (!masks_dir.empty() &&
      (!fs::is_directory(masks_dir) || !fs::exists(masks_dir))) {
    spdlog::error(
        "The \"masks_dir\" {} is not a directory or does not exist", masks_dir);
    return false;
  }
  for (const auto& path : border_paths) {
    if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
            path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY)); !dataset) {
      spdlog::error("Opening {} failed", path);
      return false;
    } else {
      auto layer(dataset->GetLayer(0));
      for (const auto& feature : layer) {
        border_infos.emplace_back(
            feature->StealGeometry(), feature->GetFieldAsString(0));
      }
    }
  }
  if (threshold < 0) {
    spdlog::error("The \"threshold\" {} must be non-negative", threshold);
    return false;
  }
  return true;
}

GDALDatasetUniquePtr CombineLaterVectors(
    const std::vector<GDALDatasetUniquePtr>& later_datasets) {
  GDALDatasetUniquePtr output_dataset(
      GetGDALDriverManager()->GetDriverByName("Memory")->CreateCopy(
          "", later_datasets[0].get(), true, nullptr, nullptr, nullptr));
  auto output_layer(output_dataset->GetLayer(0));
  for (auto i(1); i < later_datasets.size(); ++i) {
    OGRGeometryUniquePtr covered_border(new OGRMultiPolygon);
    auto _covered_border(covered_border->toMultiPolygon());
    auto later_layer(later_datasets[i]->GetLayer(0));
    for (const auto& later_feature : later_layer) {
      if (OGRGeometryUniquePtr later_geometry(later_feature->StealGeometry());
          !later_geometry->Intersect(_covered_border)) {
        _covered_border->addGeometryDirectly(later_geometry.release());
      } else {
        for (auto j(_covered_border->getNumGeometries() - 1); j >= 0; --j) {
          if (auto covered_polygon(_covered_border->getGeometryRef(j));
              later_geometry->Intersect(covered_polygon)) {
            later_geometry.reset(later_geometry->Union(covered_polygon));
            _covered_border->removeGeometry(j);
          }
        }
        switch (later_geometry->getGeometryType()) {
          case wkbPolygon: {
            auto _later_geometry(later_geometry->toPolygon());
            for (auto j(_later_geometry->getNumInteriorRings() - 1);
                j >= 0; --j) {
              if (_later_geometry->getInteriorRing(j)->get_Area() == 0.0)
                _later_geometry->removeRing(j + 1);
            }
            break;
          }
          case wkbMultiPolygon: {
            utils::ExtractBiggestPolygon(later_geometry);
          }
        }
        _covered_border->addGeometryDirectly(later_geometry.release());
      }
    }
    for (const auto& covered_polygon : _covered_border) {
      for (auto j(output_layer->GetFeatureCount() - 1); j >= 0; --j) {
        OGRFeatureUniquePtr output_feature(output_layer->GetFeature(j));
        if (auto output_geometry(output_feature->GetGeometryRef());
            covered_polygon->Intersect(output_geometry)) {
          output_feature->SetGeometryDirectly(output_geometry->Difference(
              covered_polygon));
          output_layer->SetFeature(output_feature.get());
        } else if (covered_polygon->Contains(output_geometry)) {
          auto last_fid(output_layer->GetFeatureCount() - 1);
          output_layer->DeleteFeature(j);
          if (j != last_fid) {
            output_feature.reset(output_layer->GetFeature(last_fid));
            output_feature->SetFID(j);
            output_layer->DeleteFeature(last_fid);
            output_layer->SetFeature(output_feature.get());
          }
        }
      }
    }
    auto output_layer_defn(output_layer->GetLayerDefn());
    for (const auto& later_feature : later_layer) {
      OGRFeatureUniquePtr output_feature(OGRFeature::CreateFeature(
          output_layer_defn));
      output_feature->SetFID(output_layer->GetFeatureCount());
      output_feature->SetGeometryDirectly(later_feature->StealGeometry());
      for (auto j(0); j < output_layer_defn->GetFieldCount(); ++j) {
        if (auto idx(later_layer->GetLayerDefn()->GetFieldIndex(
                output_layer_defn->GetFieldDefn(j)->GetNameRef())); idx != -1) {
          output_feature->SetField(j, later_feature->GetRawFieldRef(idx));
        }
      }
      output_layer->CreateFeature(output_feature.get());
    }
  }
  return output_dataset;
}

void CreateFields(
    OGRLayer* later_layer,
    OGRLayer* former_layer,
    const std::vector<std::pair<std::string, std::string>>& later_to_output,
    const std::vector<std::pair<std::string, std::string>>& former_to_output,
    const std::vector<int>& later_fields_idx,
    const std::vector<int>& former_fields_idx,
    OGRLayer* output_layer) {
  for (auto i(0); i < later_fields_idx.size(); ++i) {
    OGRFieldDefn field(later_layer->GetLayerDefn()->GetFieldDefn(
        later_fields_idx[i]));
    field.SetName(utils::CreateCplString(CPLRecode(
        later_to_output[i].second.c_str(), "CP936", "UTF-8")).get());
    output_layer->CreateField(&field);
    spdlog::info(
        "Copying the field {} -> {} from the later layer",
        later_to_output[i].first, later_to_output[i].second);
  }
  for (auto i(0); i < former_fields_idx.size(); ++i) {
    OGRFieldDefn field(former_layer->GetLayerDefn()->GetFieldDefn(
        former_fields_idx[i]));
    field.SetName(utils::CreateCplString(CPLRecode(
        former_to_output[i].second.c_str(), "CP936", "UTF-8")).get());
    output_layer->CreateField(&field);
    spdlog::info(
        "Creating the field {} -> {} from the former layer",
        former_to_output[i].first, former_to_output[i].second);
  }
}

void CreateFeatures(
    OGRLayer* later_layer,
    OGRLayer* former_layer,
    const std::vector<int>& later_fields_idx,
    const std::vector<int>& former_fields_idx,
    OGRLayer* output_layer) {
  spdlog::info("Creating output features in all intersecting cases");
  for (auto i(0); i < later_layer->GetFeatureCount(); ++i) {
    OGRFeatureUniquePtr later_feature(later_layer->GetFeature(i));
    auto later_geometry(later_feature->GetGeometryRef());
    former_layer->SetSpatialFilter(later_geometry);
    for (const auto& former_feature : former_layer) {
      std::vector<OGRGeometry*> polygons;
      switch (OGRGeometryUniquePtr geometry(later_geometry->Intersection(
              former_feature->GetGeometryRef())); geometry->getGeometryType()) {
        case wkbPolygon: {
          polygons.push_back(geometry.release());
          break;
        }
        case wkbMultiPolygon: {
          for (const auto& polygon : geometry->toMultiPolygon())
            polygons.push_back(polygon->clone());
        }
      }
      for (auto& polygon : polygons) {
        OGRFeatureUniquePtr output_feature(OGRFeature::CreateFeature(
            output_layer->GetLayerDefn()));
        output_feature->SetFID(output_layer->GetFeatureCount());
        output_feature->SetGeometryDirectly(polygon);
        for (auto j(0); j < later_fields_idx.size(); ++j) {
          output_feature->SetField(
              j, later_feature->GetRawFieldRef(later_fields_idx[j]));
        }
        for (auto j(0); j < former_fields_idx.size(); ++j) {
          output_feature->SetField(
              j + static_cast<int>(later_fields_idx.size()),
              former_feature->GetRawFieldRef(former_fields_idx[j]));
        }
        output_layer->CreateFeature(output_feature.get());
      }
    }
    spdlog::info(
        "---------- {}/{} - done ----------", i + 1,
        later_layer->GetFeatureCount());
  }
  former_layer->SetSpatialFilter(nullptr);
  spdlog::info("Creating output features in all intersecting cases - done");
}

std::vector<std::string> ExportMasksPath(
    const std::string& masks_dir,
    const std::string& mask_suffix) {
  std::vector<std::string> masks_path;
  std::regex reg(".*" + mask_suffix);
  for (const auto& entry : fs::directory_iterator(masks_dir)) {
    if (auto path(entry.path().string()); std::regex_match(path, reg)) {
      if (GDALDatasetUniquePtr dataset(GDALDataset::Open(
              path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY)); dataset) {
        masks_path.push_back(path);
      }
    }
  }
  return masks_path;
}

void RemoveMasks(
    const std::string& mask_suffix,
    const std::vector<std::string>& masks_path,
    int later_raster_name_idx,
    double threshold,
    OGRLayer* output_layer,
    std::vector<std::pair<OGRGeometryUniquePtr, int>>& mask_infos) {
  spdlog::info("Removing masks for the output layer");
  auto count(static_cast<int>(output_layer->GetFeatureCount()));
  for (auto i(count - 1); i >= 0; --i) {
    OGRFeatureUniquePtr output_feature(output_layer->GetFeature(i));
    auto mask_name(
        fs::path(utils::CreateCplString(CPLRecode(
            output_feature->GetFieldAsString(later_raster_name_idx), "UTF-8",
            "CP936")).get()).stem().string() + mask_suffix);
    auto it(std::find_if(
        masks_path.begin(), masks_path.end(),
        [&](const std::string& s) {
          return fs::path(s).filename() == mask_name;
        }));
    if (it == masks_path.end()) {
      spdlog::info(
          "Skipping removing masks "
          "since the corresponding mask file not found");
      continue;
    }
    GDALDatasetUniquePtr mask_dataset(GDALDataset::Open(
        it->c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
    auto mask_layer(mask_dataset->GetLayer(0));
    OGRGeometryUniquePtr output_geometry(output_feature->StealGeometry());
    mask_layer->SetSpatialFilter(output_geometry.get());
    auto remained(true);
    for (const auto& mask_feature : mask_layer) {
      if (auto mask_geometry(mask_feature->GetGeometryRef());
          mask_geometry->Contains(output_geometry.get())) {
        auto last_fid(output_layer->GetFeatureCount() - 1);
        output_layer->DeleteFeature(i);
        if (i != last_fid) {
          output_feature.reset(output_layer->GetFeature(last_fid));
          output_feature->SetFID(i);
          output_layer->DeleteFeature(last_fid);
          output_layer->SetFeature(output_feature.get());
        }
        remained = false;
        break;
      } else {
        OGRGeometryUniquePtr intersecting_geometry(
            output_geometry->Intersection(mask_geometry));
        auto area(0.0);
        switch (intersecting_geometry->getGeometryType()) {
          case wkbPolygon: {
            area = intersecting_geometry->toPolygon()->get_Area();
            break;
          }
          case wkbMultiPolygon: {
            area = intersecting_geometry->toMultiPolygon()->get_Area();
          }
        }
        if (area > threshold) {
          output_geometry.reset(output_geometry->Difference(mask_geometry));
          mask_infos.emplace_back(intersecting_geometry.release(), i);
        }
      }
    }
    if (remained) {
      output_feature->SetGeometryDirectly(output_geometry.release());
      output_layer->SetFeature(output_feature.get());
    }
    spdlog::info("---------- {}/{} - done ----------", count - i, count);
  }
  spdlog::info("Removing masks for the output layer - done");
}

void FillMasks(
    const std::string& mask_suffix,
    const std::vector<std::string>& masks_path,
    const std::vector<std::pair<OGRGeometryUniquePtr, std::string>>&
        border_infos,
    int later_raster_name_idx,
    OGRLayer* output_layer,
    std::vector<std::pair<OGRGeometryUniquePtr, int>>& mask_infos) {
  spdlog::info("Filling masks for the output layer");
  auto count(0);
  for (auto& mask_info : mask_infos) {
    for (const auto& border_info : border_infos) {
      if (!border_info.first->Contains(mask_info.first.get())) continue;
      auto mask_name(
          fs::path(utils::CreateCplString(CPLRecode(
              border_info.second.c_str(), "UTF-8", "CP936")).get())
                  .stem().string() + mask_suffix);
      auto it(std::find_if(
          masks_path.begin(), masks_path.end(),
          [&](const std::string& s) {
            return fs::path(s).filename() == mask_name;
          }));
      if (it == masks_path.end()) continue;
      GDALDatasetUniquePtr mask_dataset(GDALDataset::Open(
          it->c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY));
      auto mask_layer(mask_dataset->GetLayer(0));
      mask_layer->SetSpatialFilter(mask_info.first.get());
      if (mask_layer->GetFeatureCount() == 0) {
        OGRFeatureUniquePtr output_feature(output_layer->GetFeature(
            mask_info.second));
        output_feature->SetFID(output_layer->GetFeatureCount());
        output_feature->SetGeometryDirectly(mask_info.first.release());
        output_feature->SetField(
            later_raster_name_idx, border_info.second.c_str());
        output_layer->CreateFeature(output_feature.get());
        ++count;
        break;
      }
    }
  }
  spdlog::info("Filling {} masks for the output layer - done", count);
}

bool Run2022GuangdongJcsTask(
    const std::vector<std::string>& later_paths,
    const std::string& former_path,
    const std::string& output_path,
    const std::vector<std::pair<std::string, std::string>>& later_to_output,
    const std::vector<std::pair<std::string, std::string>>& former_to_output,
    const std::string& later_raster_name_field_name,
    const std::string& masks_dir,
    const std::string& mask_suffix,
    const std::vector<std::string>& border_paths,
    double threshold,
    const std::string& encoding) {
  std::vector<GDALDatasetUniquePtr> later_datasets;
  GDALDatasetUniquePtr former_dataset(nullptr), output_dataset(nullptr);
  std::vector<int> later_fields_idx, former_fields_idx;
  auto later_raster_name_idx(-1);
  std::vector<std::pair<OGRGeometryUniquePtr, std::string>> border_infos;
  if (!Run2022GuangdongJcsTaskCheck(
          later_paths, former_path, output_path, later_to_output,
          former_to_output, later_raster_name_field_name, masks_dir,
          mask_suffix, border_paths, threshold, encoding, later_datasets,
          former_dataset, output_dataset, later_fields_idx, former_fields_idx,
          later_raster_name_idx, border_infos)) {
    return false;
  }
  std::string string(
      "Running a 2022 GuangDong JCS task to {} with\n"
      " - Later vector(s') path: \n");
  for (const auto& path : later_paths)
    string.append(path).append("\n");
  string.append(
      " - Former vector path: {}\n"
      " - Later raster name field name: {}\n"
      " - Masks' directory: {}\n"
      " - Mask suffix: {}\n"
      " - Border(s') path: \n");
  for (const auto& path : border_paths)
    string.append(path).append("\n");
  string.append(" - Threshold: {}\n - Encoding: {}");
  spdlog::info(
      string, output_path, former_path, later_raster_name_field_name, masks_dir,
      mask_suffix, threshold, encoding);

  GDALDatasetUniquePtr later_dataset(nullptr);
  OGRLayer* later_layer(nullptr);
  if (later_datasets.size() == 1) {
    later_layer = later_datasets[0]->GetLayer(0);
  } else {
    later_dataset = CombineLaterVectors(later_datasets);
    later_layer = later_dataset->GetLayer(0);
  }
  std::unique_ptr<char*, void(*)(char**)> options(nullptr, nullptr);
  if (!output_path.empty()) {
    options = utils::CreateCslStringList(CSLSetNameValue(
        nullptr, "ENCODING", encoding.c_str()));
  }
  auto former_layer(former_dataset->GetLayer(0)),
      output_layer(output_dataset->CreateLayer(
          "", later_layer->GetSpatialRef(), wkbPolygon, options.get()));
  CreateFields(
      later_layer, former_layer, later_to_output, former_to_output,
      later_fields_idx, former_fields_idx, output_layer);
  CreateFeatures(
      later_layer, former_layer, later_fields_idx, former_fields_idx,
      output_layer);
  if (!masks_dir.empty()) {
    auto masks_path(ExportMasksPath(masks_dir, mask_suffix));
    std::vector<std::pair<OGRGeometryUniquePtr, int>> mask_infos;
    RemoveMasks(
        mask_suffix, masks_path, later_raster_name_idx, threshold, output_layer,
        mask_infos);
    FillMasks(
        mask_suffix, masks_path, border_infos, later_raster_name_idx,
        output_layer, mask_infos);
  }
  spdlog::info("Running a 2022 GuangDong JCS task - done");
  return true;
}

int main(int argc, char* argv[]) {
  cxxopts::Options options("RS-Toolset_app_2022_guangdong_jcs");
  options.add_options()
      ("h,help", "Print usage")
      ("later", "The input later vectors' path", cxxopts::value<std::vector<std::string>>())
      ("former", "The input former vector path", cxxopts::value<std::string>())
      ("output", "The output vector path", cxxopts::value<std::string>())
      ("masks-dir", "The masks' directory", cxxopts::value<std::string>())
      ("mask-suffix", "The mask file suffix", cxxopts::value<std::string>()->default_value("_cloud.shp"))
      ("border", "The border vectors' path", cxxopts::value<std::vector<std::string>>())
      ("threshold", "The area threshold for masks", cxxopts::value<double>()->default_value("1e-4"))
      ("encoding", "Available encoding methods are \"UTF-8\" and \"CP936\"", cxxopts::value<std::string>()->default_value("UTF-8"))
      ("log-level", "Available log levels are \"trace(t)\", \"debug(d)\", \"info(i)\", \"warn(w)\", \"err(e)\", \"critical(c)\" and \"off(o)\"", cxxopts::value<std::string>()->default_value("i"));
  auto result(options.parse(argc, argv));

  utils::InitGdal(argv[0]);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  if (!result.count("former") || !result.count("later") ||
      !result.count("output")) {
    std::cout << "Lack of the \"former\", the \"later\" or the \"output\""
        << std::endl;
    return -1;
  }
  std::string masks_dir;
  if (result.count("masks-dir")) {
    if (masks_dir = result["masks-dir"].as<std::string>(); masks_dir.empty()) {
      std::cout << "The \"masks_dir\" is empty" << std::endl;
      return -1;
    }
  }
  std::vector<std::string> borders_path;
  if (result.count("border")) {
    if (borders_path = result["border"].as<std::vector<std::string>>();
        borders_path.empty()) {
      std::cout << "The \"border\" is empty" << std::endl;
      return -1;
    }
  }
  auto encoding(result["encoding"].as<std::string>());
  if (encoding != "UTF-8" && encoding != "CP936") {
    std::cout << "Available encoding methods are \"UTF-8\" and \"CP936\""
        << std::endl;
    return -1;
  }
  auto log_level(result["log-level"].as<std::string>());
  if (log_level == "t" || log_level == "trace") {
    utils::InitSpdlog("2022 GuangDong JCS", spdlog::level::level_enum::trace);
  } else if (log_level == "d" || log_level == "debug") {
    utils::InitSpdlog("2022 GuangDong JCS", spdlog::level::level_enum::debug);
  } else if (log_level == "i" || log_level == "info") {
    utils::InitSpdlog("2022 GuangDong JCS", spdlog::level::level_enum::info);
  } else if (log_level == "w" || log_level == "warn") {
    utils::InitSpdlog("2022 GuangDong JCS", spdlog::level::level_enum::warn);
  } else if (log_level == "e" || log_level == "err") {
    utils::InitSpdlog("2022 GuangDong JCS", spdlog::level::level_enum::err);
  } else if (log_level == "c" || log_level == "critical") {
    utils::InitSpdlog(
        "2022 GuangDong JCS", spdlog::level::level_enum::critical);
  } else if (log_level == "o" || log_level == "off") {
    utils::InitSpdlog("2022 GuangDong JCS", spdlog::level::level_enum::off);
  } else {
    std::cout << "Available log levels are \"trace(t)\", \"debug(d)\", "
        "\"info(i)\", \"warn(w)\", \"err(e)\", \"critical(c)\" and \"off(o)\""
        << std::endl;
    return -1;
  }

  auto config_path((fs::path(argv[0]).parent_path() / "config.json").string());
  std::ifstream config_file(config_path);
  if (!config_file) {
    std::cout << "Opening " << config_path << " failed" << std::endl;
    return -1;
  }
  auto config(nlohmann::json::parse(config_file));
  std::vector<std::pair<std::string, std::string>> 
      later_to_output, former_to_output;
  for (const auto& info : config.at("later_to_output")) {
    later_to_output.emplace_back(
        utils::CreateCplString(CPLRecode(
            info[0].get<std::string>().c_str(), "UTF-8", "CP936")).get(),
        utils::CreateCplString(CPLRecode(
            info[1].get<std::string>().c_str(), "UTF-8", "CP936")).get());
  }
  for (const auto& info : config.at("former_to_output")) {
    former_to_output.emplace_back(
        utils::CreateCplString(CPLRecode(
            info[0].get<std::string>().c_str(), "UTF-8", "CP936")).get(),
        utils::CreateCplString(CPLRecode(
            info[1].get<std::string>().c_str(), "UTF-8", "CP936")).get());
  }
  return Run2022GuangdongJcsTask(
      result["later"].as<std::vector<std::string>>(),
      result["former"].as<std::string>(), result["output"].as<std::string>(),
      later_to_output, former_to_output,
      utils::CreateCplString(CPLRecode(
          config.at("later_raster_name_field_name").get<std::string>().c_str(),
          "UTF-8", "CP936")).get(), masks_dir,
      result["mask-suffix"].as<std::string>(), borders_path,
      result["threshold"].as<double>(), encoding) ? 0 : -1;
}