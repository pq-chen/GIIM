cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

find_package(LCMake REQUIRED)

LinKIn_InitVcpkg(
  ROOT D:/vcpkg 
  TRIPLET x64-windows)

project(RS-Toolset VERSION 1.0.0)
message(STATUS "[Project] ${PROJECT_NAME}")

LinKIn_InitProject()

find_package(cxxopts REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(GDAL REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)
find_package(spdlog REQUIRED)
find_package(suitesparse REQUIRED)

LinKIn_AddSubdirectories(src/core/utils)
LinKIn_AddSubdirectories(src/core/stretching)
LinKIn_AddSubdirectories(src/core)
LinKIn_AddSubdirectories(src/app)