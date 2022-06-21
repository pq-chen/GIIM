LinKIn_AddTarget(
  MODE SHARED
  CXX_STANDARD 17
  INCLUDE_PRIVATE
    "include"
  INCLUDE_PUBLIC
    ${GDAL_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
  LIBRARY_PUBLIC
    ${GDAL_LIBRARIES}
    ${OpenCV_LIBS}
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
    spdlog::spdlog
)