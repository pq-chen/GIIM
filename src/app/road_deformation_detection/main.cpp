#include "road_deformation_detection.h"

#include <rs-toolset/utils.hpp>


using namespace rs_toolset;

int main(int argc, char* argv[]) {
  utils::InitGdal(argv[0]);

  auto rdd(road_deformation_detection::RoadDeformationDetection::Create(
      R"(D:\Data\road_deformation_detection_test\line.shp)", "",
      R"(D:\Data\road_deformation_detection_test\n22_e114_1arc_v3.tif)",
      0.5 , -0.01, 1));
  rdd->Run(
      R"(D:\Data\road_deformation_detection_test\test.tif)",
      R"(D:\Data\road_deformation_detection_test\SV1-01_20171030_L2A0000192764_1109170051408_01-PAN.rpb)",
      R"(C:\Users\LinKInLeLe43\Desktop\1_3.shp)");
}