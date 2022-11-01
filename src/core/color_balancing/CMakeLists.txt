LinKIn_AddTarget(
  MODE SHARED
  CXX_STANDARD 17
  INCLUDE_PRIVATE
    "include"
    ${SUITESPARSE_INCLUDE_DIRS}
  LIBRARY_PRIVATE
    ${SUITESPARSE_LIBRARIES}
    Eigen3::Eigen
    LinKIn::RS-Toolset_core_utils
)