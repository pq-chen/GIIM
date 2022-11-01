LinKIn_AddTarget(
  MODE SHARED
  CXX_STANDARD 17
  INCLUDE_PRIVATE
    "include"
  LIBRARY_PRIVATE
    Eigen3::Eigen
    LinKIn::RS-Toolset_core_stretching
    LinKIn::RS-Toolset_core_utils
  COMPILE_OPTION_PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
)