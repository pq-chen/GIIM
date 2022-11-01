LinKIn_AddTarget(
  MODE EXE
  CXX_STANDARD 17
  INCLUDE_PRIVATE
    "include"
  LIBRARY_PRIVATE
    cxxopts::cxxopts
    LinKIn::RS-Toolset_core_mosaicking
    LinKIn::RS-Toolset_core_color_balancing
    LinKIn::RS-Toolset_core_utils
)