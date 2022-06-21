#ifndef RS_TOOLSET_INCLUDE_RS_TOOLSET_PANSHARPENING_H_
#define RS_TOOLSET_INCLUDE_RS_TOOLSET_PANSHARPENING_H_

#ifdef _WIN32

#ifdef RS_TOOLSET_EXPORTS
#define RS_TOOLSET_API __declspec(dllexport)
#else // LCMAKE_EXPORTS
#define RS_TOOLSET_API __declspec(dllimport)
#endif // LCMAKE_EXPORTS

#else // _WIN32
#define RS_TOOLSET_API
#endif // _WIN32

#include <memory>
#include <string>
#include <vector>


namespace rs_toolset {
namespace pansharpening {

class RS_TOOLSET_API PansharpeningInterface {
 public:
  virtual bool Run(
      const std::string& pan_path,
      const std::string& ms_path,
      const std::string& pansharpened_path,
      bool use_rpc = true,
      bool use_stretch = true,
      const std::vector<int>&pansharpened_bands_map = {}) = 0;
};

class RS_TOOLSET_API GramSchmidt : virtual public PansharpeningInterface {
 public:
  static std::shared_ptr<GramSchmidt> Create(int block_size = 16384);
};

class RS_TOOLSET_API GramSchmidtAdaptive 
    : virtual public PansharpeningInterface {
 public:
  static std::shared_ptr<GramSchmidtAdaptive> Create(int block_size = 16384);
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_PANSHARPENING_H_