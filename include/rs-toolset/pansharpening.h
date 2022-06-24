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

/// <summary>
/// Abstract pansharpening class
/// </summary>
class RS_TOOLSET_API PansharpeningInterface {
 public:
  /// <summary>
  /// Run pansharpening on the PAN raster path and MS raster path
  /// </summary>
  /// <param name="pan_path">The PAN raster path</param>
  /// <param name="ms_path">The PAN raster path</param>
  /// <param name="pansharpened_path">The output pansharpened raster path</param>
  /// <param name="use_rpc">Whether uses the RPC information</param>
  /// <param name="use_stretch">Whether uses stretch for the pansharpened raster</param>
  /// <param name="pansharpened_bands_map">The Bands' map, default is empty means all bands</param>
  /// <returns>Running state</returns>
  virtual bool Run(
      const std::string& pan_path,
      const std::string& ms_path,
      const std::string& pansharpened_path,
      bool use_rpc = true,
      bool use_stretch = true,
      const std::vector<int>& pansharpened_bands_map = {}) = 0;
};

/// <summary>
/// Gram-Schmidt pansharpening class
/// </summary>
class RS_TOOLSET_API GramSchmidt : virtual public PansharpeningInterface {
 public:
  /// <summary>
  /// Create the Gram-Schmidt pansharpening shared pointer
  /// </summary>
  /// <param name="block_size">The block size per operation</param>
  /// <returns>The output Gram-Schmidt pansharpening shared pointer</returns>
  static std::shared_ptr<GramSchmidt> Create(int block_size = 16384);
};

/// <summary>
/// Gram-Schmidt adaptive pansharpening class
/// </summary>
class RS_TOOLSET_API GramSchmidtAdaptive 
    : virtual public PansharpeningInterface {
 public:
  /// <summary>
  /// Create the Gram-Schmidt adaptive pansharpening shared pointer
  /// </summary>
  /// <param name="block_size">The block size per operation</param>
  /// <returns>The output Gram-Schmidt adaptive pansharpening shared pointer</returns>
  static std::shared_ptr<GramSchmidtAdaptive> Create(int block_size = 16384);
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_PANSHARPENING_H_