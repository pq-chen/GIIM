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

/** @brief Abstract pansharpening class */
class RS_TOOLSET_API PansharpeningInterface {
 public:
  /**
   * @brief Run the pansharpening algorithm on the PAN raster and MS raster
   * @param[in] pan_path The PAN raster path
   * @param[in] ms_path The MS raster path
   * @param[out] pansharpened_path The output pansharpened raster path
   * @param[in] use_rpc Whether uses the RPC information
   * @param[in] use_stretch Whether uses the stretch algorithm on the pansharpened raster
   * @param[in] pansharpened_bands_map The bands' map, default is empty means all bands
   * @return Running state
  */
  virtual bool Run(
      const std::string& pan_path,
      const std::string& ms_path,
      const std::string& pansharpened_path,
      bool use_rpc = true,
      bool use_stretch = true,
      const std::vector<int>& pansharpened_bands_map = {}) = 0;
};

/** @brief Gram-Schmidt pansharpening class implementing the Gram-Schmidt pansharpening algorithm */
class RS_TOOLSET_API GramSchmidt : virtual public PansharpeningInterface {
 public:
  /**
   * @brief Create a Gram-Schmidt pansharpening shared pointer
   * @param[in] block_size The block size per operation
   * @return The output Gram-Schmidt pansharpening shared pointer
  */
  static std::shared_ptr<GramSchmidt> Create(int block_size = 16384);
};

/** @brief Gram-Schmidt adaptive pansharpening class implementing the Gram-Schmidt adaptive pansharpening algorithm */
class RS_TOOLSET_API GramSchmidtAdaptive
    : virtual public PansharpeningInterface {
 public:
  /**
   * @brief Create a Gram-Schmidt adaptive pansharpening shared pointer
   * @param[in] block_size The block size per operation
   * @return The output Gram-Schmidt adaptive pansharpening shared pointer
  */
  static std::shared_ptr<GramSchmidtAdaptive> Create(int block_size = 16384);
};

} // pansharpening
} // rs_toolset

#endif // RS_TOOLSET_INCLUDE_RS_TOOLSET_PANSHARPENING_H_