#pragma once

//* Volumetric null-scattering mis bdpt

#include "camera.h"
#include "integrator.h"
#include "interaction.h"
#include "light.h"
#include "pbrt.h"
#include "reflection.h"
#include "sampling.h"
#include "scene.h"
#include <unordered_map>

namespace pbrt
{

struct VN_Vertex
{
  enum class Type
  {
    Camera,
    Light,
    Surface,
    Medium
  } type;

  Spectrum beta;
  bool     delta = false;

  static VN_Vertex createCamera(const Camera* camera, const Ray& ray, const Spectrum& beta);
  static VN_Vertex createLight(const Light* light, const Ray& ray, const Normal3f& n_light,
                               const Spectrum& Le, Float pdf);
  static VN_Vertex createMedium(const MediumInteraction& mi, const Spectrum& beta, Float pdf,
                                const VN_Vertex& prev);
  static VN_Vertex createSurface(const SurfaceInteraction& si, const Spectrum& beta, Float pdf,
                                 const VN_Vertex& prev);
};

int VN_GenerateCameraSubpath(const Scene& scene, Sampler& sampler, MemoryArena& arena, int maxDepth,
                             const Camera& camera, const Point2f& p_film,
                             VN_Vertex* camera_subpath);

int VN_GenerateLightSubpath(const Scene& scene, Sampler& sampler, MemoryArena& arena, int maxDepth,
                            Float time, const Distribution1D& light_distribution,
                            const std::unordered_map<const Light*, size_t>& light_to_index,
                            VN_Vertex*                                      light_subpath);

Spectrum VN_ConnectBDPT(const Scene& scene, VN_Vertex* light_subpath, VN_Vertex* camera_subpath,
                        int n_lv, int n_cv, const Distribution1D& light_distr,
                        const std::unordered_map<const Light*, size_t>& light_to_index,
                        const Camera& camera, Sampler& sampler, Point2f* p_raster,
                        Float* misw = nullptr);

int VN_RandomWalk(const Scene& scene, RayDifferential ray, Sampler& sampler, MemoryArena& arena,
                  Spectrum beta, Float pdf, int maxDepth, TransportMode mode, VN_Vertex* subpath);

class VNBDPTIntegrator : public Integrator
{
  public:
  VNBDPTIntegrator(std::shared_ptr<Sampler> sampler, std::shared_ptr<const Camera> camera,
                   int maxDepth, bool visualizeStrategies, bool visualizeWeights,
                   const Bounds2i& pixelBounds, const std::string& lightSampleStrategy = "power");

  void Render(const Scene& scene);

  private:
  std::shared_ptr<Sampler>      sampler;
  std::shared_ptr<const Camera> camera;

  const int         maxDepth;
  const bool        visualizeStrategies;
  const bool        visualizeWeights;
  const Bounds2i    pixelBounds;
  const std::string lightSampleStrategy;
};
} // namespace pbrt