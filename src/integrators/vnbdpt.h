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

namespace pbrt {

struct VN_EndPointInteraction : Interaction {
  union {
    const Camera *camera;
    const Light  *light;
  };

  VN_EndPointInteraction() : Interaction() {}

  VN_EndPointInteraction(const Interaction &it, const Camera *camera)
      : Interaction(it), camera(camera) {}

  VN_EndPointInteraction(const Camera *camera, const Ray &ray)
      : Interaction(ray.o, ray.time, ray.medium), camera(camera) {}

  VN_EndPointInteraction(const Light *light, const Ray &ray, const Normal3f &nl)
      : Interaction(ray.o, ray.time, ray.medium), light(light) {
    n = nl;
  }

  VN_EndPointInteraction(const Interaction &it, const Light *light)
      : Interaction(it), light(light) {}

  VN_EndPointInteraction(const Ray &ray)
      : Interaction(ray(1), ray.time, ray.medium), light(nullptr) {
    n = Normal3f(-ray.d);
  }
};

struct VN_Vertex {

  //* Static function
  static VN_Vertex createCamera(const Camera *camera, const Ray &ray, const Spectrum &beta);
  static VN_Vertex createCamera(const Camera *camera, const Interaction &it, const Spectrum &beta);
  static VN_Vertex createLight(const Light *light, const Ray &ray, const Normal3f &n_light,
                               const Spectrum &Le, Float pdf);
  static VN_Vertex createLight(const VN_EndPointInteraction &ei, const Spectrum &beta, Float pdf);
  static VN_Vertex createMedium(const MediumInteraction &mi, const Spectrum &beta,
                                const VN_Vertex &prev);
  static VN_Vertex createSurface(const SurfaceInteraction &si, const Spectrum &beta,
                                 const VN_Vertex &prev);

  //* Member function
  bool               IsLight() const;
  bool               IsDeltaLight() const;
  bool               IsConnectible() const;
  bool               IsOnSurface() const;
  const Interaction &GetInteraction() const;
  Spectrum           f(const VN_Vertex &next, TransportMode mode) const;
  Spectrum           Le(const Scene &scene, const VN_Vertex *towards) const;
  const Normal3f    &ns() const;
  const Normal3f    &ng() const;
  Point3f            p() const;

  //* Data field
  enum class Type { Camera, Light, Surface, Medium } type;
  Spectrum beta;

  union {
    VN_EndPointInteraction ei;
    MediumInteraction      mi;
    SurfaceInteraction     si;
  };

  bool delta = false;
  // The pdf of sampling towards this vertex by forward strategy in solid angle measure
  Float pdf_forward;
  // The pdf of sampling towards this vertex by inverse strategy in solid angle measure
  Float pdf_inverse;

  VN_Vertex() : ei() {}
  VN_Vertex(Type type, const VN_EndPointInteraction &ei, const Spectrum &beta)
      : type(type), beta(beta), ei(ei) {}
  VN_Vertex(const SurfaceInteraction &si, const Spectrum &beta)
      : type(Type::Surface), beta(beta), si(si) {}

  VN_Vertex(const VN_Vertex &v) { memcpy(this, &v, sizeof(VN_Vertex)); }
  VN_Vertex &operator=(const VN_Vertex &v) {
    memcpy(this, &v, sizeof(VN_Vertex));
    return *this;
  }
};

// Return the pdf in path space measure
Float    PDFDensityConvert(const VN_Vertex *from, const VN_Vertex *to, Float pdf_dir);
Float    PDF(const VN_Vertex *V_minus, const VN_Vertex *V, const VN_Vertex *V_plus,
             const Scene &scene);
Float    PDFLightOrigin(const Scene &scene, const VN_Vertex *shading_p, const VN_Vertex *light_v,
                        const Distribution1D                            &lightDistr,
                        const std::unordered_map<const Light *, size_t> &lightToDistrIndex);
Float    PDFLight(const VN_Vertex *light, const VN_Vertex *v, const Scene &scene);
Spectrum VN_G(const Scene &scene, Sampler &sampler, const VN_Vertex *v0, const VN_Vertex *v1);

int VN_GenerateCameraSubpath(const Scene &scene, Sampler &sampler, MemoryArena &arena, int maxDepth,
                             const Camera &camera, const Point2f &p_film,
                             VN_Vertex *camera_subpath);

int VN_GenerateLightSubpath(const Scene &scene, Sampler &sampler, MemoryArena &arena, int maxDepth,
                            Float time, const Distribution1D &light_distribution,
                            const std::unordered_map<const Light *, size_t> &light_to_index,
                            VN_Vertex                                       *light_subpath);

Spectrum VN_ConnectBDPT(const Scene &scene, VN_Vertex *light_subpath, VN_Vertex *camera_subpath,
                        int n_lv, int n_cv, const Distribution1D &light_distr,
                        const std::unordered_map<const Light *, size_t> &light_to_index,
                        const Camera &camera, Sampler &sampler, Point2f *p_raster,
                        Float *misw = nullptr);

int VN_RandomWalk(const Scene &scene, RayDifferential ray, Sampler &sampler, MemoryArena &arena,
                  Spectrum beta, Float pdf, int maxDepth, TransportMode mode, VN_Vertex *subpath);

Float VN_BDPTMIS(const Scene &scene, VN_Vertex *light_subpath, VN_Vertex *camera_subpath,
                 const VN_Vertex &sampled, int n_lv, int n_cv, const Distribution1D &lightPdf,
                 const std::unordered_map<const Light *, size_t> &lightToIndex);

class VNBDPTIntegrator : public Integrator {
public:
  VNBDPTIntegrator(std::shared_ptr<Sampler> sampler, std::shared_ptr<const Camera> camera,
                   int maxDepth, bool visualizeStrategies, bool visualizeWeights,
                   const Bounds2i &pixelBounds, const std::string &lightSampleStrategy = "power");

  void Render(const Scene &scene);

private:
  std::shared_ptr<Sampler>      sampler;
  std::shared_ptr<const Camera> camera;

  const int         maxDepth;
  const bool        visualizeStrategies;
  const bool        visualizeWeights;
  const Bounds2i    pixelBounds;
  const std::string lightSampleStrategy;
};

VNBDPTIntegrator *CreateVNBDPTIntegrator(const ParamSet &params, std::shared_ptr<Sampler> sampler,
                                         std::shared_ptr<const Camera> camera);
} // namespace pbrt