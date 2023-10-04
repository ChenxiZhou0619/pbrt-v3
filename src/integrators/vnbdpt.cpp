#include "vnbdpt.h"
#include "film.h"
#include "filters/box.h"
#include "integrator.h"
#include "lightdistrib.h"
#include "paramset.h"
#include "progressreporter.h"
#include "sampler.h"
#include "stats.h"

namespace pbrt {

template <typename Type> class TemporalProxy {
public:
  // TemporalProxy Public Methods
  TemporalProxy(Type *target = nullptr, Type value = Type()) : target(target) {
    if (target) {
      backup  = *target;
      *target = value;
    }
  }
  ~TemporalProxy() {
    if (target) *target = backup;
  }
  TemporalProxy(const TemporalProxy &)            = delete;
  TemporalProxy &operator=(const TemporalProxy &) = delete;
  TemporalProxy &operator=(TemporalProxy &&other) {
    if (target) *target = backup;
    target       = other.target;
    backup       = other.backup;
    other.target = nullptr;
    return *this;
  }

private:
  Type *target, backup;
};

VNBDPTIntegrator::VNBDPTIntegrator(
    std::shared_ptr<Sampler> sampler, std::shared_ptr<const Camera> camera, int maxDepth,
    bool visualizeStrategies, bool visualizeWeights, const Bounds2i &pixelBounds,
    const std::vector<std::shared_ptr<VolumetricLight>> &volumetric_lights,
    const std::string                                   &lightSampleStrategy)
    : sampler(sampler), camera(camera), maxDepth(maxDepth),
      visualizeStrategies(visualizeStrategies), visualizeWeights(visualizeWeights),
      pixelBounds(pixelBounds), lightSampleStrategy(lightSampleStrategy),
      volumetric_lights(volumetric_lights) {}

int VN_RandomWalk(const Scene &scene, RayDifferential ray, Sampler &sampler, MemoryArena &arena,
                  Spectrum beta, Float pdf, int maxDepth, TransportMode mode, VN_Vertex *subpath) {
  if (maxDepth == 0) return 0;

  int   bounces         = 0;
  Float pdf_forward_dir = pdf;
  Float pdf_reverse_dir = .0f;

  while (true) {
    SurfaceInteraction isect;
    bool               found_intersection = scene.Intersect(ray, &isect);

    // Found an intersection, create vertex
    VN_Vertex *cur = subpath + bounces;
    VN_Vertex *prv = subpath + bounces - 1;

    if (ray.medium) {
      bool  scattered  = false;
      bool  terminated = false;
      Float tmax       = found_intersection ? ray.tMax : Infinity;

      RNG rng(rand());

      Spectrum T_maj = SampleT_maj(ray, tmax, rng, arena, [&](MajorantSampleRecord maj_rec) {
        // Terminate if no contribution
        if (beta.IsBlack()) {
          terminated = true;
          return false;
        }

        // TODO Handle emission

        Spectrum T_maj     = maj_rec.T_maj;
        Spectrum sigma_a   = maj_rec.sigma_a;
        Spectrum sigma_s   = maj_rec.sigma_s;
        Spectrum sigma_n   = maj_rec.sigma_n;
        Spectrum sigma_maj = sigma_a + sigma_s + sigma_n;

        // Sample a collision type and update beta
        Float pdf_absorb  = sigma_a[0] / sigma_maj[0];
        Float pdf_scatter = sigma_s[0] / sigma_maj[0];

        auto SampleScatterType = [](Float um, Float p_absorb, Float p_scatter) {
          if (um < p_absorb)
            return 0;
          else if (um < p_absorb + p_scatter)
            return 1;
          return 2;
        };

        int mode = SampleScatterType(rng.UniformFloat(), pdf_absorb, pdf_scatter);

        switch (mode) {
        case 0 /*Absorb*/:
          terminated = true;
          return false;

        case 2 /* Null scattering*/:
          beta *= T_maj * sigma_n / (T_maj[0] * sigma_n[0]);
          return true;

        case 1 /*Scatter*/:
          MediumInteraction mi(maj_rec.p, -ray.d, ray.time, ray.medium, maj_rec.phase);

          beta *= T_maj * sigma_s / (T_maj[0] * sigma_s[0]);

          *cur = VN_Vertex::createMedium(mi, beta, *prv);
          //* Convert the measure from solid angle measure to pathspace measure
          Float pdf_prv2cur_pathspace = PDFDensityConvert(prv, cur, pdf_forward_dir);
          cur->pdf_forward            = pdf_prv2cur_pathspace;

          if (++bounces == maxDepth) {
            terminated = true;
            return false;
          }

          // sample phase and spwan the path
          Vector3f wi;
          Vector3f wo     = Normalize(-ray.d);
          pdf_forward_dir = maj_rec.phase->Sample_p(wo, &wi, sampler.Get2D());

          ray             = mi.SpawnRay(wi);
          scattered       = true;
          pdf_reverse_dir = maj_rec.phase->p(wi, wo);
          //* Convert the measure from solid angle measure to pathspace measure
          Float pdf_cur2prv_pathspace = PDFDensityConvert(cur, prv, pdf_reverse_dir);
          prv->pdf_inverse            = pdf_cur2prv_pathspace;
          return false;
        }
        return false;
      });

      if (terminated || beta.IsBlack()) break;
      if (scattered) continue;

      beta *= T_maj / T_maj[0];
    }

    if (!found_intersection) break;

    isect.ComputeScatteringFunctions(ray, arena, true, mode);
    if (!isect.bsdf) {
      ray = isect.SpawnRay(ray.d);
      continue;
    }

    *cur = VN_Vertex::createSurface(isect, beta, *prv);

    //* Convert the measure from solid angle measure to pathspace measure
    Float pdf_prv2cur_pathspace = PDFDensityConvert(prv, cur, pdf_forward_dir);
    cur->pdf_forward            = pdf_prv2cur_pathspace;

    if (++bounces == maxDepth) break;

    // Sample bsdf, continue random walk
    Vector3f wi;
    Vector3f wo = isect.wo;
    BxDFType type;
    Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf_forward_dir, BSDF_ALL, &type);

    if (pdf_forward_dir == .0f || f.IsBlack()) break;

    // update beta
    beta *= f * AbsDot(wi, isect.shading.n) / pdf_forward_dir;
    pdf_reverse_dir = isect.bsdf->Pdf(wi, wo, BSDF_ALL);
    if (type & BSDF_SPECULAR) {
      cur->delta      = true;
      pdf_forward_dir = pdf_reverse_dir = 0;
    }
    ray = isect.SpawnRay(wi);

    //* Convert the measure from solid angle measure to pathspace measure
    Float pdf_cur2prv_pathspace = PDFDensityConvert(cur, prv, pdf_reverse_dir);
    prv->pdf_inverse            = pdf_cur2prv_pathspace;
  }

  return bounces;
}

//* Generate up to maxDepth vertices
int VN_GenerateCameraSubpath(const Scene &scene, Sampler &sampler, MemoryArena &arena, int maxDepth,
                             const Camera &camera, const Point2f &p_film,
                             VN_Vertex *camera_subpath) {
  if (maxDepth == 0) return 0;

  // Generate camera ray
  CameraSample camera_sample;
  camera_sample.pFilm = p_film;
  camera_sample.time  = sampler.Get1D();
  camera_sample.pLens = sampler.Get2D();
  RayDifferential ray;
  Spectrum        beta = camera.GenerateRayDifferential(camera_sample, &ray);
  ray.ScaleDifferentials(1 / std::sqrt(sampler.samplesPerPixel));

  Float pdf_position, pdf_direction;
  camera_subpath[0] = VN_Vertex::createCamera(&camera, ray, beta);
  camera.Pdf_We(ray, &pdf_position, &pdf_direction);
  return VN_RandomWalk(scene, ray, sampler, arena, beta, pdf_direction, maxDepth - 1,
                       TransportMode::Radiance, camera_subpath + 1) +
         1;
}

int VN_GenerateLightSubpath(const Scene &scene, Sampler &sampler, MemoryArena &arena, int maxDepth,
                            Float time, const Distribution1D &light_distribution,
                            const std::unordered_map<const Light *, size_t> &light_to_index,
                            VN_Vertex                                       *light_subpath) {
  if (maxDepth == 0) return 0;
  Float           light_pdf;
  int             light_idx = light_distribution.SampleDiscrete(sampler.Get1D(), &light_pdf);
  const Light    *light     = scene.lights[light_idx].get();
  RayDifferential ray;
  Normal3f        n_light;
  Float           pdf_position, pdf_direction;
  Spectrum        Le = light->Sample_Le(sampler.Get2D(), sampler.Get2D(), time, &ray, &n_light,
                                        &pdf_position, &pdf_direction);

  if (pdf_position == 0 || pdf_direction == 0 || Le.IsBlack()) return 0;

  light_subpath[0] = VN_Vertex::createLight(light, ray, n_light, Le, pdf_position * light_pdf);
  Spectrum beta    = Le * AbsDot(n_light, ray.d) / (light_pdf * pdf_position * pdf_direction);

  return VN_RandomWalk(scene, ray, sampler, arena, beta, pdf_direction, maxDepth - 1,
                       TransportMode::Importance, light_subpath + 1) +
         1;
}

Spectrum VN_ConnectBDPT(const Scene &scene, VN_Vertex *light_subpath, VN_Vertex *camera_subpath,
                        int n_lv, int n_cv, const Distribution1D &light_distr,
                        const std::unordered_map<const Light *, size_t> &light_to_index,
                        const Camera &camera, Sampler &sampler, Point2f *p_raster, Float *misw) {
  Spectrum L(.0f);

  if (n_cv > 1 && n_lv != 0 && camera_subpath[n_cv - 1].type == VN_Vertex::Type::Light)
    return Spectrum(.0f);

  //* A new vertex may be sampled during connection
  VN_Vertex        sampled;
  const VN_Vertex *c_vtx = camera_subpath + n_cv - 1;
  const VN_Vertex *l_vtx = light_subpath + n_lv - 1;

  if (n_lv == 0 /* The path totally generated by forward sampling*/) {
    Spectrum Le   = c_vtx->Le(scene, c_vtx - 1);
    Spectrum beta = c_vtx->beta;
    L             = beta * Le;

  } else if (n_cv == 1 /* Connect a light path to the camera*/) {
    if (l_vtx->IsConnectible()) {
      VisibilityTester vis;
      Vector3f         wi;
      Float            pdf;
      Spectrum         Wi =
          camera.Sample_Wi(l_vtx->GetInteraction(), sampler.Get2D(), &wi, &pdf, p_raster, &vis);

      if (pdf > 0 && !Wi.IsBlack()) {
        sampled       = VN_Vertex::createCamera(&camera, vis.P1(), Wi / pdf);
        Spectrum beta = l_vtx->beta * sampled.beta;
        Spectrum f    = l_vtx->f(sampled, TransportMode::Importance);
        L             = beta * f;
        if (l_vtx->IsOnSurface()) L *= AbsDot(wi, l_vtx->ns());
        if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
      }
    }
  } else if (n_lv == 1 /*NEE estimation for camara subpath*/) {
    if (c_vtx->IsConnectible()) {
      Float            light_pdf;
      VisibilityTester vis;
      Vector3f         wi;
      Float            pdf;
      int              light_idx = light_distr.SampleDiscrete(sampler.Get1D(), &light_pdf);
      const Light     *light     = scene.lights[light_idx].get();
      Spectrum         light_weight =
          light->Sample_Li(c_vtx->GetInteraction(), sampler.Get2D(), &wi, &pdf, &vis);

      if (pdf > .0f && !light_weight.IsBlack()) {
        VN_EndPointInteraction ei(vis.P1(), light);
        sampled             = VN_Vertex::createLight(ei, light_weight / (pdf * light_pdf), 0);
        sampled.pdf_forward = PDFLightOrigin(scene, c_vtx, &sampled, light_distr, light_to_index);

        L = c_vtx->beta * c_vtx->f(sampled, TransportMode::Radiance) * sampled.beta;

        if (c_vtx->IsOnSurface()) L *= AbsDot(wi, c_vtx->ns());
        if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
      }
    }
  } else {
    if (c_vtx->IsConnectible() && l_vtx->IsConnectible()) {
      Spectrum beta = l_vtx->beta * c_vtx->beta;
      Spectrum f =
          l_vtx->f(*c_vtx, TransportMode::Importance) * c_vtx->f(*l_vtx, TransportMode::Radiance);
      L = beta * f;
      if (!L.IsBlack()) {
        Spectrum G = VN_G(scene, sampler, l_vtx, c_vtx);
        L *= G;
      }
    }
  }

  Float mis_weight = L.IsBlack() ? .0f
                                 : VN_BDPTMIS(scene, light_subpath, camera_subpath, sampled, n_lv,
                                              n_cv, light_distr, light_to_index);
  L *= mis_weight;
  return L;
}

Float VN_BDPTMIS(const Scene &scene, VN_Vertex *light_subpath, VN_Vertex *camera_subpath,
                 const VN_Vertex &sampled, int n_lv, int n_cv, const Distribution1D &lightPdf,
                 const std::unordered_map<const Light *, size_t> &lightToIndex) {
  if (n_lv + n_cv == 2) return 1;

  Float sum_ri = .0f;
  auto  remap0 = [](Float f) { return f != 0 ? f : 1; };

  VN_Vertex *l_vtx     = n_lv > 0 ? light_subpath + n_lv - 1 : nullptr;
  VN_Vertex *c_vtx     = n_cv > 0 ? camera_subpath + n_cv - 1 : nullptr;
  VN_Vertex *l_vtx_prv = n_lv > 1 ? l_vtx - 1 : nullptr;
  VN_Vertex *c_vtx_prv = n_cv > 1 ? c_vtx - 1 : nullptr;

  //  Temporal replace vertex if sampled in connection
  TemporalProxy<VN_Vertex> a1;
  if (n_lv == 1)
    a1 = {l_vtx, sampled};
  else if (n_cv == 1)
    a1 = {c_vtx, sampled};

  TemporalProxy<bool> a2, a3;
  if (l_vtx) a2 = {&l_vtx->delta, false};
  if (c_vtx) a3 = {&c_vtx->delta, false};

  TemporalProxy<Float> a4;
  if (c_vtx)
    a4 = {&c_vtx->pdf_inverse,
          n_lv > 0 ? PDF(l_vtx_prv, l_vtx, c_vtx, scene)
                   : PDFLightOrigin(scene, c_vtx_prv, c_vtx, lightPdf, lightToIndex)};

  TemporalProxy<Float> a5;
  if (c_vtx_prv)
    a5 = {&c_vtx_prv->pdf_inverse,
          n_lv > 0 ? PDF(l_vtx, c_vtx, c_vtx_prv, scene) : PDFLight(c_vtx, c_vtx_prv, scene)};

  TemporalProxy<Float> a6;
  if (l_vtx) a6 = {&l_vtx->pdf_inverse, PDF(c_vtx_prv, c_vtx, l_vtx, scene)};

  TemporalProxy<Float> a7;
  if (l_vtx_prv) a7 = {&l_vtx_prv->pdf_inverse, PDF(c_vtx, l_vtx, l_vtx_prv, scene)};

  Float ri = 1;
  for (int i = n_cv - 1; i > 0; --i) {
    ri *= remap0(camera_subpath[i].pdf_inverse) / remap0(camera_subpath[i].pdf_forward);
    if (!camera_subpath[i].delta && !camera_subpath[i - 1].delta) sum_ri += ri;
  }

  ri = 1;
  for (int i = n_lv - 1; i >= 0; --i) {
    ri *= remap0(light_subpath[i].pdf_inverse) / remap0(light_subpath[i].pdf_forward);
    bool delta_light_vertex = i > 0 ? light_subpath[i - 1].delta : light_subpath[0].IsDeltaLight();

    if (!light_subpath[i].delta && !delta_light_vertex) sum_ri += ri;
  }

  return 1.f / (1.f + sum_ri);

} // namespace pbrt

void VNBDPTIntegrator::Render(const Scene &scene) {
  // Construct the light distribution
  std::unique_ptr<LightDistribution> light_distribution =
      CreateLightSampleDistribution(lightSampleStrategy, scene);
  std::unordered_map<const Light *, size_t> light_to_index;
  for (size_t i = 0; i < scene.lights.size(); ++i)
    light_to_index[scene.lights[i].get()] = i;

  // Construct the volumetric light distribution
  // TODO Initialize volumetric_lights_distrib
  // TODO Initialize vlight_to_index
  const Distribution1D *vlight_distr = volumetric_lights_distrib.get();
  std::unordered_map<const VolumetricLight *, size_t> vlight_to_index;

  // Partition the image into tiles
  Film          *film          = camera->film;
  const Bounds2i sample_bounds = film->GetSampleBounds();
  const Vector2i sample_extent = sample_bounds.Diagonal();
  constexpr int  tile_size     = 16;
  const int      N_xtiles      = (sample_extent.x + tile_size - 1) / tile_size;
  const int      N_ytiles      = (sample_extent.y + tile_size - 1) / tile_size;

  ProgressReporter reporter(N_xtiles * N_ytiles, "Rendering");

  // Light Distribution is based on power
  const Distribution1D *light_distr = light_distribution->Lookup(Point3f(.0f, .0f, .0f));

  ParallelFor2D(
      [&](const Point2i tile) {
        MemoryArena arena;
        int         seed = tile.y * N_xtiles + tile.x;

        int      x0 = sample_bounds.pMin.x + tile.x * tile_size;
        int      x1 = std::min(x0 + tile_size, sample_bounds.pMax.x);
        int      y0 = sample_bounds.pMin.y + tile.y * tile_size;
        int      y1 = std::min(y0 + tile_size, sample_bounds.pMax.y);
        Bounds2i tile_bounds(Point2i(x0, y0), Point2i(x1, y1));

        std::unique_ptr<FilmTile> film_tile    = camera->film->GetFilmTile(tile_bounds);
        std::unique_ptr<Sampler>  tile_sampler = sampler->Clone(seed);

        for (Point2i p_pixel : tile_bounds) {
          tile_sampler->StartPixel(p_pixel);
          if (!InsideExclusive(p_pixel, pixelBounds)) continue;

          // Sample and estimate the pixel integral
          do {
            // Sample a location in the pixel
            Point2f p_film = (Point2f)p_pixel + tile_sampler->Get2D();

            // Alloc the vertices
            VN_Vertex *camera_subpath = arena.Alloc<VN_Vertex>(maxDepth + 2);
            VN_Vertex *light_subpath  = arena.Alloc<VN_Vertex>(maxDepth + 1);

            int N_camera_vertices = VN_GenerateCameraSubpath(
                scene, *tile_sampler, arena, maxDepth + 2, *camera, p_film, camera_subpath);
            int N_light_vertices =
                VN_GenerateLightSubpath(scene, *tile_sampler, arena, maxDepth + 1, 0, *light_distr,
                                        light_to_index, light_subpath);

            // Connect the subpaths
            Spectrum L(.0f);
            for (int n_cv = 1; n_cv <= N_camera_vertices; ++n_cv) {
              for (int n_lv = 0; n_lv <= N_light_vertices; ++n_lv) {
                // Compute the path construct by connecting camera_subpath[t]
                // and light_subpath[s]
                int depth = n_cv + n_lv - 2;
                if ((n_cv == 1 && n_lv == 1) || depth < 0 || depth > maxDepth) continue;
                // The pixel location for n_cv == 1
                Point2f p_film_new = p_film;
                Float   mis_weight = .0f;

                Spectrum L_path = VN_ConnectBDPT(scene, light_subpath, camera_subpath, n_lv, n_cv,
                                                 *light_distr, light_to_index, *camera,
                                                 *tile_sampler, &p_film_new, &mis_weight);
                if (n_cv != 1)
                  L += L_path;
                else
                  film->AddSplat(p_film_new, L_path);
              }
            }

            film_tile->AddSample(p_film, L);
            arena.Reset();
          } while (tile_sampler->StartNextSample());
        }

        film->MergeFilmTile(std::move(film_tile));
        reporter.Update();
      },
      Point2i(N_xtiles, N_ytiles));
  reporter.Done();
  film->WriteImage(1.0f / sampler->samplesPerPixel);
}

VN_Vertex VN_Vertex::createCamera(const Camera *camera, const Ray &ray, const Spectrum &beta) {
  return VN_Vertex{VN_Vertex::Type::Camera, VN_EndPointInteraction(camera, ray), beta};
}

VN_Vertex VN_Vertex::createCamera(const Camera *camera, const Interaction &it,
                                  const Spectrum &beta) {
  return VN_Vertex{VN_Vertex::Type::Camera, VN_EndPointInteraction(it, camera), beta};
}

VN_Vertex VN_Vertex::createLight(const Light *light, const Ray &ray, const Normal3f &n_light,
                                 const Spectrum &Le, Float pdf) {
  VN_Vertex v{VN_Vertex::Type::Light, VN_EndPointInteraction(light, ray, n_light), Le};
  v.pdf_forward = pdf;
  return v;
}

VN_Vertex VN_Vertex::createLight(const VN_EndPointInteraction &ei, const Spectrum &beta,
                                 Float pdf) {
  VN_Vertex v{VN_Vertex::Type::Light, ei, beta};
  v.pdf_forward = pdf;
  return v;
}

VN_Vertex VN_Vertex::createSurface(const SurfaceInteraction &si, const Spectrum &beta,
                                   const VN_Vertex &prev) {
  VN_Vertex v{si, beta};
  return v;
}

VN_Vertex VN_Vertex::createMedium(const MediumInteraction &mi, const Spectrum &beta,
                                  const VN_Vertex &prev) {
  VN_Vertex v(mi, beta);
  return v;
}

bool VN_Vertex::IsLight() const {
  return type == Type::Light || (type == Type::Surface && si.primitive->GetAreaLight());
}

bool VN_Vertex::IsDeltaLight() const {
  return type == Type::Light && ei.light && pbrt::IsDeltaLight(ei.light->flags);
}
bool VN_Vertex::IsConnectible() const {
  switch (type) {
  case Type::Medium:
  case Type::Camera:
    return true;
  case Type::Light:
    return (ei.light->flags & (int)LightFlags::DeltaDirection) == 0;
  case Type::Surface:
    return si.bsdf->NumComponents(
               BxDFType(BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)) > 0;
  }
}

bool VN_Vertex::IsOnSurface() const { return ng() != Normal3f(); }

const Interaction &VN_Vertex::GetInteraction() const {
  switch (type) {
  case Type::Medium:
    return mi;
  case Type::Surface:
    return si;
  default:
    return ei;
  }
}

Spectrum VN_Vertex::f(const VN_Vertex &next, TransportMode mode) const {
  Vector3f wi = next.p() - p();
  if (wi.LengthSquared() == 0) return 0;
  wi = Normalize(wi);
  switch (type) {
  case Type::Surface:
    return si.bsdf->f(si.wo, wi);
  case Type::Medium:
    return mi.phase->p(mi.wo, wi);
  default:
    return Spectrum(.0f);
  }
}

Spectrum VN_Vertex::Le(const Scene &scene, const VN_Vertex *towards) const {
  if (!IsLight()) return Spectrum(.0f);

  Vector3f w = towards->p() - p();
  if (w.LengthSquared() == 0) return Spectrum(.0f);
  w = Normalize(w);

  const AreaLight *light = si.primitive->GetAreaLight();
  return light->L(si, w);
}

const Normal3f &VN_Vertex::ns() const {
  return type == Type::Surface ? si.shading.n : GetInteraction().n;
}

const Normal3f &VN_Vertex::ng() const { return GetInteraction().n; }

Point3f VN_Vertex::p() const { return GetInteraction().p; }

Float PDFDensityConvert(const VN_Vertex *from, const VN_Vertex *to, Float pdf_dir) {

  Vector3f w = to->p() - from->p();
  if (w.LengthSquared() == 0) return 0;
  Float inv_distsqr = 1 / w.LengthSquared();
  w                 = Normalize(w);

  Float pdf = pdf_dir * inv_distsqr;
  pdf *= to->IsOnSurface() ? AbsDot(to->ng(), w) : 1.f;
  return pdf;
}

Float PDF(const VN_Vertex *V_minus, const VN_Vertex *V, const VN_Vertex *V_plus,
          const Scene &scene) {
  if (V->type == VN_Vertex::Type::Light) return PDFLight(V, V_plus, scene);

  Vector3f wi = V_plus->p() - V->p();
  if (wi.LengthSquared() == 0) return 0;
  wi = Normalize(wi);

  Vector3f wo;
  if (V_minus) {
    wo = V_minus->p() - V->p();
    if (wo.LengthSquared() == 0) return 0;
    wo = Normalize(wo);
  }

  Float pdf, un_used;
  if (V->type == VN_Vertex::Type::Camera)
    V->ei.camera->Pdf_We(V->ei.SpawnRay(wi), &un_used, &pdf);
  else if (V->type == VN_Vertex::Type::Surface)
    pdf = V->si.bsdf->Pdf(wo, wi);
  else if (V->type == VN_Vertex::Type::Medium)
    pdf = V->mi.phase->p(wo, wi);

  return PDFDensityConvert(V, V_plus, pdf);
}

Float PDFLightOrigin(const Scene &scene, const VN_Vertex *shading_p, const VN_Vertex *light_v,
                     const Distribution1D                            &lightDistr,
                     const std::unordered_map<const Light *, size_t> &lightToDistrIndex) {
  Vector3f w = shading_p->p() - light_v->p();
  if (w.LengthSquared() == 0) return 0;
  w = Normalize(w);

  Float        pdf_position, pdf_direction, pdf_choice;
  const Light *light = light_v->type == VN_Vertex::Type::Light
                           ? light_v->ei.light
                           : light_v->si.primitive->GetAreaLight();

  size_t idx = lightToDistrIndex.find(light)->second;
  pdf_choice = lightDistr.DiscretePDF(idx);

  light->Pdf_Le(Ray(light_v->p(), w, Infinity), light_v->ng(), &pdf_position, &pdf_direction);
  return pdf_position * pdf_choice;
}

Float PDFLight(const VN_Vertex *light_v, const VN_Vertex *v, const Scene &scene) {
  Vector3f w           = v->p() - light_v->p();
  Float    inv_distsqr = 1.f / w.LengthSquared();
  w                    = Normalize(w);
  Float pdf;

  const Light *light = light_v->type == VN_Vertex::Type::Light
                           ? light_v->ei.light
                           : light_v->si.primitive->GetAreaLight();
  Float        pdf_position, pdf_direction;
  light->Pdf_Le(Ray(light_v->p(), w, Infinity, 0), light_v->ng(), &pdf_position, &pdf_direction);
  pdf = pdf_direction * inv_distsqr;

  if (v->IsOnSurface()) pdf *= AbsDot(v->ng(), w);
  return pdf;
}

Spectrum VN_G(const Scene &scene, Sampler &sampler, const VN_Vertex *v0, const VN_Vertex *v1) {
  Vector3f d = v0->p() - v1->p();
  Float    g = 1 / d.LengthSquared();
  d *= std::sqrt(g);
  if (v0->IsOnSurface()) g *= AbsDot(v0->ns(), d);
  if (v1->IsOnSurface()) g *= AbsDot(v1->ns(), d);
  VisibilityTester vis(v0->GetInteraction(), v1->GetInteraction());
  return g * vis.Tr(scene, sampler);
}

VNBDPTIntegrator *
CreateVNBDPTIntegrator(const ParamSet &params, std::shared_ptr<Sampler> sampler,
                       std::shared_ptr<const Camera>               camera,
                       const std::vector<std::shared_ptr<Medium>> &emission_mediums) {
  int        maxDepth = params.FindOneInt("maxdepth", 5);
  int        np;
  const int *pb          = params.FindInt("pixelbounds", &np);
  Bounds2i   pixelBounds = camera->film->GetSampleBounds();
  if (pb) {
    if (np != 4)
      Error("Expected four values for \"pixelbounds\" parameter. Got %d.", np);
    else {
      pixelBounds = Intersect(pixelBounds, Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
      if (pixelBounds.Area() == 0) Error("Degenerate \"pixelbounds\" specified.");
    }
  }

  std::string lightStrategy = params.FindOneString("lightsamplestrategy", "power");

  //* Construct volumetric light sample structures
  std::vector<std::shared_ptr<VolumetricLight>> volumetric_lights;
  for (const auto emission_medium : emission_mediums) {
    std::shared_ptr<VolumetricLight> vol_light =
        std::make_shared<VolumetricLight>(emission_medium.get());
    vol_light->Initialize();
    volumetric_lights.emplace_back(vol_light);
  }

  return new VNBDPTIntegrator(sampler, camera, maxDepth, false, false, pixelBounds,
                              volumetric_lights, lightStrategy);
}

} // namespace pbrt