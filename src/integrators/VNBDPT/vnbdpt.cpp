#include "vnbdpt.h"
#include "film.h"
#include "filters/box.h"
#include "integrator.h"
#include "lightdistrib.h"
#include "paramset.h"
#include "progressreporter.h"
#include "sampler.h"
#include "stats.h"

namespace pbrt
{

VNBDPTIntegrator::VNBDPTIntegrator(std::shared_ptr<Sampler>      sampler,
                                   std::shared_ptr<const Camera> camera, int maxDepth,
                                   bool visualizeStrategies, bool visualizeWeights,
                                   const Bounds2i&    pixelBounds,
                                   const std::string& lightSampleStrategy)
    : sampler(sampler), camera(camera), maxDepth(maxDepth),
      visualizeStrategies(visualizeStrategies), visualizeWeights(visualizeWeights),
      pixelBounds(pixelBounds), lightSampleStrategy(lightSampleStrategy)
{
  // TODO
}

int VN_RandomWalk(const Scene& scene, RayDifferential ray, Sampler& sampler, MemoryArena& arena,
                  Spectrum beta, Float pdf, int maxDepth, TransportMode mode, VN_Vertex* subpath)
{
  if (maxDepth == 0)
    return 0;

  int   bounces               = 0;
  Float pdf_forward_tocur_dir = pdf;
  Float pdf_inverse_toprv_dir = 0;

  while (true)
  {
    SurfaceInteraction isect;
    bool               found_intersection = scene.Intersect(ray, &isect);

    if (beta.IsBlack())
      break;

    VN_Vertex& cur_vertex = subpath[bounces];
    VN_Vertex& prv_vertex = subpath[bounces - 1];

    if (!found_intersection)
    {
      if (mode == TransportMode::Radiance /* The subpath is construct from camera*/)
      {
        // TODO Create an infinite light vertex
      }
      break;
    }

    isect.ComputeScatteringFunctions(ray, arena, true, mode);
    if (!isect.bsdf)
    {
      ray = isect.SpawnRay(ray.d);
      continue;
    }

    // Create a vertex at current shading point
    cur_vertex = VN_Vertex::createSurface(isect, beta, pdf_forward_tocur_dir, prv_vertex);

    if (++bounces == maxDepth /* Reach the upper limit*/)
      break;

    // Sampling bsdf at current vertex and compute the inverse pdf of sampling prev vertex
    Vector3f wo = isect.wo;
    Vector3f wi;
    BxDFType type;
    Spectrum f =
        isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf_forward_tocur_dir, BSDF_ALL, &type);

    if (f.IsBlack() || pdf_forward_tocur_dir == .0f)
      break;
    beta *= f * AbsDot(wi, isect.shading.n) / pdf_forward_tocur_dir;
    pdf_inverse_toprv_dir = isect.bsdf->Pdf(wi, wo, BSDF_ALL);
    if (type & BSDF_SPECULAR)
    {
      cur_vertex.delta      = true;
      pdf_forward_tocur_dir = pdf_inverse_toprv_dir = .0f;
    }
    ray = isect.SpawnRay(wi);

    prv_vertex.pdf_inverse_dir = pdf_inverse_toprv_dir;
  }

  return bounces;
}

//* Generate up to maxDepth vertices
int VN_GenerateCameraSubpath(const Scene& scene, Sampler& sampler, MemoryArena& arena, int maxDepth,
                             const Camera& camera, const Point2f& p_film, VN_Vertex* camera_subpath)
{
  if (maxDepth == 0)
    return 0;

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

int VN_GenerateLightSubpath(const Scene& scene, Sampler& sampler, MemoryArena& arena, int maxDepth,
                            Float time, const Distribution1D& light_distribution,
                            const std::unordered_map<const Light*, size_t>& light_to_index,
                            VN_Vertex*                                      light_subpath)
{
  if (maxDepth == 0)
    return 0;
  Float           light_pdf;
  int             light_idx = light_distribution.SampleDiscrete(sampler.Get1D(), &light_pdf);
  const Light*    light     = scene.lights[light_idx].get();
  RayDifferential ray;
  Normal3f        n_light;
  Float           pdf_position, pdf_direction;
  Spectrum        Le = light->Sample_Le(sampler.Get2D(), sampler.Get2D(), time, &ray, &n_light,
                                        &pdf_position, &pdf_direction);

  if (pdf_position == 0 || pdf_direction == 0 || Le.IsBlack())
    return 0;

  light_subpath[0] = VN_Vertex::createLight(light, ray, n_light, Le, pdf_position * light_pdf);
  Spectrum beta    = Le * AbsDot(n_light, ray.d) / (light_pdf * pdf_position * pdf_direction);

  // TODO Handle infinite light
  return VN_RandomWalk(scene, ray, sampler, arena, beta, pdf_direction, maxDepth - 1,
                       TransportMode::Importance, light_subpath + 1) +
         1;
}

Spectrum VN_ConnectBDPT(const Scene& scene, VN_Vertex* light_subpath, VN_Vertex* camera_subpath,
                        int n_lv, int n_cv, const Distribution1D& light_distr,
                        const std::unordered_map<const Light*, size_t>& light_to_index,
                        const Camera& camera, Sampler& sampler, Point2f* p_raster)
{
  Spectrum L(.0f);

  if (n_cv > 1 && n_lv != 0 && camera_subpath[n_cv - 1].type == VN_Vertex::Type::Light)
    return L;

  VN_Vertex sampled;
  if (n_lv == 0 /* The path totally generated by forward sampling*/)
  {
    const VN_Vertex& path_end = camera_subpath[n_cv - 1];
    if (path_end.IsLight())
      L = path_end.Le(scene, camera_subpath[n_cv - 2]) * path_end.beta;
  }
  else if (n_cv == 1 /* Connect a light path to the camera*/)
  {
    const VN_Vertex& light_vertex = light_subpath[n_lv - 1];
    if (light_vertex.IsConnectible())
    {
      VisibilityTester vis;
      Vector3f         wi;
      Float            pdf;
      Spectrum Wi = camera.Sample_Wi(light_vertex.GetInteraction(), sampler.Get2D(), &wi, &pdf,
                                     p_raster, &vis);

      if (pdf > 0 && !Wi.IsBlack())
      {
        // TODO
        sampled = VN_Vertex::createCamera(&camera, vis.P1(), Wi / pdf);
        L = light_vertex.beta * light_vertex.f(sampled, TransportMode::Importance) * sampled.beta;

        if (light_vertex.IsOnSurface())
          L *= AbsDot(wi, light_vertex.ns());

        if (!L.IsBlack())
          L *= vis.Tr(scene, sampler);
      }
    }
  }
  else if (n_lv == 1 /*Nee estimation for camara subpath*/)
  {
    const VN_Vertex& camera_vertex = camera_subpath[n_cv - 1];
    if (camera_vertex.IsConnectible())
    {
      Float            light_pdf;
      VisibilityTester vis;
      Vector3f         wi;
      Float            pdf;
      int              light_idx = light_distr.SampleDiscrete(sampler.Get1D(), &light_pdf);
      const Light*     light     = scene.lights[light_idx].get();
      Spectrum         light_weight =
          light->Sample_Li(camera_subpath->GetInteraction(), sampler.Get2D(), &wi, &pdf, &vis);

      if (pdf > .0f && !light_weight.IsBlack())
      {
        // TODO create light_vertex
        sampled;

        L = camera_vertex.beta * camera_vertex.f(sampled, TransportMode::Radiance) * sampled.beta;

        if (camera_vertex.IsOnSurface())
          L *= AbsDot(wi, camera_vertex.ns());

        if (!L.IsBlack())
          L *= vis.Tr(scene, sampler);
      }
    }
  }
  else
  {
    const VN_Vertex& camera_vertex = camera_subpath[n_cv - 1];
    const VN_Vertex& light_vertex  = light_subpath[n_lv - 1];
    if (camera_vertex.IsConnectible() && light_vertex.IsConnectible())
    {
      L = camera_vertex.beta * camera_vertex.f(light_vertex, TransportMode::Radiance) *
          light_vertex.f(camera_vertex, TransportMode::Importance) * light_vertex.beta;
      if (!L.IsBlack())
      {
        // TODO G is geometry term between cv and lv
        Float G;
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

void VNBDPTIntegrator::Render(const Scene& scene)
{
  // Construct the light distribution
  std::unique_ptr<LightDistribution> light_distribution =
      CreateLightSampleDistribution(lightSampleStrategy, scene);
  std::unordered_map<const Light*, size_t> light_to_index;
  for (size_t i = 0; i < scene.lights.size(); ++i)
    light_to_index[scene.lights[i].get()] = i;

  // Partition the image into tiles
  Film*          film          = camera->film;
  const Bounds2i sample_bounds = film->GetSampleBounds();
  const Vector2i sample_extent = sample_bounds.Diagonal();
  constexpr int  tile_size     = 16;
  const int      N_xtiles      = (sample_extent.x + tile_size - 1) / tile_size;
  const int      N_ytiles      = (sample_extent.y + tile_size - 1) / tile_size;

  ProgressReporter reporter(N_xtiles * N_ytiles, "Rendering");

  // Light Distribution is based on power
  const Distribution1D* light_distr = light_distribution->Lookup(Point3f(.0f, .0f, .0f));

  ParallelFor2D(
      [&](const Point2i tile)
      {
        MemoryArena arena;
        int         seed = tile.y * N_xtiles + tile.x;

        int      x0 = sample_bounds.pMin.x + tile.x * tile_size;
        int      x1 = std::min(x0 + tile_size, sample_bounds.pMax.x);
        int      y0 = sample_bounds.pMin.y + tile.y * tile_size;
        int      y1 = std::min(y0 + tile_size, sample_bounds.pMax.y);
        Bounds2i tile_bounds(Point2i(x0, y0), Point2i(x1, y1));

        std::unique_ptr<FilmTile> film_tile    = camera->film->GetFilmTile(tile_bounds);
        std::unique_ptr<Sampler>  tile_sampler = sampler->Clone(seed);

        for (Point2i p_pixel : tile_bounds)
        {
          tile_sampler->StartPixel(p_pixel);
          if (!InsideExclusive(p_pixel, pixelBounds))
            continue;

          // Sample and estimate the pixel integral
          do
          {
            // Sample a location in the pixel
            Point2f p_film = (Point2f)p_pixel + tile_sampler->Get2D();

            // Alloc the vertices
            VN_Vertex* camera_subpath = arena.Alloc<VN_Vertex>(maxDepth + 2);
            VN_Vertex* light_subpath  = arena.Alloc<VN_Vertex>(maxDepth + 1);

            int N_camera_vertices = VN_GenerateCameraSubpath(
                scene, *tile_sampler, arena, maxDepth + 2, *camera, p_film, camera_subpath);
            int N_light_vertices =
                VN_GenerateLightSubpath(scene, *tile_sampler, arena, maxDepth + 1, 0, *light_distr,
                                        light_to_index, light_subpath);

            // Connect the subpaths
            Spectrum L(.0f);
            for (int n_cv = 1; n_cv < N_camera_vertices; ++n_cv)
            {
              for (int n_lv = 0; n_lv < N_light_vertices; ++n_lv)
              {
                // Compute the path construct by connecting camera_subpath[t] and light_subpath[s]
                int depth = n_cv + n_lv - 2;
                if ((n_cv == 1 && n_lv == 1) || depth < 0 || depth > maxDepth)
                  continue;
                //* The pixel location for n_cv == 1
                Point2f p_film_new = p_film;
                Float   mis_weight = .0f;

                Spectrum L_path =
                    VN_ConnectBDPT(scene, light_subpath, camera_subpath, n_lv, n_cv, *light_distr,
                                   light_to_index, *camera, *sampler, &p_film_new, &mis_weight);
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

} // namespace pbrt