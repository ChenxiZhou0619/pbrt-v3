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