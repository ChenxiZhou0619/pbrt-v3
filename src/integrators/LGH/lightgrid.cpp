#include "lightgrid.h"
#include <media/nanovdbmedium.h>
#include <parallel.h>
namespace pbrt {

//*--------------------------------------------------------------------------------
LightGrid::LightGrid(Float voxel_size, Vector3i resolution, Point3f p_min)
    : voxel_size(voxel_size), resolution(resolution) {
  int N_vertex = (resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1);
  vertices     = std::vector<GridVertex>(N_vertex);

  for (int x = 0; x <= resolution[0]; ++x) {
    for (int y = 0; y <= resolution[1]; ++y) {
      for (int z = 0; z <= resolution[2]; ++z) {
        Point3f vertex_position           = p_min + Vector3f(x, y, z) * voxel_size;
        this->at(x, y, z).vertex_position = vertex_position;
      }
    }
  }
}

GridVertex &LightGrid::at(int x, int y, int z) {
  int offset = x * (resolution[1] + 1) * (resolution[2] + 1) + y * (resolution[2] + 1) + z;
  return vertices[offset];
}

const GridVertex &LightGrid::at(int x, int y, int z) const {
  int offset = x * (resolution[1] + 1) * (resolution[2] + 1) + y * (resolution[2] + 1) + z;
  return vertices[offset];
}

//*--------------------------------------------------------------------------------
LightGridHierarchy::LightGridHierarchy(int N_hierarchies) : N_hierarchies(N_hierarchies) {
  grids = std::vector<std::unique_ptr<LightGrid>>(N_hierarchies);
}

std::unique_ptr<LightGridHierarchy> CreateLGH(std::shared_ptr<NanovdbMedium> media,
                                              int                            N_hierarchies) {

  Bounds3f medium_worldbound = media->medium_worldbound;

  Point3f  center      = medium_worldbound.Lerp(Point3f(.5f, .5f, .5f));
  Vector3f diagnol     = medium_worldbound.Diagonal();
  Float    diagnol_max = std::max(diagnol.x, std::max(diagnol.y, diagnol.z));

  Point3f p_min = center - Vector3f(1.f, 1.f, 1.f) * diagnol_max * .5f;
  Point3f p_max = center + Vector3f(1.f, 1.f, 1.f) * diagnol_max * .5f;

  std::unique_ptr<LightGridHierarchy> lgh = std::make_unique<LightGridHierarchy>(N_hierarchies);

  //* The resolution of level i is 2^i
  Float    cell_size = diagnol_max;
  Vector3i resolution(1, 1, 1);
  for (int level = 0; level < N_hierarchies; ++level) {
    lgh->grids[level] = std::make_unique<LightGrid>(cell_size, resolution, p_min);

    cell_size *= 0.5f;
    resolution *= 2;
  }

  //* Building the hierarchy, distributing the energy into grid verticies

  int voxel_res_x = media->maxIndex[0] - media->minIndex[0] + 1;
  int voxel_res_y = media->maxIndex[1] - media->minIndex[1] + 1;
  int voxel_res_z = media->maxIndex[2] - media->minIndex[2] + 1;

  int N_voxel = voxel_res_x * voxel_res_y * voxel_res_z;

  ParallelFor(
      [&](int emission_voxel_index) {
        int voxel_index_x = emission_voxel_index / (voxel_res_y * voxel_res_z) + media->minIndex[0];
        int voxel_index_y = (emission_voxel_index / voxel_res_z) % voxel_res_y + media->minIndex[1];
        int voxel_index_z = emission_voxel_index % voxel_res_z + media->minIndex[2];

        Point3f voxel_center_world = media->indexToWorld(
            Point3f(voxel_index_x + .5f, voxel_index_y + .5f, voxel_index_z + .5f));

        Spectrum sigma_a = media->sampleDensity(voxel_center_world) * media->sigma_a;
        Spectrum Le      = media->Le(voxel_center_world);

        if ((sigma_a * Le).IsBlack()) return;

        Float voxel_volume = std::pow(media->voxel_size, 3);

        //* A voxel will updates 8 verticies in each hierarchy
        Vector3f diag = (voxel_center_world - p_min);

        for (int level = 0; level < N_hierarchies; ++level) {
          Float    cell_size  = lgh->grids[level]->voxel_size;
          Vector3f cell_index = diag / cell_size;

          int ix = cell_index.x;
          int iy = cell_index.y;
          int iz = cell_index.z;

          for (int order = 0; order < 8; ++order) {
            int x = ix + ((order & 0b0100) ? 1 : 0); //! operator priority
            int y = iy + ((order & 0b0010) ? 1 : 0);
            int z = iz + ((order & 0b0001) ? 1 : 0);

            GridVertex &vtx = lgh->grids[level]->at(x, y, z);

            Float trilinear_weight = 1.f;
            trilinear_weight *= 1.f - std::abs(cell_index.x - x);
            trilinear_weight *= 1.f - std::abs(cell_index.y - y);
            trilinear_weight *= 1.f - std::abs(cell_index.z - z);

            //* Update the intensity and illumination_center of the vertex
            vtx.mtx.lock();
            Spectrum weighted_intensity  = trilinear_weight * sigma_a * Le * voxel_volume;
            Float    weight              = weighted_intensity.y();
            Vector3f illumination_offset = weight * (voxel_center_world - vtx.vertex_position);

            vtx.intensity += weighted_intensity;
            vtx.illumination_offset_weight += weight;
            vtx.illumination_offset += weight * illumination_offset;

            vtx.mtx.unlock();
          }
        }
      },
      N_voxel);

  for (int level = 0; level < N_hierarchies; ++level) {
    int N_vertex = lgh->grids[level]->vertices.size();
    int X        = lgh->grids[level]->resolution[0] + 1;
    int Y        = lgh->grids[level]->resolution[1] + 1;
    int Z        = lgh->grids[level]->resolution[2] + 1;

    ParallelFor(
        [&](int vertex_index) {
          int x = vertex_index / (Y * Z);
          int y = (vertex_index / Z) % Y;
          int z = vertex_index % Z;

          GridVertex &vtx         = lgh->grids[level]->at(x, y, z);
          vtx.illumination_center = vtx.vertex_position;
          if (vtx.illumination_offset_weight != .0f)
            vtx.illumination_center += vtx.illumination_offset / vtx.illumination_offset_weight;

          //* Construct deep shadowmap for every vertex with non-zero intensity
          if (!vtx.intensity.IsBlack()) {
            // TODO compute the resolution and r_l
            int   resolution;
            Float r_l;

            //* Just malloc the memory
            vtx.deepshadowmap =
                std::make_unique<DeepShadowMap>(resolution, r_l, vtx.illumination_center);
          }
        },
        N_vertex);
  }

  return std::move(lgh);
}

void InitializeDeepShadowmap_(const NanovdbMedium &media, LightGridHierarchy &lgh) {
  int N = lgh.N_hierarchies;
  for (int level = 0; level < N; ++level) {
    int N_vertex = lgh.grids[level]->vertices.size();
    int X        = lgh.grids[level]->resolution[0] + 1;
    int Y        = lgh.grids[level]->resolution[1] + 1;
    int Z        = lgh.grids[level]->resolution[2] + 1;

    // TODO
    auto queryFilterdDensity = [&](Point3f p_world, Float filter_size) -> Float { return .0f; };
    auto accmulateDensity    = [&](Point3f from, Vector3f direction, Float tMax,
                                Float accmulated_densities[4]) {};

    ParallelFor(
        [&](int vertex_index) {
          int x = vertex_index / (Y * Z);
          int y = (vertex_index / Z) % Y;
          int z = vertex_index % Z;

          GridVertex &vtx = lgh.grids[level]->at(x, y, z);

          if (!vtx.deepshadowmap) return;

          //* Compute accumulated filtered density towards each texel
          Point3f from       = vtx.illumination_center;
          int     resolution = vtx.deepshadowmap->resolution;
          // TODO compute tmax
          Float tMax;
          for (int face = 0; face < 6; ++face) {
            for (int u = 0; u < resolution; ++u) {
              for (int v = 0; v < resolution; ++v) {
                //* Compute the direction towards each texel
                Vector3f direction = vtx.deepshadowmap->texelToDirection(face, u, v);
                Float    accmulated_densities[4];
                accmulateDensity(from, direction, tMax, accmulated_densities);
                //                vtx.deepshadowmap->setCubeMap(accmulated_densities, u, v);
              }
            }
          }
        },
        N);
  }
}
}; // namespace pbrt