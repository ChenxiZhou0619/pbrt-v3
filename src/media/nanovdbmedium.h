#pragma once

#include <NanoVDB.h>
#include <util/GridHandle.h>
#include <util/IO.h>
#include <util/SampleFromVoxels.h>

#include <optional>

#include "medium.h"
#include "parallel.h"
#include "transform.h"

//* GridDensityMedium based on nanovdb
//* No scalation now

namespace pbrt {

struct MajorantGrid {
  public:
    MajorantGrid(Bounds3f world_bound,
                 Vector3i resolution = Vector3i{32, 32, 32})
        : world_bound(world_bound), resolution(resolution) {
        int capacity = resolution.x * resolution.y * resolution.z;
        majorants.resize(capacity);

        for (int i = 0; i < 3; ++i) {
            voxel_size_w[i] =
                (world_bound.pMax[i] - world_bound.pMin[i]) / resolution[i];
        }
    }

    Float at(int x, int y, int z) const;

    Float &at(int x, int y, int z);

    int size() const;

    Point3f toIndex(Point3f p_world) const;

  public:
    Bounds3f world_bound;
    Vector3i resolution;
    std::vector<Float> majorants;
    Vector3f voxel_size_w;
};

struct MajorantSeg {
    //
};

class DDATracker {
  public:
    DDATracker(bool terminate) : terminate(terminate) {}

    DDATracker(Vector3i resolution, Point3f pIndex, Vector3f direction,
               Vector3f voxel_size_w, Float cur_t_w, Float t_max_w)
        : terminate(false),
          voxel_size_w(voxel_size_w),
          cur_t_w(cur_t_w),
          t_max_w(t_max_w) {
        for (int axis = 0; axis < 3; ++axis) {
            cur_index[axis] = Clamp(pIndex[axis], 0, resolution[axis] - 1);
            delta_t_w[axis] = voxel_size_w[axis] / std::abs(direction[axis]);

            if (direction[axis] == -.0) direction[axis] = .0;

            if (direction[axis] >= .0) {
                next_crossingT_w[axis] =
                    ((cur_index[axis] + 1) - pIndex[axis]) / direction[axis];
                step[axis] = 1;
                voxel_limit[axis] = resolution[axis];
            } else {
                next_crossingT_w[axis] =
                    (cur_index[axis] - pIndex[axis]) / direction[axis];
                step[axis] = -1;
                voxel_limit[axis] = -1;
            }
        }

        for (int i = 0; i < 3; ++i) {
            next_crossingT_w[i] =
                next_crossingT_w[i] * voxel_size_w[i] + cur_t_w;
        }
    }

  public:
    bool terminate;

  private:
    Vector3f voxel_size_w;      // voxel size in world sapce
    Float cur_t_w, t_max_w;     // current t and tmax in world space
    Float next_crossingT_w[3];  // t for ray reach next boundary in world space
    Float delta_t_w[3];         // t for ray to cross a voxel in world space

    int step[3], voxel_limit[3], cur_index[3];
};

struct TrackSegment {
    int voxelIndex[3];
    Float dt;  // in world space
    Float t;   // in world space
};

//* Track a local ray in density grid space
class RegularTracker {
  public:
    RegularTracker(bool terminate) : terminate(terminate){};

    RegularTracker(const int minIndex[3], const int maxIndex[3], Point3f origin,
                   Vector3f direction, Float cur_t, Float tmax,
                   Float voxel_size);

    bool track(TrackSegment *seg);

    Float get_world_t() const;

  public:
    bool terminate;  // reach the grid bound or tmax

  private:
    Float voxel_size;         // the size of a voxel (in world space)
    Float cur_t, t_max;       // in grid space
    Float next_crossingT[3],  // t for ray to reach next boundary
        delta_t[3];           // t for ray to cross a voxel

    int step[3],         // update the index when tracking
        voxel_limit[3],  // the termination index
        cur_index[3];    // the index of current voxel
};

using BufferT = nanovdb::HostBuffer;

class NanovdbMedium : public Medium {
  public:
    NanovdbMedium(const std::string &vdbfilename, Float density_scale,
                  Spectrum sigma_a, Spectrum sigma_s, Float g,
                  const Transform &medium_transform)
        : g(g),
          density_scale(density_scale),
          sigma_a(sigma_a),
          sigma_s(sigma_s),
          sigma_maj(sigma_a + sigma_s),
          medium_transform(medium_transform) {
        // Read the nanovdb file

        densityGrid = nanovdb::io::readGrid(vdbfilename, "density", 0);
        densityFloatGrid = densityGrid.grid<Float>();

        if (!densityGrid) {
            std::cerr << ".nvdb file must contains density grid!\n";
            exit(1);
        }

        temperatureGrid = nanovdb::io::readGrid(vdbfilename, "temperature", 0);
        temperatureFloatGrid = temperatureGrid.grid<Float>();

        minIndex[0] = densityFloatGrid->indexBBox().min().x();
        minIndex[1] = densityFloatGrid->indexBBox().min().y();
        minIndex[2] = densityFloatGrid->indexBBox().min().z();

        maxIndex[0] = densityFloatGrid->indexBBox().max().x();
        maxIndex[1] = densityFloatGrid->indexBBox().max().y();
        maxIndex[2] = densityFloatGrid->indexBBox().max().z();

        auto vs = densityFloatGrid->voxelSize();

        if (vs[0] != vs[1] || vs[0] != vs[2]) {
            std::cerr << "Only support cube voxel!\n";
            exit(1);
        }
        voxel_size = vs[0];

        // Initialize maj_grid
        Point3f bound_min = indexToWorld(
                    Point3f(minIndex[0], minIndex[1], minIndex[2])),
                bound_max = indexToWorld(Point3f(
                    maxIndex[0] + 1.0, maxIndex[1] + 1.0, maxIndex[2] + 1.0));
        Bounds3f medium_worldbound{bound_min, bound_max};
        maj_grid = std::make_unique<MajorantGrid>(medium_worldbound);

        {
            int X = maj_grid->resolution.x, Y = maj_grid->resolution.y,
                Z = maj_grid->resolution.z;
            ParallelFor(
                [&](int index) {
                    int z = index % Z, y = (index / Z) % Y, x = index / (Z * Y);

                    Bounds3f wb(medium_worldbound.Lerp(Point3f(
                                    (float)x / X, (float)y / Y, (float)z / Z)),
                                medium_worldbound.Lerp(Point3f(
                                    (float)(x + 1) / X, (float)(y + 1) / Y,
                                    (float)(z + 1) / Z)));

                    auto i_min = worldToIndex(wb.pMin);
                    auto i_max = worldToIndex(wb.pMax);

                    int nx0 = std::max((int)i_min[0] - 1, minIndex[0]),
                        nx1 = std::min((int)i_max[0] + 1, maxIndex[0]),
                        ny0 = std::max((int)i_min[1] - 1, minIndex[1]),
                        ny1 = std::min((int)i_max[1] + 1, maxIndex[1]),
                        nz0 = std::max((int)i_min[2] - 1, minIndex[2]),
                        nz1 = std::min((int)i_max[2] + 1, maxIndex[2]);

                    Float max_density = .0;
                    auto accessor = densityFloatGrid->getAccessor();
                    for (int i = nx0; i <= nx1; ++i)
                        for (int j = ny0; j <= ny1; ++j)
                            for (int k = nz0; k <= nz1; ++k)
                                max_density = std::max(
                                    max_density, accessor.getValue({i, j, k}));
                    maj_grid->at(x, y, z) = max_density * density_scale;
                },
                maj_grid->size());
        }
    }

    Spectrum Tr(const Ray &ray, Sampler &sampler) const;

    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;

    bool SampleT_maj(const RayDifferential &ray, Float u, MemoryArena &arena,
                     MajorantSampleRecord *maj_record) const;

  protected:
    Float sampleDensity(int voxel_index[3]) const;

    // Transform a world point to medium local coordinate
    Point3f worldToIndex(Point3f p_world) const;

    // Transform a world vector to local coordinate
    Vector3f worldToIndex(Vector3f d_world) const;

    // Transform a index point to world space
    Point3f indexToWorld(Point3f p_index) const;

    RegularTracker get_regular_tracker(Ray ray_world) const;

    DDATracker get_dda_tracker(Ray ray_world) const;

  private:
    Float g;

    Float density_scale;
    Spectrum sigma_a, sigma_s, sigma_maj;
    Transform medium_transform;

    Float voxel_size;

    nanovdb::GridHandle<BufferT> densityGrid;
    nanovdb::GridHandle<BufferT> temperatureGrid;
    const nanovdb::FloatGrid *densityFloatGrid = nullptr;
    const nanovdb::FloatGrid *temperatureFloatGrid = nullptr;
    Float Le_scale, temperature_offset, temperature_scale;

    int minIndex[3], maxIndex[3];

    std::unique_ptr<MajorantGrid> maj_grid;
};

}  // namespace pbrt