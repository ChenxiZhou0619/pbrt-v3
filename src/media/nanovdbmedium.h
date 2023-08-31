#pragma once

#include <NanoVDB.h>
#include <util/GridHandle.h>
#include <util/IO.h>
#include <util/SampleFromVoxels.h>

#include <optional>

#include "medium.h"
#include "transform.h"

//* GridDensityMedium based on nanovdb
//* No scalation now

namespace pbrt {

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

    RegularTracker get_regular_tracker(Ray ray_world) const;

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
};

}  // namespace pbrt