#pragma once

#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/GridHandle.h>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>

#include <unordered_map>

#include "medium.h"
#include "parallel.h"
#include "sampling.h"
#include "transform.h"

//* GridDensityMedium based on nanovdb
//* No scalation now

namespace pbrt {

struct MajorantGrid {
public:
  MajorantGrid(Bounds3f world_bound, Vector3i resolution = Vector3i{32, 32, 32});

  Float at(int x, int y, int z) const;

  Float &at(int x, int y, int z);

  int size() const;

  Point3f toIndex(Point3f p_world) const;

public:
  Bounds3f           world_bound;
  Vector3i           resolution;
  std::vector<Float> majorants;
  Vector3f           voxel_size_w;
};

struct MajorantSeg {
  Float    dt_w;
  Vector3i index;
};

class DDATracker {
public:
  DDATracker(bool terminate) : terminate(terminate) {}

  DDATracker(Vector3i resolution, Point3f pIndex, Vector3f direction, Vector3f voxel_size_w,
             Float cur_t_w, Float t_max_w);

  bool track(MajorantSeg *seg);

  // still in same voxel
  void march(float t_w);

  // to next voxel
  void next();

  Float get_world_t() const { return cur_t_w; }

public:
  bool terminate;

private:
  Vector3f voxel_size_w;        // voxel size in world sapce
  Float    cur_t_w, t_max_w;    // current t and tmax in world space
  Float    next_crossingT_w[3]; // t for ray reach next boundary in world space
  Float    delta_t_w[3];        // t for ray to cross a voxel in world space

  int step[3], voxel_limit[3], cur_index[3];

  int step_axis;
};

class EmissionGrid {
public:
  EmissionGrid() = default;

  void emplace_back(Vector3i index, Float weight) {
    voxelIndexMap[index] = voxels.size();
    voxels.emplace_back(index);
    weights.emplace_back(weight);
  }

  void build() {
    if (!voxels.empty()) {
      size          = voxels.size();
      voxel_distrib = std::make_unique<Distribution1D>(weights.data(), size);
    }
  }

  Point3f sampleVoxel(Float u, Float *pdf = nullptr) const;

  Float P_voxel(Vector3i voxel_idx) const;

private:
  int                               size;
  std::vector<Vector3i>             voxels;
  std::vector<Float>                weights;
  std::unique_ptr<Distribution1D>   voxel_distrib;
  std::unordered_map<Vector3i, int> voxelIndexMap;
};

using BufferT = nanovdb::HostBuffer;
class LightGridHierarchy;

class NanovdbMedium : public Medium {
public:
  NanovdbMedium(const std::string &vdbfilename, Float density_scale, Spectrum sigma_a,
                Spectrum sigma_s, Float g, const Transform &medium_transform,
                std::string density_name, std::string temperature_name, Float LeScale,
                Float temperatureOffset, Float temperatureScale, bool sampleLe);

  Spectrum Tr(const Ray &ray, Sampler &sampler) const;

  virtual Spectrum Tr_coarse(const Ray &ray, Sampler &sampler) const;

  Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                  MediumInteraction *mi) const;

  bool SampleT_maj(const RayDifferential &ray, Float u_t, Float u_channel, MemoryArena &arena,
                   MajorantSampleRecord *maj_record) const;

  virtual bool SampleEmissionPoint(Float u, Vector3f u3, VolumetricEmissionPoint *res,
                                   Float *pdf = nullptr) const;

  virtual Float pdf_emissionP(Point3f p_world) const;

protected:
  Float sampleDensity(Point3f p_world) const;

  Float sampleTemperature(Point3f p_world) const;

  // Transform a world point to medium local coordinate
  Point3f worldToIndex(Point3f p_world) const;

  // Transform a world vector to local coordinate
  Vector3f worldToIndex(Vector3f d_world) const;

  // Transform a index point to world space
  Point3f indexToWorld(Point3f p_index) const;

  DDATracker get_dda_tracker(Ray ray_world) const;

  Spectrum Le(Point3f p_world) const;

  //* Friend declarations
  friend std::unique_ptr<LightGridHierarchy> CreateLGH(std::shared_ptr<NanovdbMedium> media,
                                                       int N_hierarchies);

  friend void InitializeDeepShadowmap(const NanovdbMedium &media, LightGridHierarchy &lgh);

private:
  Float g;

  Float     density_scale;
  Spectrum  sigma_a, sigma_s, sigma_t;
  Transform medium_transform;

  nanovdb::GridHandle<BufferT> densityGrid;
  nanovdb::GridHandle<BufferT> temperatureGrid;
  const nanovdb::FloatGrid    *densityFloatGrid     = nullptr;
  const nanovdb::FloatGrid    *temperatureFloatGrid = nullptr;
  Float                        Le_scale, temperature_offset, temperature_scale;

  int      minIndex[3], maxIndex[3];
  Bounds3f medium_worldbound;

  std::unique_ptr<MajorantGrid> maj_grid;
  std::unique_ptr<EmissionGrid> emission_grid;
  std::unique_ptr<MajorantGrid> coarse_grid;

  Float voxel_size;
};

} // namespace pbrt