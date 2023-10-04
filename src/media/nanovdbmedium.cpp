#include "nanovdbmedium.h"

#include "memory.h"
#include "sampler.h"
namespace pbrt {

MajorantGrid::MajorantGrid(Bounds3f world_bound, Vector3i resolution)
    : world_bound(world_bound), resolution(resolution) {
  int capacity = resolution.x * resolution.y * resolution.z;
  majorants.resize(capacity);

  for (int i = 0; i < 3; ++i) {
    voxel_size_w[i] = (world_bound.pMax[i] - world_bound.pMin[i]) / resolution[i];
  }
}

Float MajorantGrid::at(int x, int y, int z) const {
  int idx = x * resolution.y * resolution.z + y * resolution.z + z;
  return majorants[idx];
}

Float &MajorantGrid::at(int x, int y, int z) {
  int idx = x * resolution.y * resolution.z + y * resolution.z + z;
  return majorants[idx];
}

int MajorantGrid::size() const { return resolution.x * resolution.y * resolution.z; }

Point3f MajorantGrid::toIndex(Point3f p_world) const {
  Float ix_f = (p_world[0] - world_bound.pMin[0]) / voxel_size_w[0],
        iy_f = (p_world[1] - world_bound.pMin[1]) / voxel_size_w[1],
        iz_f = (p_world[2] - world_bound.pMin[2]) / voxel_size_w[2];

  return Point3f(ix_f, iy_f, iz_f);
}

DDATracker::DDATracker(Vector3i resolution, Point3f pIndex, Vector3f direction,
                       Vector3f voxel_size_w, Float cur_t_w, Float t_max_w)
    : terminate(false), voxel_size_w(voxel_size_w), cur_t_w(cur_t_w), t_max_w(t_max_w) {
  for (int axis = 0; axis < 3; ++axis) {
    cur_index[axis] = Clamp(pIndex[axis], 0, resolution[axis] - 1);
    delta_t_w[axis] = voxel_size_w[axis] / std::abs(direction[axis]);

    if (direction[axis] == -.0) direction[axis] = .0;

    if (direction[axis] >= .0) {
      next_crossingT_w[axis] = ((cur_index[axis] + 1) - pIndex[axis]) / direction[axis];
      step[axis]             = 1;
      voxel_limit[axis]      = resolution[axis];
    } else {
      next_crossingT_w[axis] = (cur_index[axis] - pIndex[axis]) / direction[axis];
      step[axis]             = -1;
      voxel_limit[axis]      = -1;
    }
  }

  for (int i = 0; i < 3; ++i) {
    next_crossingT_w[i] = next_crossingT_w[i] * voxel_size_w[i] + cur_t_w;
  }
}

bool DDATracker::track(MajorantSeg *seg) {
  if (terminate) return false;

  step_axis = -1;
  if (next_crossingT_w[0] < next_crossingT_w[1] && next_crossingT_w[0] < next_crossingT_w[2]) {
    step_axis = 0;
  } else if (next_crossingT_w[1] < next_crossingT_w[2]) {
    step_axis = 1;
  } else
    step_axis = 2;

  Float dt_w;
  if (next_crossingT_w[step_axis] > t_max_w) {
    /* Terminate in current voxel */
    //        terminate = true;
    dt_w = t_max_w - cur_t_w;

  } else {
    dt_w = next_crossingT_w[step_axis] - cur_t_w;
  }

  seg->index[0] = cur_index[0];
  seg->index[1] = cur_index[1];
  seg->index[2] = cur_index[2];
  seg->dt_w     = dt_w;

  return true;
}

void DDATracker::march(float t_w) { cur_t_w += t_w; }

void DDATracker::next() {
  if (next_crossingT_w[step_axis] > t_max_w) {
    terminate = true;
    cur_t_w   = t_max_w;
  } else {
    cur_t_w = next_crossingT_w[step_axis];
  }

  cur_index[step_axis] += step[step_axis];
  next_crossingT_w[step_axis] += delta_t_w[step_axis];

  if (cur_index[step_axis] == voxel_limit[step_axis]) terminate = true;
}

NanovdbMedium::NanovdbMedium(const std::string &vdbfilename, Float density_scale, Spectrum sigma_a,
                             Spectrum sigma_s, Float g, const Transform &medium_transform,
                             std::string density_name, std::string temperature_name, Float LeScale,
                             Float temperatureOffset, Float temperatureScale, bool sampleLe)
    : g(g), density_scale(density_scale), sigma_a(sigma_a), sigma_s(sigma_s),
      sigma_t(sigma_a + sigma_s), medium_transform(medium_transform), Le_scale(LeScale),
      temperature_offset(temperatureOffset), temperature_scale(temperatureScale) {
  // Read the nanovdb file

  densityGrid      = nanovdb::io::readGrid(vdbfilename, density_name, 1);
  densityFloatGrid = densityGrid.grid<Float>();

  if (!densityGrid) {
    std::cerr << ".nvdb file must contains density grid!\n";
    exit(1);
  }

  temperatureGrid      = nanovdb::io::readGrid(vdbfilename, temperature_name, 1);
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

  // Initialize maj_grid
  Point3f bound_min = indexToWorld(Point3f(minIndex[0], minIndex[1], minIndex[2])),
          bound_max =
              indexToWorld(Point3f(maxIndex[0] + 1.0, maxIndex[1] + 1.0, maxIndex[2] + 1.0));
  medium_worldbound = Bounds3f{bound_min, bound_max};

  printf("The bounding of medium "
         "\nmin:%.2f,%.2f,%.2f\nmax:%.2f,%.2f,%.2f\n",
         bound_min.x, bound_min.y, bound_min.z, bound_max.x, bound_max.y, bound_max.z);

  maj_grid = std::make_unique<MajorantGrid>(medium_worldbound, Vector3i{64, 64, 64});

  {
    int X = maj_grid->resolution.x, Y = maj_grid->resolution.y, Z = maj_grid->resolution.z;
    ParallelFor(
        [&](int index) {
          int z = index % Z, y = (index / Z) % Y, x = index / (Z * Y);

          Bounds3f wb(medium_worldbound.Lerp(Point3f((float)x / X, (float)y / Y, (float)z / Z)),
                      medium_worldbound.Lerp(
                          Point3f((float)(x + 1) / X, (float)(y + 1) / Y, (float)(z + 1) / Z)));

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
              for (int k = nz0; k <= nz1; ++k) {
                Float density = accessor.getValue({i, j, k});
                max_density   = std::max(max_density, density);
              }

          maj_grid->at(x, y, z) = max_density * density_scale;
        },
        maj_grid->size());
  }

  coarse_grid = std::make_unique<MajorantGrid>(medium_worldbound, Vector3i{16, 16, 16});
  {
    int X = coarse_grid->resolution.x, Y = coarse_grid->resolution.y, Z = coarse_grid->resolution.z;

    ParallelFor(
        [&](int index) {
          int z = index % Z, y = (index / Z) % Y, x = index / (Z * Y);

          Bounds3f wb(medium_worldbound.Lerp(Point3f((float)x / X, (float)y / Y, (float)z / Z)),
                      medium_worldbound.Lerp(
                          Point3f((float)(x + 1) / X, (float)(y + 1) / Y, (float)(z + 1) / Z)));

          auto i_min = worldToIndex(wb.pMin);
          auto i_max = worldToIndex(wb.pMax);

          int nx0 = std::max((int)i_min[0] - 1, minIndex[0]),
              nx1 = std::min((int)i_max[0] + 1, maxIndex[0]),
              ny0 = std::max((int)i_min[1] - 1, minIndex[1]),
              ny1 = std::min((int)i_max[1] + 1, maxIndex[1]),
              nz0 = std::max((int)i_min[2] - 1, minIndex[2]),
              nz1 = std::min((int)i_max[2] + 1, maxIndex[2]);

          Float density_sum = .0;
          Float weight      = .0;
          auto  accessor    = densityFloatGrid->getAccessor();
          for (int i = nx0; i <= nx1; ++i)
            for (int j = ny0; j <= ny1; ++j)
              for (int k = nz0; k <= nz1; ++k) {
                Float density = accessor.getValue({i, j, k});
                density_sum += density;
                weight += 1;
              }
          coarse_grid->at(x, y, z) = density_sum / weight * density_scale;
        },
        coarse_grid->size());
  }

  sample_volumetric_emission = sampleLe;

  if (sample_volumetric_emission && !temperatureGrid) {
    std::cerr << "Nanovdbmedium without temperatureGrid cann't sample Le\n";
    sample_volumetric_emission = false;
  }

  if (sample_volumetric_emission) {
    emission_grid = std::make_unique<EmissionGrid>();

    // traverse the density grid
    auto density_accessor = densityFloatGrid->getAccessor();

    for (int i = minIndex[0]; i < maxIndex[0]; ++i)
      for (int j = minIndex[1]; j < maxIndex[1]; ++j)
        for (int k = minIndex[2]; k < maxIndex[2]; ++k) {
          Float    d       = density_accessor.getValue({i, j, k});
          Point3f  p_world = indexToWorld(Point3f(i, j, k));
          Spectrum le      = Le(p_world);

          if (!le.IsBlack()) {
            Float weight = AverageRGB(le) * d;
            emission_grid->emplace_back(Vector3i(i, j, k), weight);
          }
        }

    emission_grid->build();
  }

  {
    // Compute voxel size
    Vector3f unit_v(1.0, .0, .0);
    Vector3f index_v = worldToIndex(unit_v);
    voxel_size       = 1.0 / index_v.Length();
  }
}

Spectrum NanovdbMedium::Tr(const Ray &_ray, Sampler &sampler) const {
  int channel = SampleChannel(sampler.Get1D());
  Ray ray     = _ray;
  ray.tMax *= ray.d.Length();
  ray.d = Normalize(ray.d);

  DDATracker tracker = get_dda_tracker(ray);
  Spectrum   tr(1.f);

  Float    thick_bound = -std::log(1 - sampler.Get1D());
  Spectrum sum(.0);
  Float    t_w = tracker.get_world_t();

  Spectrum weight(1.0); // mis weight

  MajorantSeg seg;
  while (tracker.track(&seg)) {
    int      ix = seg.index[0], iy = seg.index[1], iz = seg.index[2];
    Float    maj_density = maj_grid->at(ix, iy, iz);
    Float    dt          = seg.dt_w;
    Spectrum sigma_maj   = maj_density * sigma_t;

    if (sum[channel] + sigma_maj[channel] * seg.dt_w > thick_bound) {
      dt = (thick_bound - sum[channel]) / sigma_maj[channel];
      t_w += dt;
      sum += dt * sigma_maj;

      Float density = sampleDensity(ray(t_w));

      Spectrum sigma_n = (maj_density - density) * sigma_t, T_maj = Exp(-sum);

      Float pdf = T_maj[channel] * sigma_maj[channel];

      tr *= T_maj * sigma_n / pdf;
      weight *= T_maj * sigma_maj / pdf;

      sum         = Spectrum(.0f);
      thick_bound = -std::log(1 - sampler.Get1D());
      tracker.march(dt);

      if (tr.IsBlack()) return Spectrum(.0);

      if (tr.MaxComponentValue() < 0.1) {
        if (sampler.Get1D() < 0.75) {
          return Spectrum(.0);
        }
        tr /= 1 - 0.75;
      }

      continue;
    }
    sum += dt * sigma_maj;
    t_w += dt;

    tracker.next();
  }

  Spectrum T_maj = Exp(-sum);
  Float    pdf   = T_maj[channel];

  tr *= T_maj / pdf;
  weight *= T_maj / pdf;

  return tr / AverageRGB(weight);
}

Spectrum NanovdbMedium::Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                               MediumInteraction *mi) const {
  // TODO
}

Float NanovdbMedium::sampleDensity(Point3f p_world) const {
  using Sampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
  auto p_index  = worldToIndex(p_world);
  return Sampler(densityFloatGrid->tree())(p_index) * density_scale;
}

Float NanovdbMedium::sampleTemperature(Point3f p_world) const {
  if (!temperatureFloatGrid) return .0;

  auto worldToTemperatureIndex = [&](Point3f p_world) {
    p_world = Inverse(medium_transform)(p_world);
    auto p_index =
        temperatureFloatGrid->worldToIndexF(nanovdb::Vec3f(p_world[0], p_world[1], p_world[2]));
    return Point3f{p_index[0], p_index[1], p_index[2]};
  };
  using Sampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
  auto  p_index = worldToIndex(p_world);
  Float temp    = Sampler(temperatureFloatGrid->tree())(p_index);

  temp = (temp - temperature_offset) * temperature_scale;
  return temp;
}

bool NanovdbMedium::SampleT_maj(const RayDifferential &_ray, Float u_t, Float u_channel,
                                MemoryArena &arena, MajorantSampleRecord *maj_record) const {
  // \sum dt * sigma_maj = - log(1 - u)

  Ray ray = _ray;
  ray.tMax *= ray.d.Length();
  ray.d = Normalize(ray.d);

  DDATracker tracker = get_dda_tracker(ray);

  Float    thick_bound = -std::log(1 - u_t);
  Spectrum sum(.0f);
  Float    t_w = tracker.get_world_t();
  Float    maj_density;

  int channel         = SampleChannel(u_channel);
  maj_record->channel = channel;

  MajorantSeg seg;

  while (tracker.track(&seg)) {
    int ix = seg.index[0], iy = seg.index[1], iz = seg.index[2];
    maj_density        = maj_grid->at(ix, iy, iz);
    Float    dt        = seg.dt_w;
    Spectrum sigma_maj = maj_density * sigma_t;

    if (sum[channel] + sigma_maj[channel] * seg.dt_w > thick_bound) {
      dt = (thick_bound - sum[channel]) / sigma_maj[channel];

      //* Sample the information at sampled point
      t_w += dt;
      sum += dt * maj_density * sigma_t;

      maj_record->T_maj = Exp(-sum);
      maj_record->p     = ray(t_w);
      Float density     = sampleDensity(maj_record->p);
      //            maj_record->Le = Spectrum(.0f);
      maj_record->Le      = Le(maj_record->p);
      maj_record->phase   = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
      maj_record->sigma_a = density * sigma_a;
      maj_record->sigma_s = density * sigma_s;
      maj_record->sigma_n = (maj_density - density) * sigma_t;
      maj_record->t       = t_w;
      return false; // not escape from current medium
    }

    sum += dt * sigma_maj;
    t_w += dt;

    tracker.next();
  }

  maj_record->T_maj = Exp(-sum);
  return true; // escape the current medium
}

DDATracker NanovdbMedium::get_dda_tracker(Ray ray_world) const {
  // Check if ray_world intersect the maj_grid
  Point3f  pmin = maj_grid->world_bound.pMin, pmax = maj_grid->world_bound.pMax;
  Float    t_min = -Infinity, t_max = Infinity;
  Point3f  origin    = ray_world.o;
  Vector3f direction = ray_world.d;

  //    ray_world.tMax *= ray_world.d.Length();
  //    Vector3f direction = Normalize(ray_world.d);

  for (int axis = 0; axis < 3; ++axis) {
    if (direction[axis] == 0) continue;

    Float t_0 = (pmin[axis] - origin[axis]) / direction[axis],
          t_1 = (pmax[axis] - origin[axis]) / direction[axis];
    if (t_0 > t_1) std::swap(t_0, t_1);
    t_min = std::max(t_min, t_0);
    t_max = std::min(t_max, t_1);

    if (t_min > t_max || t_max < 0) return DDATracker(true); // Just terminate
  }

  Float   cur_t_w = t_min > 0 ? t_min : 0;
  Point3f pIndex  = maj_grid->toIndex(origin + direction * cur_t_w);

  return DDATracker(maj_grid->resolution, pIndex, direction, maj_grid->voxel_size_w, cur_t_w,
                    ray_world.tMax);
}

Point3f NanovdbMedium::worldToIndex(Point3f p_world) const {
  p_world = Inverse(medium_transform)(p_world);
  auto p_index =
      densityFloatGrid->worldToIndexF(nanovdb::Vec3f(p_world[0], p_world[1], p_world[2]));
  return Point3f{p_index[0], p_index[1], p_index[2]};
}

Vector3f NanovdbMedium::worldToIndex(Vector3f d_world) const {
  d_world = Inverse(medium_transform)(d_world);
  auto d_index =
      densityFloatGrid->worldToIndexDirF(nanovdb::Vec3f(d_world[0], d_world[1], d_world[2]));
  return Vector3f{d_index[0], d_index[1], d_index[2]};
}

Point3f NanovdbMedium::indexToWorld(Point3f p_index) const {
  auto p_world =
      densityFloatGrid->indexToWorldF(nanovdb::Vec3f(p_index[0], p_index[1], p_index[2]));
  return medium_transform(Point3f(p_world[0], p_world[1], p_world[2]));
}

Spectrum NanovdbMedium::Le(Point3f p_world) const {
  Float temperature = sampleTemperature(p_world);

  if (temperature <= 100.0) return Spectrum(.0);

  // Compute blackbody emission at 12 sampled wavelength
  SampledSpectrum blackbody_emission_spectrum;

  constexpr int   N_wavelengths = 12;
  constexpr Float sampled_wavelength[]{400.0,  427.27, 454.54, 481.82, 509.10, 536.36,
                                       563.64, 590.91, 618.18, 645.45, 672.72, 700.0};
  Float           Le_lambda[12];

  BlackbodyNormalized(sampled_wavelength, N_wavelengths, temperature, Le_lambda);

  blackbody_emission_spectrum =
      SampledSpectrum::FromSampled(sampled_wavelength, Le_lambda, N_wavelengths);

  auto rgb = blackbody_emission_spectrum.ToRGBSpectrum();
  for (int i = 0; i < 3; ++i)
    rgb[i] = std::max(.0f, rgb[i]);
  return rgb * Le_scale;
}

Point3f EmissionGrid::sampleVoxel(Float u, Float *pdf) const {
  int idx = voxel_distrib->SampleDiscrete(u, pdf);
  return Point3f(voxels[idx][0], voxels[idx][1], voxels[idx][2]);
}

Float EmissionGrid::P_voxel(Vector3i voxel_idx) const {
  if (voxelIndexMap.count(voxel_idx) == 0) return .0;
  auto voxel_info = voxelIndexMap.find(voxel_idx);
  return voxel_distrib->DiscretePDF(voxel_info->second);
}

Float NanovdbMedium::pdf_emissionP(Point3f p_world) const {
  Point3f p_index = worldToIndex(p_world);

  Vector3i idx_i(p_index[0], p_index[1], p_index[2]);
  Float    p = emission_grid->P_voxel(idx_i) / (voxel_size * voxel_size * voxel_size);
  return p;
}

bool NanovdbMedium::SampleEmissionPoint(Float u, Vector3f u3, VolumetricEmissionPoint *res,
                                        Float *pdf) const {
  if (!sample_volumetric_emission) return false;

  // Sample a point in the emission volume
  Float   p_voxel;
  Point3f p_index = emission_grid->sampleVoxel(u, &p_voxel) + u3;

  Point3f p_world = indexToWorld(p_index);

  res->p       = p_world;
  res->sigma_a = sampleDensity(p_world) * sigma_a;
  res->Le      = Le(p_world);
  res->pdf     = p_voxel / (voxel_size * voxel_size * voxel_size);

  if (pdf) *pdf = res->pdf;
  return true;
}

Spectrum NanovdbMedium::Tr_coarse(const Ray &_ray, Sampler &sampler) const {
  auto get_coarse_tracker = [&](Ray ray_world) {
    // Check if ray_world intersect the coarse_grid
    Point3f  pmin = coarse_grid->world_bound.pMin, pmax = coarse_grid->world_bound.pMax;
    Float    t_min = -Infinity, t_max = Infinity;
    Point3f  origin    = ray_world.o;
    Vector3f direction = ray_world.d;

    //    ray_world.tMax *= ray_world.d.Length();
    //    Vector3f direction = Normalize(ray_world.d);

    for (int axis = 0; axis < 3; ++axis) {
      if (direction[axis] == 0) continue;

      Float t_0 = (pmin[axis] - origin[axis]) / direction[axis],
            t_1 = (pmax[axis] - origin[axis]) / direction[axis];
      if (t_0 > t_1) std::swap(t_0, t_1);
      t_min = std::max(t_min, t_0);
      t_max = std::min(t_max, t_1);

      if (t_min > t_max || t_max < 0) return DDATracker(true); // Just terminate
    }

    Float   cur_t_w = t_min > 0 ? t_min : 0;
    Point3f pIndex  = coarse_grid->toIndex(origin + direction * cur_t_w);

    return DDATracker(coarse_grid->resolution, pIndex, direction, coarse_grid->voxel_size_w,
                      cur_t_w, ray_world.tMax);
  };

  Ray ray = _ray;
  ray.tMax *= ray.d.Length();
  ray.d = Normalize(ray.d);

  DDATracker tracker = get_coarse_tracker(ray);

  Spectrum sum(.0);

  MajorantSeg seg;
  while (tracker.track(&seg)) {
    int   ix = seg.index[0], iy = seg.index[1], iz = seg.index[2];
    Float density = coarse_grid->at(ix, iy, iz);

    sum += density * sigma_t * seg.dt_w;
    tracker.next();
  }

  return Exp(-sum);
}

} // namespace pbrt