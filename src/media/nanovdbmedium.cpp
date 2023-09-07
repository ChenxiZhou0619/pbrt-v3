#include "nanovdbmedium.h"

#include "memory.h"
#include "sampler.h"
namespace pbrt {

Float MajorantGrid::at(int x, int y, int z) const {
    int idx = x * resolution.y * resolution.z + y * resolution.z + z;
    return majorants[idx];
}

Float &MajorantGrid::at(int x, int y, int z) {
    int idx = x * resolution.y * resolution.z + y * resolution.z + z;
    return majorants[idx];
}

int MajorantGrid::size() const {
    return resolution.x * resolution.y * resolution.z;
}

Point3f MajorantGrid::toIndex(Point3f p_world) const {
    Float ix_f = (p_world[0] - world_bound.pMin[0]) / voxel_size_w[0],
          iy_f = (p_world[1] - world_bound.pMin[1]) / voxel_size_w[1],
          iz_f = (p_world[2] - world_bound.pMin[2]) / voxel_size_w[2];

    return Point3f(ix_f, iy_f, iz_f);
}

bool DDATracker::track(MajorantSeg *seg) {
    if (terminate) return false;

    step_axis = -1;
    if (next_crossingT_w[0] < next_crossingT_w[1] &&
        next_crossingT_w[0] < next_crossingT_w[2]) {
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
    seg->dt_w = dt_w;

    return true;
}

void DDATracker::march(float t_w) { cur_t_w += t_w; }

void DDATracker::next() {
    if (next_crossingT_w[step_axis] > t_max_w) {
        terminate = true;
        cur_t_w = t_max_w;
    } else {
        cur_t_w = next_crossingT_w[step_axis];
    }

    cur_index[step_axis] += step[step_axis];
    next_crossingT_w[step_axis] += delta_t_w[step_axis];

    if (cur_index[step_axis] == voxel_limit[step_axis]) terminate = true;
}

RegularTracker::RegularTracker(const int minIndex[3], const int maxIndex[3],
                               Point3f origin, Vector3f direction, Float cur_t,
                               Float tmax, Float voxel_size)
    : terminate(false), voxel_size(voxel_size), cur_t(cur_t), t_max(tmax) {
    for (int axis = 0; axis < 3; ++axis) {
        cur_index[axis] = Clamp(origin[axis], minIndex[axis], maxIndex[axis]);
        delta_t[axis] = 1.0 / std::abs(direction[axis]);

        if (direction[axis] == -.0) direction[axis] = .0;
        if (direction[axis] >= .0) {
            next_crossingT[axis] =
                (cur_index[axis] + 1 - origin[axis]) / direction[axis];
            step[axis] = 1;
            voxel_limit[axis] = maxIndex[axis] + 1;
        } else {
            next_crossingT[axis] =
                (cur_index[axis] - origin[axis]) / direction[axis];
            step[axis] = -1;
            voxel_limit[axis] = minIndex[axis] - 1;
        }
    }
    for (int i = 0; i < 3; ++i) next_crossingT[i] += cur_t;
}

Float RegularTracker::get_world_t() const { return cur_t * voxel_size; }

bool RegularTracker::track(TrackSegment *seg) {
    if (terminate) return false;

    int step_axis = -1;
    if (next_crossingT[0] < next_crossingT[1] &&
        next_crossingT[0] < next_crossingT[2]) {
        // Step along x axis
        step_axis = 0;
    } else if (next_crossingT[1] < next_crossingT[2]) {
        // step along y axis
        step_axis = 1;
    } else
        step_axis = 2;

    Float dt_grid;
    if (next_crossingT[step_axis] > t_max /* Terminate in current voxel */) {
        terminate = true;
        dt_grid = t_max - cur_t;
        cur_t = t_max;
    } else { /* step into next voxel */
        dt_grid = next_crossingT[step_axis] - cur_t;
        cur_t = next_crossingT[step_axis];
    }

    seg->voxelIndex[0] = cur_index[0];
    seg->voxelIndex[1] = cur_index[1];
    seg->voxelIndex[2] = cur_index[2];
    seg->dt = dt_grid * voxel_size;
    seg->t = cur_t * voxel_size;

    // Update the current voxel
    cur_index[step_axis] += step[step_axis];
    next_crossingT[step_axis] += delta_t[step_axis];

    if (cur_index[step_axis] == voxel_limit[step_axis]) terminate = true;
    return true;
}

Spectrum NanovdbMedium::Tr(const Ray &_ray, Sampler &sampler) const {
    int channel = SampleChannel(sampler.Get1D());
    Ray ray = _ray;
    ray.tMax *= ray.d.Length();
    ray.d = Normalize(ray.d);

    DDATracker tracker = get_dda_tracker(ray);
    Spectrum tr(1.f);

    Float thick_bound = -std::log(1 - sampler.Get1D());
    Spectrum sum(.0);
    Float t_w = tracker.get_world_t();

    Spectrum weight(1.0);  // mis weight

    MajorantSeg seg;
    while (tracker.track(&seg)) {
        int ix = seg.index[0], iy = seg.index[1], iz = seg.index[2];
        Float maj_density = maj_grid->at(ix, iy, iz);
        Float dt = seg.dt_w;
        Spectrum sigma_maj = maj_density * sigma_t;

        if (sum[channel] + sigma_maj[channel] * seg.dt_w > thick_bound) {
            dt = (thick_bound - sum[channel]) / sigma_maj[channel];
            t_w += dt;
            sum += dt * sigma_maj;

            Float density = sampleDensity(ray(t_w));

            Spectrum sigma_n = (maj_density - density) * sigma_t,
                     T_maj = Exp(-sum);

            Float pdf = T_maj[channel] * sigma_maj[channel];

            tr *= T_maj * sigma_n / pdf;
            weight *= T_maj * sigma_maj / pdf;

            sum = Spectrum(.0f);
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
    Float pdf = T_maj[channel];

    tr *= T_maj / pdf;
    weight *= T_maj / pdf;

    return tr / AverageRGB(weight);
}

Spectrum NanovdbMedium::Sample(const Ray &ray, Sampler &sampler,
                               MemoryArena &arena,
                               MediumInteraction *mi) const {
    // TODO
}

Float NanovdbMedium::sampleDensity(Point3f p_world) const {
    using Sampler =
        nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
    auto p_index = worldToIndex(p_world);
    return Sampler(densityFloatGrid->tree())(p_index) * density_scale;
}

Float NanovdbMedium::sampleTemperature(Point3f p_world) const {
    if (!temperatureFloatGrid) return .0;

    auto worldToTemperatureIndex = [&](Point3f p_world) {
        p_world = Inverse(medium_transform)(p_world);
        auto p_index = temperatureFloatGrid->worldToIndexF(
            nanovdb::Vec3f(p_world[0], p_world[1], p_world[2]));
        return Point3f{p_index[0], p_index[1], p_index[2]};
    };
    using Sampler =
        nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
    auto p_index = worldToIndex(p_world);
    Float temp = Sampler(temperatureFloatGrid->tree())(p_index);

    temp = (temp - temperature_offset) * temperature_scale;
    return temp;
}

bool NanovdbMedium::SampleT_maj(const RayDifferential &_ray, Float u_t,
                                Float u_channel, MemoryArena &arena,
                                MajorantSampleRecord *maj_record) const {
    // \sum dt * sigma_maj = - log(1 - u)

    Ray ray = _ray;
    ray.tMax *= ray.d.Length();
    ray.d = Normalize(ray.d);

    DDATracker tracker = get_dda_tracker(ray);

    Float thick_bound = -std::log(1 - u_t);
    Spectrum sum(.0f);
    Float t_w = tracker.get_world_t();
    Float maj_density;

    int channel = SampleChannel(u_channel);
    maj_record->channel = channel;

    MajorantSeg seg;

    while (tracker.track(&seg)) {
        int ix = seg.index[0], iy = seg.index[1], iz = seg.index[2];
        maj_density = maj_grid->at(ix, iy, iz);
        Float dt = seg.dt_w;
        Spectrum sigma_maj = maj_density * sigma_t;

        if (sum[channel] + sigma_maj[channel] * seg.dt_w > thick_bound) {
            dt = (thick_bound - sum[channel]) / sigma_maj[channel];

            //* Sample the information at sampled point
            t_w += dt;
            sum += dt * maj_density * sigma_t;

            maj_record->T_maj = Exp(-sum);
            maj_record->p = ray(t_w);
            Float density = sampleDensity(maj_record->p);
            //            maj_record->Le = Spectrum(.0f);
            maj_record->Le = Le(maj_record->p);
            maj_record->phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
            maj_record->sigma_a = density * sigma_a;
            maj_record->sigma_s = density * sigma_s;
            maj_record->sigma_n = (maj_density - density) * sigma_t;
            maj_record->t = t_w;
            return false;  // not escape from current medium
        }

        sum += dt * sigma_maj;
        t_w += dt;

        tracker.next();
    }

    maj_record->T_maj = Exp(-sum);
    return true;  // escape the current medium
}

RegularTracker NanovdbMedium::get_regular_tracker(Ray ray_world) const {
    // Transform the ray in world space to grid index space
    Point3f o_grid = worldToIndex(ray_world.o);
    Vector3f d_grid = worldToIndex(ray_world.d);
    Float scalation = 1.0 / d_grid.Length();
    Float tmax_grid = ray_world.tMax / scalation;
    d_grid = Normalize(d_grid);

    // Reset the ray origin at the boundary at grid(if intersect)

    // Check if ray intersect the grid, if not, set terminate = true
    Float t_min = -Infinity, t_max = Infinity;

    for (int axis = 0; axis < 3; ++axis) {
        if (d_grid[axis] == 0) continue;

        Float t_0 = (minIndex[axis] - o_grid[axis]) / d_grid[axis],
              t_1 = (maxIndex[axis] + 1 - o_grid[axis]) / d_grid[axis];
        if (t_0 > t_1) std::swap(t_0, t_1);
        t_min = std::max(t_0, t_min);
        t_max = std::min(t_1, t_max);

        if (t_min > t_max || t_max < 0) {
            return RegularTracker(true);  // Just terminate the tracking
        }
    }

    Float scaled_voxel_size = scalation;
    Float cur_t = t_min > 0 ? t_min : 0;
    o_grid += cur_t * d_grid;

    return RegularTracker(minIndex, maxIndex, o_grid, d_grid, cur_t, tmax_grid,
                          scaled_voxel_size);
}

DDATracker NanovdbMedium::get_dda_tracker(Ray ray_world) const {
    // Check if ray_world intersect the maj_grid
    Point3f pmin = maj_grid->world_bound.pMin,
            pmax = maj_grid->world_bound.pMax;
    Float t_min = -Infinity, t_max = Infinity;
    Point3f origin = ray_world.o;
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

        if (t_min > t_max || t_max < 0)
            return DDATracker(true);  // Just terminate
    }

    Float cur_t_w = t_min > 0 ? t_min : 0;
    Point3f pIndex = maj_grid->toIndex(origin + direction * cur_t_w);

    return DDATracker(maj_grid->resolution, pIndex, direction,
                      maj_grid->voxel_size_w, cur_t_w, ray_world.tMax);
}

Point3f NanovdbMedium::worldToIndex(Point3f p_world) const {
    p_world = Inverse(medium_transform)(p_world);
    auto p_index = densityFloatGrid->worldToIndexF(
        nanovdb::Vec3f(p_world[0], p_world[1], p_world[2]));
    return Point3f{p_index[0], p_index[1], p_index[2]};
}

Vector3f NanovdbMedium::worldToIndex(Vector3f d_world) const {
    d_world = Inverse(medium_transform)(d_world);
    auto d_index = densityFloatGrid->worldToIndexDirF(
        nanovdb::Vec3f(d_world[0], d_world[1], d_world[2]));
    return Vector3f{d_index[0], d_index[1], d_index[2]};
}

Point3f NanovdbMedium::indexToWorld(Point3f p_index) const {
    auto p_world = densityFloatGrid->indexToWorldF(
        nanovdb::Vec3f(p_index[0], p_index[1], p_index[2]));
    return medium_transform(Point3f(p_world[0], p_world[1], p_world[2]));
}

Spectrum NanovdbMedium::Le(Point3f p_world) const {
    Float temperature = sampleTemperature(p_world);

    // TODO temperature operation
    if (temperature <= 100.0) return Spectrum(.0);

    // Compute blackbody emission at 12 sampled wavelength
    SampledSpectrum blackbody_emission_spectrum;

    constexpr int N_wavelengths = 12;
    constexpr Float sampled_wavelength[]{400.0,  427.27, 454.54, 481.82,
                                         509.10, 536.36, 563.64, 590.91,
                                         618.18, 645.45, 672.72, 700.0};
    Float Le_lambda[12];

    BlackbodyNormalized(sampled_wavelength, N_wavelengths, temperature,
                        Le_lambda);

    blackbody_emission_spectrum = SampledSpectrum::FromSampled(
        sampled_wavelength, Le_lambda, N_wavelengths);

    auto rgb = blackbody_emission_spectrum.ToRGBSpectrum();
    for (int i = 0; i < 3; ++i) rgb[i] = std::max(.0f, rgb[i]);
    return rgb * Le_scale;
}

Point3f EmissionGrid::sampleVoxel(Float u, Float *pdf) const {
    int idx = voxel_distrib->SampleDiscrete(u, pdf);
    return Point3f(voxels[idx][0], voxels[idx][1], voxels[idx][2]);
}

}  // namespace pbrt