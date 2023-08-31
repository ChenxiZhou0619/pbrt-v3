#include "nanovdbmedium.h"

#include "memory.h"

namespace pbrt {

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

Spectrum NanovdbMedium::Tr(const Ray &ray, Sampler &sampler) const {
    RegularTracker tracker = get_regular_tracker(ray);
    Spectrum thick(.0f);
    TrackSegment seg;
    while (tracker.track(&seg)) {
        Float density = sampleDensity(seg.voxelIndex);
        thick += density * sigma_maj * seg.dt;
    }
    return Exp(-thick);
}

Spectrum NanovdbMedium::Sample(const Ray &ray, Sampler &sampler,
                               MemoryArena &arena,
                               MediumInteraction *mi) const {
    // TODO
}

Float NanovdbMedium::sampleDensity(int voxel_index[3]) const {
    using Sampler =
        nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
    nanovdb::Vec3R sample_loc{voxel_index[0] + .5, voxel_index[1] + .5,
                              voxel_index[2] + .5};
    return Sampler(densityFloatGrid->tree())(sample_loc) * density_scale;
}

bool NanovdbMedium::SampleT_maj(const RayDifferential &ray, Float u,
                                MemoryArena &arena,
                                MajorantSampleRecord *maj_record) const {
    // \sum dt * sigma_maj = - log(1 - u)
    RegularTracker tracker = get_regular_tracker(ray);
    Float thick_bound = -std::log(1 - u);

    bool sampled = false;
    Spectrum sum(.0f);
    Float t = tracker.get_world_t();

    Float density;

    TrackSegment seg;
    while (tracker.track(&seg)) {
        density = sampleDensity(seg.voxelIndex);
        Float maj = density * sigma_maj[0];  // TODO Chromatic
        Float dt = seg.dt;

        if (sum[0] + maj * seg.dt > thick_bound) {
            Float dt = (thick_bound - sum[0]) / maj;
            sampled = true;
        }

        t += dt;
        sum += dt * density * sigma_maj;

        if (sampled) break;
    }

    maj_record->T_maj = Exp(-sum);
    if (sampled) {
        maj_record->p = ray(t);
        maj_record->Le = Spectrum(.0f);
        maj_record->phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
        maj_record->sigma_a = density * sigma_a;
        maj_record->sigma_s = density * sigma_s;
        maj_record->sigma_n = Spectrum(.0f);  // TODO
        maj_record->t = t;
    }

    return !sampled;
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

}  // namespace pbrt