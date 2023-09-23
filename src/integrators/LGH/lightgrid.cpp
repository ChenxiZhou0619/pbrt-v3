#include "lightgrid.h"

namespace pbrt {
LightGrid::LightGrid(Float voxel_size, Vector3i resolution)
    : voxel_size(voxel_size), resolution(resolution) {
    int N_vertex =
        (resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1);
    grid.reserve(N_vertex);
}

std::unique_ptr<LightGridHierarchy> CreateLGH(
    std::shared_ptr<NanovdbMedium> media) {
    //* Create the light grid hierarchy

    // get the base resolution
}

};  // namespace pbrt