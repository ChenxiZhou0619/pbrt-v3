#pragma once
#include <medium.h>
namespace pbrt {

class NanovdbMedium;

struct GridVertex {
    Point3f vertex_position;
    Point3f illumination_center;
    Spectrum intensity;
};

class LightGrid {
  public:
    LightGrid() = default;

    LightGrid(Float voxel_size, Vector3i resolution);

  private:
    Float voxel_size;
    Vector3i resolution;
    std::vector<GridVertex> grid;
};

class LightGridHierarchy {
  public:
    LightGridHierarchy() = default;

  private:
    std::vector<LightGrid> hierarchy;
};

std::unique_ptr<LightGridHierarchy> CreateLGH(
    std::shared_ptr<NanovdbMedium> media);

}  // namespace pbrt
