#pragma once
#include <medium.h>
#include <mutex>
namespace pbrt
{

class NanovdbMedium;

struct GridVertex
{
  Point3f  vertex_position;
  Spectrum intensity = Spectrum(.0f);

  Point3f  illumination_center;
  Vector3f illumination_offset        = Vector3f(.0f, .0f, .0f);
  Float    illumination_offset_weight = .0f;

  std::mutex mtx;
};

class LightGrid
{
  public:
  LightGrid() = default;

  LightGrid(Float voxel_size, Vector3i resolution, Point3f p_min);

  GridVertex& at(int x, int y, int z);

  const GridVertex& at(int x, int y, int z) const;

  Float                   voxel_size;
  Vector3i                resolution;
  std::vector<GridVertex> vertices;
};

class LightGridHierarchy
{
  public:
  LightGridHierarchy(int N_hierarchies);

  int                                     N_hierarchies;
  std::vector<std::unique_ptr<LightGrid>> grids;

  //* Friend declarations
  friend std::unique_ptr<LightGridHierarchy> CreateLGH(std::shared_ptr<NanovdbMedium> media,
                                                       int N_hierarchies);
};

std::unique_ptr<LightGridHierarchy> CreateLGH(std::shared_ptr<NanovdbMedium> media,
                                              int                            N_hierarchies);

} // namespace pbrt
