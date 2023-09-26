#pragma once

#include <medium.h>
#include <texture.h>
namespace pbrt
{

//* Used to store the accumulated density
struct FloatTexture2D
{
  public:
  FloatTexture2D(Vector2i resolution);

  Float at(Point2f uv) const;

  private:
  std::vector<Float> data;
};

#define X_NEGATIVE 0
#define X_POSITIVE 1
#define Y_NEGATIVE 2
#define Y_POSITIVE 3
#define Z_NEGATIVE 4
#define Z_POSITIVE 5

class NanovdbMedium;
class LightGridHierarchy;

// TODO pre-defined distance
class DeepShadowMap
{
  public:
  //* The cube-map contains resolution ^ 2 * 6 texels
  DeepShadowMap(int resolution, Float r_l, Point3f center);

  Spectrum query(Vector3f direction, Float distance) const;

  Vector3f texelToDirection(int face, int u, int v) const;

  public:
  int     resolution;
  Point3f center;

  private:
  std::unique_ptr<FloatTexture2D> cubemap_t0[6];
  std::unique_ptr<FloatTexture2D> cubemap_t1[6];
  std::unique_ptr<FloatTexture2D> cubemap_t2[6];
  std::unique_ptr<FloatTexture2D> cubemap_t3[6];

  // TODO
  void setCubeMap(Float accmulated_densities[4], int u, int v);

  friend void InitializeDeepShadowmap(const NanovdbMedium& media, LightGridHierarchy& lgh);
};
} // namespace pbrt
