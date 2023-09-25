#include <api.h>
#include <integrators/LGH/lightgrid.h>
#include <media/nanovdbmedium.h>
using namespace pbrt;

int main()
{
  pbrtInit(Options());

  std::string vdbfilename =
      "/home/chenxizhou/Master/Programming/pbrt-v3/scenes/explosion/gasoline.nvdb";
  Float       density_scale     = 1.f;
  Float       g                 = 0.;
  Spectrum    sigam_a           = 1.f;
  Spectrum    sigma_s           = 1.f;
  Transform   medium_transform  = Transform();
  std::string density_name      = "density";
  std::string temperature_name  = "temperature";
  Float       LeScale           = 8.f;
  Float       temperatureOffset = .0f;
  Float       temperatureScale  = 8000;

  auto nanovdb_media = std::make_shared<NanovdbMedium>(
      vdbfilename, density_scale, sigam_a, sigma_s, g, medium_transform, density_name,
      temperature_name, LeScale, temperatureOffset, temperatureScale, false);

  auto lgh = CreateLGH(nanovdb_media, 8);

  for (int level = 0; level < lgh->N_hierarchies; ++level)
  {
    Spectrum sum = .0f;

    for (const auto& vtx : lgh->grids[level]->vertices)
    {
      sum += vtx.intensity;
    }

    std::cout << "Level " << level << "'s sum :" << sum << std::endl;
  }

  for (const auto& vtx : lgh->grids[0]->vertices)
  {
    std::cout << "Vertex[0] position is : " << vtx.vertex_position;
    std::cout << ", illumination center is : " << vtx.illumination_center << std::endl;
  }

  pbrtCleanup();
  return 0;
}