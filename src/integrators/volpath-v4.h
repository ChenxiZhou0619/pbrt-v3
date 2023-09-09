#pragma once

#include "integrator.h"
#include "lightdistrib.h"
#include "pbrt.h"

namespace pbrt {
class VolPathIntegratorV4 : public SamplerIntegrator {
  public:
    VolPathIntegratorV4(int maxDepth, std::shared_ptr<const Camera> camera,
                        std::shared_ptr<Sampler> sampler,
                        std::vector<std::shared_ptr<Medium>> emissionMediums,
                        const Bounds2i &pixelBounds, bool sampleLe,
                        Float rrThreshold = 1,
                        const std::string &lightSampleStrategy = "spatial")
        : SamplerIntegrator(camera, sampler, pixelBounds),
          maxDepth(maxDepth),
          rrThreshold(rrThreshold),
          lightSampleStrategy(lightSampleStrategy),
          emissionMediums(emissionMediums),
          sampleLe(sampleLe) {}

    void Preprocess(const Scene &scene, Sampler &sampler);

    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;
    std::vector<std::shared_ptr<Medium>> emissionMediums;

    bool sampleLe = false;
};

VolPathIntegratorV4 *CreateVolPathIntegratorV4(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera,
    std::vector<std::shared_ptr<Medium>> emissionMediums);
};  // namespace pbrt