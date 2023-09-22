#pragma once
#include "integrator.h"
#include "lightdistrib.h"
#include "pbrt.h"

namespace pbrt {
// Helper functions for volumetric path tracer
int SampleCollisionType(Float u, Float p_absorb, Float p_scatter);

std::shared_ptr<Medium> UniformChooseOneEmissionMedium(
    Float u, const std::vector<std::shared_ptr<Medium>> &emission_mediums,
    Float *pdf);

Vector3f SampleDirectionTowardsEmitter(std::shared_ptr<Medium> emitter,
                                       const Interaction &it, Vector2f u2,
                                       Float *p_dir);

/*
-------------------------------------------------------------------
*/

// Weighted-delta tracking volumetric path tracer
class VolpathIntegrator_WeightedDelta : public SamplerIntegrator {
  public:
    VolpathIntegrator_WeightedDelta(
        int maxDepth, std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler, const Bounds2i &pixelBounds,
        Float rrThreshold = 1,
        const std::string &lightSampleStrategy = "sptial");

    void Preprocess(const Scene &scene, Sampler &sampler);

    virtual Spectrum Li(const RayDifferential &ray, const Scene &scene,
                        Sampler &sampler, MemoryArena &arena, int depth) const;

  protected:
  protected:
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;
};

VolpathIntegrator_WeightedDelta *CreateVolPathIntegrator_WeightedDelta(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

// This integrator just trate emission volume as classic emitter
// However, when computing direct lighting, we gather the volumetric
// emission along the shadowray (when ratio tracking)
class VolPathIntegrator_LineIntegration : public SamplerIntegrator {
  public:
    VolPathIntegrator_LineIntegration(
        int maxDepth, std::shared_ptr<const Camera> camera,
        std::shared_ptr<Sampler> sampler,
        std::vector<std::shared_ptr<Medium>> emission_mediums,
        const Bounds2i &pixelBounds, bool sample_le, Float rrThreshold = 1,
        const std::string &lightSampleStrategy = "spatial");

    void Preprocess(const Scene &scene, Sampler &sampler);

    Spectrum IntegrateLe(const Scene &scene, const RayDifferential &ray,
                         Sampler &sampler, MemoryArena &arena) const;

    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;

    bool sample_le;
    std::vector<std::shared_ptr<Medium>> emission_mediums;
};

VolPathIntegrator_LineIntegration *CreateVolPathIntegrator_LineIntegration(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera,
    std::vector<std::shared_ptr<Medium>> emission_mediums);
}  // namespace pbrt