#include "volintegrators.h"

#include "scene.h"

namespace pbrt {

int SampleCollisionType(Float u, Float p_absorb, Float p_scatter) {
    if (u < p_absorb)
        return 0;  // Sample absorb
    else if (u < p_absorb + p_scatter)
        return 1;  // Sample scatter
    return 2;      // Sample null scatter
}

std::shared_ptr<Medium> UniformChooseOneEmissionMedium(
    Float u, const std::vector<std::shared_ptr<Medium>> &emission_mediums,
    Float *pdf) {
    int N = emission_mediums.size();
    if (N == 0) return nullptr;
    *pdf = 1.0 / N;
    return emission_mediums[u * N];
}

Vector3f SampleDirectionTowardsEmitter(std::shared_ptr<Medium> emitter,
                                       const Interaction &it, Vector2f u2,
                                       Float *p_dir) {
    //
}

/*
-------------------------------------------------------------------
*/

VolPathIntegrator_LineIntegration::VolPathIntegrator_LineIntegration(
    int maxDepth, std::shared_ptr<const Camera> camera,
    std::shared_ptr<Sampler> sampler,
    std::vector<std::shared_ptr<Medium>> emission_mediums,
    const Bounds2i &pixelBounds, bool sample_le, Float rrThreshold,
    const std::string &lightSampleStrategy)
    : SamplerIntegrator(camera, sampler, pixelBounds),
      maxDepth(maxDepth),
      rrThreshold(rrThreshold),
      lightSampleStrategy(lightSampleStrategy),
      emission_mediums(emission_mediums),
      sample_le(sample_le) {}

void VolPathIntegrator_LineIntegration::Preprocess(const Scene &scene,
                                                   Sampler &sampler) {
    // do nothing
}

// Do ratio tracking like line integration to gather the le along given ray
Spectrum VolPathIntegrator_LineIntegration::IntegrateLe(
    const Scene &scene, const RayDifferential &_ray, Sampler &sampler) const {
    //
}

Spectrum VolPathIntegrator_LineIntegration::Li(const RayDifferential &_ray,
                                               const Scene &scene,
                                               Sampler &sampler,
                                               MemoryArena &arena,
                                               int depth) const {
    constexpr int channel = 0;

    Spectrum L(.0f), beta(1.f);
    RayDifferential ray(_ray);

    int bounces = 0;

    // Start tracing the ray
    while (true) {
        SurfaceInteraction isect;
        bool found_intersection = scene.Intersect(ray, &isect);

        const Medium *medium = ray.medium;
        if (medium) {
            // Tracing the medium

            bool scattered = false, terminated = false;
            Float t_max = found_intersection ? ray.tMax : Infinity;
            RNG rng(rand());

            // Cuz all volumetric emission term will be handled in IntegrateLe,
            // So, just tracing in this function (construct the path)
            Spectrum T_maj = SampleT_maj(
                ray, t_max, rng, arena, [&](MajorantSampleRecord maj_rec) {
                    Spectrum T_maj = maj_rec.T_maj;
                    Spectrum sigma_a = maj_rec.sigma_a;
                    Spectrum sigma_n = maj_rec.sigma_n;
                    Spectrum sigma_s = maj_rec.sigma_s;
                    Spectrum sigma_maj = sigma_a + sigma_s + sigma_n;
                    Point3f position = maj_rec.p;

                    if (beta.IsBlack()) {
                        terminated = true;
                        return false;  // No needs to continue sampling T_maj
                    }

                    // TODO no chromatic medium support
                    Float p_absorb = sigma_a[channel] / sigma_maj[channel];
                    Float p_scatter = sigma_s[channel] / sigma_maj[channel];

                    int mode = SampleCollisionType(rng.UniformFloat(), p_absorb,
                                                   p_scatter);

                    switch (mode) {
                    case 0 /* Absorb, just terminate the tracing*/:
                        terminated = true;
                        return false;
                    case 1 /* Scatter, update the ray and continue the tracing*/:
                        Float pdf = T_maj[channel] * sigma_s[channel];
                        beta *= T_maj * sigma_s / pdf;

                        if (++bounces > maxDepth || beta.IsBlack()) {
                            terminated = true;
                            return false;
                        }

                        MediumInteraction mi(position, -ray.d, ray.time, medium,
                                             maj_rec.phase);

                        // Sample a direction towards emission medium, and
                        // integrate on it
                        Float p_emitter;
                        auto emitter = UniformChooseOneEmissionMedium(
                            rng.UniformFloat(), emission_mediums, &p_emitter)
                    }
                });
        }
    }
}

}  // namespace pbrt