#include "volintegrators.h"

#include "camera.h"
#include "paramset.h"
#include "scene.h"

namespace pbrt {

#define COLLISION_ABSORB 0
#define COLLISION_SCATTER 1
#define COLLISION_NULL 2

int SampleCollisionType(Float u, Float p_absorb, Float p_scatter) {
    if (u < p_absorb)
        return COLLISION_ABSORB;  // Sample absorb
    else if (u < p_absorb + p_scatter)
        return COLLISION_SCATTER;  // Sample scatter
    return COLLISION_NULL;         // Sample null scatter
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
-------------   VolpathIntegrator_WeightedDelta   --------------
*/

VolpathIntegrator_WeightedDelta::VolpathIntegrator_WeightedDelta(
    int maxDepth, std::shared_ptr<const Camera> camera,
    std::shared_ptr<Sampler> sampler, const Bounds2i &pixelBounds,
    Float rrThreshold, const std::string &lightSampleStrategy)
    : SamplerIntegrator(camera, sampler, pixelBounds),
      maxDepth(maxDepth),
      rrThreshold(rrThreshold),
      lightSampleStrategy(lightSampleStrategy) {}

void VolpathIntegrator_WeightedDelta::Preprocess(const Scene &scene,
                                                 Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum VolpathIntegrator_WeightedDelta::Li(const RayDifferential &_ray,
                                             const Scene &scene,
                                             Sampler &sampler,
                                             MemoryArena &arena,
                                             int depth) const {
    // TODO Support chromatic media
    constexpr int channel = 0;

    Spectrum L(.0f);            //* The accumulative radiance on sampled path
    Spectrum beta(1.f);         //* The weight for current path
    RayDifferential ray(_ray);  //* The ray to be traced
    int bounces = 0;  //* How many scatter collisions exist on current path

    bool specular_bounce = false;

    while (true) {
        if (beta.IsBlack()) break;

        SurfaceInteraction isect;
        bool found_intersection = scene.Intersect(ray, &isect);

        const Medium *medium = ray.medium;
        if (medium) {
            //* Start tracking if current ray starts in medium

            bool scattered = false;   //* Whether a scatter collision happens
            bool terminated = false;  //* Whether to terminate the ray (no
                                      //* contribution afterwards or absorbed)
            Float t_max = found_intersection ? ray.tMax : Infinity;

            RNG rng(rand());  //* Random numbers for stochastic tracking

            //* Weighted delta tracking may have negative sigma_n
            //* Return false to stop majorant tracking on current direction
            SampleT_maj(
                ray, t_max, rng, arena, [&](MajorantSampleRecord maj_rec) {
                    Point3f position = maj_rec.p;
                    Spectrum T_maj = maj_rec.T_maj;
                    Spectrum sigma_a = maj_rec.sigma_a;
                    Spectrum sigma_s = maj_rec.sigma_s;
                    Spectrum sigma_n = maj_rec.sigma_n;
                    Spectrum sigma_t = sigma_a + sigma_s;
                    Spectrum sigma_maj = sigma_t + sigma_n;
                    Spectrum Le = maj_rec.Le;

                    //^ update majorant tracking weight
                    Float pdf = T_maj[channel] * sigma_maj[channel];
                    beta *= T_maj * sigma_maj / pdf;

                    //* compute the collision posibility
                    Float p_absorb =
                        sigma_a[channel] /
                        (sigma_t[channel] + std::abs(sigma_n[channel]));
                    Float p_scatter =
                        sigma_s[channel] /
                        (sigma_t[channel] + std::abs(sigma_n[channel]));
                    Float p_null = 1.f - p_absorb - p_scatter;

                    int mode = SampleCollisionType(rng.UniformFloat(), p_absorb,
                                                   p_scatter);

                    switch (mode) {
                    case COLLISION_ABSORB:
                        //^ Accumulate the self-emission and terminate
                        beta *= sigma_a / (sigma_maj * p_absorb);
                        L += beta * Le;
                        terminated = true;
                        return false;
                    case COLLISION_SCATTER:
                        //^ Sample the direct lighting and update ray
                        beta *= sigma_s / (sigma_maj * p_scatter);

                        if (++bounces > maxDepth) {
                            terminated = true;
                            return false;
                        }

                        if (!beta.IsBlack()) {
                            //^ 1. Sample Ld (mis in UniformSampleOneLight)
                            const Distribution1D *lightDistrib =
                                lightDistribution->Lookup(position);
                            MediumInteraction mi(position, -ray.d, ray.time,
                                                 medium, maj_rec.phase);
                            L += beta * UniformSampleOneLight(mi, scene, arena,
                                                              sampler, true,
                                                              lightDistrib);

                            //^ 2. Continue the ray by phase sampling
                            Vector3f wo = -ray.d, wi;
                            mi.phase->Sample_p(wo, &wi, sampler.Get2D());
                            ray = mi.SpawnRay(wi);
                            scattered = true;
                        }
                        return false;
                    case COLLISION_NULL:
                        //^ Update beta and continue the majorant tracking

                        //^ Possible alternate the sign
                        beta *= sigma_n / (sigma_maj * p_null);

                        return true;
                    }
                    return false;
                });

            if (terminated || beta.IsBlack()) break;
            if (scattered) continue;

            //* The ray go through the medium without collision
        }

        //* Handle the surface or miss interaction
        if (bounces == 0 || specular_bounce) {
            if (found_intersection)
                L += beta * isect.Le(-ray.d);
            else
                for (const auto &light : scene.infiniteLights)
                    L += beta * light->Le(ray);
        }

        if (!found_intersection) break;

        isect.ComputeScatteringFunctions(ray, arena, true);
        if (!isect.bsdf) {
            ray = isect.SpawnRay(ray.d);
            continue;
        }

        if (++bounces > maxDepth) break;
        //^ Sample direct lighting
        const Distribution1D *lightDistrib = lightDistribution->Lookup(isect.p);
        L += beta * UniformSampleOneLight(isect, scene, arena, sampler, true,
                                          lightDistrib);
        //^ Bsdf sampling
        Vector3f wo = -ray.d, wi;
        Float pdf;
        BxDFType flags;
        Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                          BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == .0f) break;
        beta *= f * AbsDot(wi, isect.shading.n) / pdf;
        specular_bounce = (flags & BSDF_SPECULAR) != 0;
        ray = isect.SpawnRay(wi);

        Spectrum rrBeta = beta;
        if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
            Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
            if (sampler.Get1D() < q) break;
            beta /= 1 - q;
            DCHECK(std::isinf(beta.y()) == false);
        }
    }

    return L;
}

VolpathIntegrator_WeightedDelta *CreateVolPathIntegrator_WeightedDelta(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string lightStrategy =
        params.FindOneString("lightsamplestrategy", "spatial");

    return new VolpathIntegrator_WeightedDelta(
        maxDepth, camera, sampler, pixelBounds, rrThreshold, lightStrategy);
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
                            rng.UniformFloat(), emission_mediums, &p_emitter);
                    }
                });
        }
    }
}

}  // namespace pbrt