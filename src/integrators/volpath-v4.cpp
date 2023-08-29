#include "integrators/volpath-v4.h"

#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

namespace pbrt {
void VolPathIntegratorV4::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum VolPathIntegratorV4::Li(const RayDifferential &r, const Scene &scene,
                                 Sampler &sampler, MemoryArena &arena,
                                 int depth) const {
    Spectrum L(.0f), beta(1.f), r_u(1.f), r_l(1.f);
    bool specular_bounce = false, any_nonspecular_bounce = false;
    Float eta_scale = 1;

    int bounces;
    RayDifferential ray(r);
    for (bounces = 0;; bounces++) {
        SurfaceInteraction isect;
        bool found_intersection = scene.Intersect(ray, &isect);

        if (ray.medium) {
            bool scattered = false, terminated = false;
            Float tmax = found_intersection ? ray.tMax : Infinity;

            RNG rng(rand());

            Spectrum T_maj = SampleT_maj(
                ray, tmax, rng, arena, [&](MajorantSampleRecord maj_rec) {
                    // Medium properties
                    Spectrum T_maj = maj_rec.T_maj, sigma_a = maj_rec.sigma_a,
                             sigma_s = maj_rec.sigma_s,
                             sigma_n = maj_rec.sigma_n,
                             sigma_maj = sigma_a + sigma_s + sigma_n,
                             Le = maj_rec.Le;
                    Point3f p = maj_rec.p;

                    if (beta.IsBlack()) {
                        terminated = true;
                        return false;
                    }

                    if (bounces < maxDepth && !maj_rec.Le.IsBlack()) {
                        Float pdf = T_maj[0] * sigma_maj[0];
                        Spectrum betap = beta * T_maj / pdf;

                        if (!betap.IsBlack()) L += betap * sigma_a * Le;
                    }

                    Float p_absorb = sigma_a[0] / sigma_maj[0];
                    Float p_scatter = sigma_s[0] / sigma_maj[0];
                    // Float p_null = std::max<Float>(0, 1 - p_absorb -
                    // p_scatter);

                    Float um = rng.UniformFloat();

                    auto SampleScatterType = [](Float um, Float p_absorb,
                                                Float p_scatter) {
                        if (um < p_absorb)
                            return 0;
                        else if (um < p_absorb + p_scatter)
                            return 1;
                        return 2;
                    };

                    int mode = SampleScatterType(um, p_absorb, p_scatter);

                    if (mode == 0 /* Absorb */) {
                        terminated = true;
                        return false;
                    } else if (mode == 1 /* Scatter */) {
                        if (depth++ > maxDepth) {
                            terminated = true;
                            return false;
                        }

                        Float pdf = T_maj[0] * sigma_s[0];
                        beta *= T_maj * sigma_s / pdf;
                        r_u *= T_maj * sigma_s / pdf;

                        if (!beta.IsBlack()) {
                            // handle real scatter
                            const Distribution1D *lightDistrib =
                                lightDistribution->Lookup(p);
                            // TODO
                        }

                        return false;
                    } else {
                        // TODO
                    }
                });

            if (terminated || beta.IsBlack()) return L;
            if (scattered) continue;

            beta *= T_maj / T_maj[0];
        }

        if (bounces == 0 || specular_bounce) {
            if (found_intersection)
                L += beta * isect.Le(-ray.d);
            else
                for (const auto &light : scene.infiniteLights)
                    L += beta * light->Le(ray);
        }

        if (!found_intersection || bounces >= maxDepth) break;

        isect.ComputeScatteringFunctions(ray, arena, true);
        if (!isect.bsdf) {
            ray = isect.SpawnRay(ray.d);
            bounces--;
            continue;
        }

        const Distribution1D *lightDistrib = lightDistribution->Lookup(isect.p);
        L += beta * UniformSampleOneLight(isect, scene, arena, sampler, true,
                                          lightDistrib);

        // Sampling volumetric emission

        Vector3f wo = -ray.d, wi;
        Float pdf;
        BxDFType flags;
        Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                          BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == .0f) break;
        beta *= f * AbsDot(wi, isect.shading.n) / pdf;
        specular_bounce = (flags & BSDF_SPECULAR) != 0;
        if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
            Float eta = isect.bsdf->eta;
            // Update the term that tracks radiance scaling for refraction
            // depending on whether the ray is entering or leaving the
            // medium.
            eta_scale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
        }
        ray = isect.SpawnRay(wi);

        // TODO Handle bssrdf

        Spectrum rrBeta = beta * eta_scale;
        if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
            Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
            if (sampler.Get1D() < q) break;
            beta /= 1 - q;
            DCHECK(std::isinf(beta.y()) == false);
        }
    }
    return L;
}

VolPathIntegratorV4 *CreateVolPathIntegratorV4(
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
    return new VolPathIntegratorV4(maxDepth, camera, sampler, pixelBounds,
                                   rrThreshold, lightStrategy);
}
}  // namespace pbrt