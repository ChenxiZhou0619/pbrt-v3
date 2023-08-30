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

Spectrum VolPathIntegratorV4::SampleVolumeLd(const Interaction &it,
                                             const Scene &scene,
                                             MemoryArena &arena,
                                             Sampler &sampler) const {
    const Point3f center{0, 0.5, 0};
    constexpr Float radius = 0.2;

    // random sample a point in sphere (respective to volume)
    auto uniform_sample_in_sphere = [&](Point3f shading_p, Point3f _center,
                                        Float _radius, Float *_pdf) -> Point3f {
        Float r_u = sampler.Get1D();
        Float phi_u = sampler.Get1D(), theta_u = sampler.Get1D();

        Float r = _radius * std::pow(r_u, 1.0 / 3.0);
        Float phi = 2 * M_PI * phi_u;
        Float cos_theta = 1 - 2 * theta_u,
              sin_theta = std::sqrt(std::max(1.0 - cos_theta * cos_theta, .0));

        Vector3f p_local{sin_theta * r * std::cos(phi), cos_theta * r,
                         sin_theta * r * std::sin(phi)};
        Point3f p_e = _center + p_local;

        // Compute the sampling pdf respective to solid angle
        Float dist_sqr = (p_e - shading_p).LengthSquared();
        *_pdf = 3.0 / (4 * M_PI * _radius * _radius * _radius) * dist_sqr;

        return p_e;
    };

    Float pdf_pe = .0;
    Point3f p_e = uniform_sample_in_sphere(it.p, center, radius, &pdf_pe);

    MediumInteraction mi;
    mi.p = p_e;
    VisibilityTester vis{it, mi};

    Vector3f wi = Normalize(p_e - it.p);

    Spectrum f, Tr;
    Float misw;

    if (it.IsSurfaceInteraction()) {
        const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
        f = isect.bsdf->f(isect.wo, wi) * AbsDot(wi, isect.shading.n);
        Tr = vis.Tr(scene, sampler);

        Float pdf_wi = isect.bsdf->Pdf(isect.wo, wi);
        Float pdf_t = Tr[0] * Spectrum(1.f)[0];
        Float pdf_1 = pdf_pe, pdf_2 = pdf_wi * pdf_t;
        misw = pdf_1 / (pdf_1 + pdf_2);
    } else {
        // TODO Volumetric scattering
    }

    // TODO Hard code here

    return f * Tr * Spectrum(1.f) * Spectrum(1.f) / pdf_pe * misw;
}

Spectrum VolPathIntegratorV4::Li(const RayDifferential &r, const Scene &scene,
                                 Sampler &sampler, MemoryArena &arena,
                                 int depth) const {
    Spectrum L(.0f), beta(1.f), r_u(1.f), r_l(1.f);
    bool specular_bounce = false, any_nonspecular_bounce = false;
    Float eta_scale = 1;

    int bounces;

    Float prev_scatter_pdf;
    Point3f prev_shading_p;

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

                        if (!betap.IsBlack()) {
                            if (bounces == 0 || specular_bounce)
                                L += betap * sigma_a * Le;
                            else {
                                // MIS
                                Float pdf_1 = prev_scatter_pdf * pdf;
                                Float pdf_2 =
                                    3 / (4 * M_PI * 0.2 * 0.2 * 0.2 *
                                         (p - prev_shading_p).LengthSquared());
                                Float misw = pdf_1 / (pdf_1 + pdf_2);
                                L += betap * sigma_a * Le * misw;
                            }
                        }
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

        L += beta * SampleVolumeLd(isect, scene, arena, sampler);

        // Sampling volumetric emission

        Vector3f wo = -ray.d, wi;
        Float pdf;
        BxDFType flags;
        Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                          BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == .0f) break;
        beta *= f * AbsDot(wi, isect.shading.n) / pdf;

        prev_scatter_pdf = pdf;
        prev_shading_p = isect.p;

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