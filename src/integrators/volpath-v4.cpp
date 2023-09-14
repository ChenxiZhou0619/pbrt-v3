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

    Float prev_scatter_pdf;
    Point3f prev_shading_p;

    RayDifferential ray(r);

    int bounces = 0;
    while (true) {
        SurfaceInteraction isect;
        bool found_intersection = scene.Intersect(ray, &isect);

        if (ray.medium) {
            bool scattered = false, terminated = false;
            Float tmax = found_intersection ? ray.tMax : Infinity;

            RNG rng(rand());

            Float pdf_t = 1.f;

            Spectrum T_maj = SampleT_maj(
                ray, tmax, rng, arena, [&](MajorantSampleRecord maj_rec) {
                    // Medium properties
                    Spectrum T_maj = maj_rec.T_maj, sigma_a = maj_rec.sigma_a,
                             sigma_s = maj_rec.sigma_s,
                             sigma_n = maj_rec.sigma_n,
                             sigma_maj = sigma_a + sigma_s + sigma_n,
                             Le = maj_rec.Le;
                    Point3f p = maj_rec.p;
                    int channel = maj_rec.channel;

                    //                    L = Spectrum(1.f);
                    //                    terminated = true;
                    //                    return false;

                    if (beta.IsBlack()) {
                        terminated = true;
                        return false;
                    }

                    //                    if (bounces <= maxDepth &&
                    //                    !maj_rec.Le.IsBlack()) {
                    //                        Float pdf = T_maj[channel] *
                    //                        sigma_maj[channel]; Spectrum betap
                    //                        = beta * T_maj / pdf,
                    //                                 r_e = sigma_maj * T_maj /
                    //                                 pdf;
                    //
                    //                        if (!sampleLe || bounces == 0 ||
                    //                        specular_bounce)
                    //                            L += betap * sigma_a * Le;
                    //                    }

                    Float p_absorb = sigma_a[channel] / sigma_maj[channel];
                    Float p_scatter = sigma_s[channel] / sigma_maj[channel];

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
                        if (bounces <= maxDepth && !maj_rec.Le.IsBlack()) {
                            Float pdf = T_maj[channel] * sigma_a[channel];

                            pdf_t *= pdf;  //!

                            Spectrum betap = beta * T_maj / pdf,
                                     r_e = sigma_maj * T_maj / pdf;
                            Float misw;
                            if (!sampleLe || bounces == 0 || specular_bounce)
                                misw = 1.f;
                            else {
                                // compute misw
                                // misw = p_u / (p_u +p_l)
                                Float pdf_u = prev_scatter_pdf * pdf_t;
                                Float pdf_l =
                                    ray.medium->pdf_emissionP(maj_rec.p) *
                                    (maj_rec.p - prev_shading_p)
                                        .LengthSquared();
                                misw = pdf_u / (pdf_u + pdf_l);
                            }

                            L += betap * sigma_a * Le * misw;
                        }

                        terminated = true;
                        return false;

                    } else if (mode == 1 /* Scatter */) {
                        if (bounces++ > maxDepth) {
                            terminated = true;
                            return false;
                        }
                        Float pdf = T_maj[channel] * sigma_s[channel];

                        pdf_t *= pdf;  //!

                        Spectrum weight = T_maj * sigma_s / pdf;
                        beta *= T_maj * sigma_s / pdf / AverageRGB(weight);

                        if (!beta.IsBlack()) {
                            // handle real scatter
                            const Distribution1D *lightDistrib =
                                lightDistribution->Lookup(p);
                            MediumInteraction mi(maj_rec.p, -ray.d, ray.time,
                                                 ray.medium, maj_rec.phase);
                            L += beta * UniformSampleOneLight(mi, scene, arena,
                                                              sampler, true,
                                                              lightDistrib);
                            Vector3f wo = -ray.d, wi;
                            mi.phase->Sample_p(wo, &wi, sampler.Get2D());
                            ray = mi.SpawnRay(wi);
                            specular_bounce = false;
                            scattered = true;
                        }

                        return false;

                    } else {
                        Float pdf = (T_maj * sigma_n)[channel];

                        pdf_t += pdf;  //!

                        beta *= T_maj * sigma_n / pdf;
                        if (pdf == 0) beta = Spectrum(.0f);
                        return !beta.IsBlack();
                    }
                });

            if (terminated || beta.IsBlack()) break;
            if (scattered) continue;

            //            beta *= T_maj / AverageRGB(T_maj);
            beta *= T_maj / AverageRGB(T_maj);
        }

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

        const Distribution1D *lightDistrib = lightDistribution->Lookup(isect.p);
        L += beta * UniformSampleOneLight(isect, scene, arena, sampler, true,
                                          lightDistrib);

        if (sampleLe) {
            Float pdf_em;
            auto em = SelectEmissionVolume(sampler.Get1D(), &pdf_em);
            L += beta * SampleEmissionVolume(isect, scene, arena, sampler, em) /
                 pdf_em;
        }

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
    std::shared_ptr<const Camera> camera,
    std::vector<std::shared_ptr<Medium>> emissionMediums) {
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
    bool sampleLe = params.FindOneBool("sampleLe", false);

    int M = params.FindOneInt("M", 1);

    return new VolPathIntegratorV4(maxDepth, camera, sampler, emissionMediums,
                                   pixelBounds, sampleLe, M, rrThreshold,
                                   lightStrategy);
}

Spectrum VolPathIntegratorV4::SampleEmissionVolume(
    const Interaction &it, const Scene &scene, MemoryArena &arena,
    Sampler &sampler, std::shared_ptr<Medium> emission_medium) const {
    if (!sampleLe || !emission_medium) return Spectrum(.0);

    VolumetricEmissionPoint vep;
    Float u = sampler.Get1D();
    Vector3f u3{sampler.Get1D(), sampler.Get1D(), sampler.Get1D()};

    bool sampled = emission_medium->SampleEmissionPoint(u, u3, &vep);

    if (!sampled) return Spectrum(.0);

    MediumInteraction mi;
    mi.p = vep.p;

    VisibilityTester vt(it, mi);
    Spectrum tr = vt.Tr(scene, sampler);

    Float r_sqr = (vep.p - it.p).LengthSquared();

    Vector3f wo = it.wo, wi = Normalize(vep.p - it.p);

    Spectrum f;
    Float pdf_w;
    if (it.IsSurfaceInteraction()) {
        const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
        f = isect.bsdf->f(wo, wi) * AbsDot(wi, isect.shading.n);
        pdf_w = isect.bsdf->Pdf(wo, wi);

    } else {
        const MediumInteraction &mi = (const MediumInteraction &)it;
        f = mi.phase->p(wo, wi);
        pdf_w = mi.phase->p(wo, wi);
    }

    Spectrum Ld = f * tr * vep.sigma_a * vep.Le / (vep.pdf * r_sqr);

    // TODO misw
    // misw = p_l / (p_l + p_u)

    Float p_l = vep.pdf * r_sqr, p_u = pdf_w * tr[0] * vep.sigma_a[0];

    return Ld * p_l / (p_l + p_u);
}

Spectrum VolPathIntegratorV4::SampleEmissionVolumeRIS(
    const Interaction &it, const Scene &scene, MemoryArena &arena,
    Sampler &sampler, std::shared_ptr<Medium> emission_medium, int M) const {
    //

    Spectrum selected_tr;
    Point3f selected_p;
    Float sum_w = .0;
    Spectrum sum_w_spectrual(.0);

    for (int i = 0; i < M; ++i) {
        VolumetricEmissionPoint vep;
        emission_medium->SampleEmissionPoint(
            sampler.Get1D(),
            {sampler.Get1D(), sampler.Get1D(), sampler.Get1D()}, &vep);

        MediumInteraction mi;
        mi.p = vep.p;
        VisibilityTester vt(it, mi);
        Spectrum tr_approx = vt.Tr_coarse(scene, sampler);
        //        Spectrum tr_approx = vt.Tr(scene, sampler);

        Float r_sqr = (vep.p - it.p).LengthSquared();
        Vector3f wo = it.wo, wi = Normalize(vep.p - it.p);

        Spectrum f;
        if (it.IsSurfaceInteraction()) {
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->f(wo, wi) * AbsDot(wi, isect.shading.n);

        } else {
            const MediumInteraction &mi = (const MediumInteraction &)it;
            f = mi.phase->p(wo, wi);
        }

        // TODO MIS
        Spectrum w = f * tr_approx * vep.sigma_a * vep.Le / (vep.pdf * r_sqr);

        sum_w += AverageRGB(w);
        sum_w_spectrual += w;

        if (i == 0) {
            selected_tr = tr_approx;
            selected_p = vep.p;
        } else {
            if (sampler.Get1D() < AverageRGB(w) / sum_w) {
                // Replace
                selected_tr = tr_approx;
                selected_p = vep.p;
            }
        }
    }

    // Compute real tau
    MediumInteraction mi;
    mi.p = selected_p;
    VisibilityTester vt(it, mi);
    Spectrum tr = vt.Tr(scene, sampler);

    Spectrum up = tr * sum_w_spectrual, down = selected_tr * M;

    if (!up.IsBlack() && !down.IsBlack()) return up / down;
    //        return sum_w_spectrual;
    return Spectrum(.0f);
}

std::shared_ptr<Medium> VolPathIntegratorV4::SelectEmissionVolume(
    Float u, Float *pdf) const {
    int n = emissionMediums.size();
    *pdf = 1.f / n;
    return emissionMediums[u * n];
}
}  // namespace pbrt