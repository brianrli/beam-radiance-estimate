#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_WHITTED_H
#define PBRT_INTEGRATORS_WHITTED_H

// integrators/whitted.h*
#include "pbrt.h"
#include "integrator.h"
#include "scene.h"

// WhittedIntegrator Declarations
class WhittedIntegrator : public SurfaceIntegrator {
public:
    // WhittedIntegrator Public Methods
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    WhittedIntegrator(int md) {
        maxDepth = md;
    }
private:
    // WhittedIntegrator Private Data
    int maxDepth;
};


WhittedIntegrator *CreateWhittedSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_WHITTED_H
