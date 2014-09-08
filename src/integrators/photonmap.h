#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PHOTONMAP_H
#define PBRT_INTEGRATORS_PHOTONMAP_H

// integrators/photonmap.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"

struct Photon;
struct RadiancePhoton;
struct ClosePhoton;
struct PhotonProcess;
struct RadiancePhotonProcess;


// PhotonIntegrator Declarations
class PhotonIntegrator : public SurfaceIntegrator {
public:
    // PhotonIntegrator Public Methods
    PhotonIntegrator(int ncaus, int nindir, int nLookup, int maxspecdepth,
        int maxphotondepth, float maxdist, bool finalGather, int gatherSamples,
        float ga);
    ~PhotonIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
private:
    // PhotonIntegrator Private Methods
    friend class PhotonShootingTask;

    // PhotonIntegrator Private Data
    uint32_t nCausticPhotonsWanted, nIndirectPhotonsWanted, nLookup;
    float maxDistSquared;
    int maxSpecularDepth, maxPhotonDepth;
    bool finalGather;
    int gatherSamples;
    float cosGatherAngle;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    BSDFSampleOffsets bsdfGatherSampleOffsets, indirGatherSampleOffsets;
    int nCausticPaths, nIndirectPaths;
    KdTree<Photon> *causticMap;
    KdTree<Photon> *indirectMap;
    KdTree<RadiancePhoton> *radianceMap;
};


PhotonIntegrator *CreatePhotonMapSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PHOTONMAP_H
