#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_VOLUMEPHOTONMAP_H
#define PBRT_INTEGRATORS_VOLUMEPHOTONMAP_H

// integrators/photonmap.h*
#include "pbrt.h"
#include "volume.h"
#include "integrator.h"
#include "core/vkdtree.h"

struct Photon;
struct RadiancePhoton;



struct PhotonProcess;
struct RadiancePhotonProcess;

struct VPhotonProcess;


// PhotonIntegrator Declarations
class VolumePhotonIntegrator : public VolumeIntegrator {

public:
    // PhotonIntegrator Public Methods
    
    VolumePhotonIntegrator(int nvol, int nl, int mphodepth, float mdist,
                           int gs, float ga, float mstep, float steps, float blur);
    
    ~VolumePhotonIntegrator(){}
    
    Spectrum Li(const Scene *scene,
                const Renderer *renderer,
                const RayDifferential &ray,
                const Sample *sample,
                RNG &rng,
                Spectrum *transmittance,
                MemoryArena &arena) const;
    
    Spectrum Transmittance(const Scene *scene,
                           const Renderer *,
                           const RayDifferential &ray,
                           const Sample *sample,
                           RNG &rng,
                           MemoryArena &arena) const;
    
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
    
private:
    // PhotonIntegrator Private Methods
    friend class VolumePhotonShootingTask;
    
    // PhotonIntegrator Private Data
    uint32_t nCausticPhotonsWanted, nIndirectPhotonsWanted, nVolumePhotonsWanted, nLookup;
    float maxDistSquared;
    int maxPhotonDepth;
    int gatherSamples;
    float cosGatherAngle;
    float marchStep;
    float stepSize;
    int tauSampleOffset, scatterSampleOffset;
    float blur; //controls radius size for beam photon
    
    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    BSDFSampleOffsets bsdfGatherSampleOffsets, indirGatherSampleOffsets;
    int nVolumePaths;

    VKdTree<Photon> *volumeMap;
};

VolumePhotonIntegrator *CreatePhotonMapVolumeIntegrator(const ParamSet &params);


#endif // PBRT_INTEGRATORS_PHOTONMAP_H
