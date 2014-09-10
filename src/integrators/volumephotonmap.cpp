// integrators/photonmap.cpp*
#include "stdafx.h"
#include "integrators/volumephotonmap.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"

struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w)
    : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
};

struct RadiancePhoton {
    RadiancePhoton(const Point &pp, const Normal &nn)
    : p(pp), n(nn), Lo(0.f) { }
    RadiancePhoton() { }
    Point p;
    Normal n;
    Spectrum Lo;
};

class VolumePhotonShootingTask : public Task {
public:
    VolumePhotonShootingTask(int tn, float ti, Mutex &m, VolumePhotonIntegrator *in,
                       ProgressReporter &prog, bool &at, int &ndp,
                       vector<Photon> &direct, vector<Photon> &indir, vector<Photon> &caustic, vector<Photon> &vol,
                       vector<RadiancePhoton> &rps, vector<Spectrum> &rpR, vector<Spectrum> &rpT,
                       uint32_t &ns, Distribution1D *distrib, const Scene *sc,
                       const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
    abortTasks(at), nDirectPaths(ndp),
    directPhotons(direct), indirectPhotons(indir), causticPhotons(caustic),
    radiancePhotons(rps), rpReflectances(rpR), rpTransmittances(rpT),
    nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr),
    volumePhotons(vol){ }
    void Run();
    
    int taskNum;
    float time;
    Mutex &mutex;
    VolumePhotonIntegrator *integrator;
    ProgressReporter &progress;
    bool &abortTasks;
    int &nDirectPaths;
    vector<Photon> &directPhotons, &indirectPhotons, &causticPhotons, &volumePhotons;
    vector<RadiancePhoton> &radiancePhotons;
    vector<Spectrum> &rpReflectances, &rpTransmittances;
    uint32_t &nshot;
    const Distribution1D *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
};

inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}

//use Hayley-Greenstein function
Vector VSampleHG(const Vector &w, float g, float u1, float u2) {
    float costheta;
    if (fabsf(g) < 1e-3)
        costheta = 1.f - 2.f * u1;
    else {
        float sqrTerm = (1.f - g * g) /
        (1.f - g + 2.f * g * u1);
        costheta = (1.f + g * g - sqrTerm * sqrTerm) / (2.f * g);
    }
    float sintheta = sqrtf(max(0.f, 1.f-costheta*costheta));
    float phi = 2.f * M_PI * u2;
    Vector v1, v2;
    CoordinateSystem(w, &v1, &v2);
    return SphericalDirection(sintheta, costheta, phi, v1, v2, w);
}

float VPhaseHG(const Vector &w, const Vector &wp, float g) {
    float costheta = Dot(w, wp);
    return 1.f / (4.f * M_PI) *
    (1.f - g*g) / powf(1.f + g*g - 2.f * g * costheta, 1.5f);
}

float SampleScattering(const Vector &wi, float u1, float u2, Vector &wo){
    wo = VSampleHG(wi, 0.9f, u1, u2);
    return VPhaseHG(wi, wo, 0.9f);
}
void VolumePhotonShootingTask::Run() {
    // Declare local variables for _VolumePhotonShootingTask_
    VolumeRegion *volume = scene->volumeRegion;
    float marchstep = integrator->marchStep;

    MemoryArena arena;
    
    RNG rng(31 * taskNum);
    vector<Photon> localVolumePhotons; //new data structure
    
    uint32_t totalPaths = 0;
    
    //set finish conditions
    bool volumeDone = (integrator->nVolumePhotonsWanted == 0);
    
    PermutedHalton halton(6, rng);
    
    while (true) {

//        const uint32_t blockSize = 4096;
        const uint32_t blockSize = 500;
        for (uint32_t i = 0; i < blockSize; ++i) {
            
            float u[6];
            halton.Sample(++totalPaths, u); //random sample
            
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];
            
            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);

            if (!alpha.IsBlack()) {

                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);

                bool hitVolumeOnce = false;

                Intersection photonIsect;
                int nIntersections = 0;
                float vt0, vt1;
                Spectrum oldTransmittance(1.f);
                
                //if ray intersects with volume
                if (volume->IntersectP(photonRay, &vt0, &vt1)) {

                    ++nIntersections;
                    float inVolume = true;

                    // Enter volume photon mode till we leave the volumeRegion
                    hitVolumeOnce = true;
                    float u2, u6;
                    
                    Ray volumeRay;
                    volumeRay.o = photonRay(vt0);
                    volumeRay.d = Normalize(photonRay.d); // Hard to know if it's normalized. Just do it again.
                    volumeRay.mint = 0;
                    
                    u2 = (rng.RandomFloat() + 0.5f) * marchstep;
                    if( u2 >= (vt1-vt0)/photonRay.d.Length()){
                        photonRay.o = volumeRay(vt1+0.00001f);
                        alpha *= Exp(-volume->tau(volumeRay));
                        break;
                    }
                    else volumeRay.maxt = u2;

                    while(inVolume)
                    {
                    
                    // Ray march and decide to scattering/absorbing or move on
                        Spectrum transmittance = Exp(-volume->tau(volumeRay)) * oldTransmittance;
                        alpha *= transmittance;

                        if(1){ // Enforce interaction
                            Point interactPoint = volumeRay(volumeRay.maxt);
                            
                            // Store a photon
                            if(!volumeDone){
                                Photon photon(interactPoint, alpha, volumeRay.d);
                                localVolumePhotons.push_back(photon);
                                progress.Update();
                            }
                            
                            // Choose scattering or absorbing
                            if(rng.RandomFloat() * volume->sigma_t(interactPoint, volumeRay.d,volumeRay.time).y() <
                               volume->sigma_s(interactPoint, volumeRay.d,volumeRay.time).y()){

                                // Scattering, sampling a new direction
                                oldTransmittance = Spectrum(1.f);
                                Vector newDirection;
                                float pdf = SampleScattering(volumeRay.d, rng.RandomFloat(), rng.RandomFloat(), newDirection);
                                alpha *= pdf;

                                // Specify new volumeRay
                                photonRay.o = interactPoint;
                                photonRay.d = newDirection;
                                u6 = rng.RandomFloat() + 0.5f;
                                photonRay.maxt = marchstep * u6; //new photonray
                                
                                if(!volume->Inside(photonRay(photonRay.maxt))){ break; }
                                //++nIntersections;
    //                                volume->IntersectP(volumeRay, &vt0, &vt1);
    //                                volumeRay.maxt = vt1;
    //                                alpha *= Exp(-volume->tau(volumeRay));
    //                                photonRay = RayDifferential(volumeRay(vt1+0.0001f), volumeRay.d);
    //                                inVolume = false;
    //                            }
                            }
                            else
                                break; // Absorbing, end of the story
                            //if (nIntersections >= integrator->maxPhotonDepth) break;
                        }
                    }
//                    else alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
                }
                PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
            }
            arena.FreeAll();
        }
        
        // Merge local photon data with data in _PhotonIntegrator_
        { MutexLock lock(mutex);
            
            // Give up if we're not storing enough photons
            if (abortTasks)
                return;
            if (nshot > 500000 && (unsuccessful(integrator->nVolumePhotonsWanted,
                                                volumePhotons.size(), blockSize))){
                Error("Unable to store enough volume photons.  Giving up.\n");
                volumePhotons.erase(volumePhotons.begin(), volumePhotons.end());
                abortTasks = true;
                return;
            }
            
            //update progress
            progress.Update(localVolumePhotons.size());
            nshot += blockSize;
            
            //Merge Volume Photons into main
            if (!volumeDone) {
                integrator->nVolumePaths += blockSize;
                for (uint32_t i = 0; i < localVolumePhotons.size(); ++i)
                    volumePhotons.push_back(localVolumePhotons[i]);
                localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
                if (volumePhotons.size() >= integrator->nVolumePhotonsWanted)
                    volumeDone = true;
            }
            
        }
        
        // Exit task if enough photons have been found
        if (volumeDone)
            break;
    }
    
}

void VolumePhotonIntegrator::Preprocess(const Scene *scene,
                                  const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    // Declare shared variables for photon shooting
    Mutex *mutex = Mutex::Create();
    int nDirectPaths = 0;
    vector<Photon> causticPhotons, directPhotons, indirectPhotons, volumePhotons;
    vector<RadiancePhoton> radiancePhotons;
    bool abortTasks = false;
    
    causticPhotons.reserve(nCausticPhotonsWanted);
    indirectPhotons.reserve(nIndirectPhotonsWanted);
    volumePhotons.reserve(nVolumePhotonsWanted);
    
    uint32_t nshot = 0;
    vector<Spectrum> rpReflectances, rpTransmittances;
    
    // Compute light power CDF for photon shooting
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
    
    // Run parallel tasks for photon shooting
    ProgressReporter progress(nCausticPhotonsWanted+nIndirectPhotonsWanted, "Shooting photons");
    vector<Task *> VolumePhotonShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        VolumePhotonShootingTasks.push_back(new VolumePhotonShootingTask(
                                                             i, camera ? camera->shutterOpen : 0.f, *mutex, this, progress, abortTasks, nDirectPaths,
                                                             directPhotons, indirectPhotons, causticPhotons, volumePhotons, radiancePhotons, 
                                                             rpReflectances, rpTransmittances,
                                                             nshot, lightDistribution, scene, renderer));
    EnqueueTasks(VolumePhotonShootingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < VolumePhotonShootingTasks.size(); ++i)
        delete VolumePhotonShootingTasks[i];
    Mutex::Destroy(mutex);
    progress.Done();
    
    // Build kd-trees for indirect and caustic photons
    KdTree<Photon> *directMap = NULL;

    if (volumePhotons.size() > 0)
        volumeMap = new KdTree<Photon>(volumePhotons);

    delete directMap;
}

struct ClosePhoton {
    // ClosePhoton Public Methods
    ClosePhoton(const Photon *p = NULL, float md2 = INFINITY)
    : photon(p), distanceSquared(md2) { }
    bool operator<(const ClosePhoton &p2) const {
        return distanceSquared == p2.distanceSquared ?
        (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
    }
    const Photon *photon;
    float distanceSquared;
};

struct VPhotonProcess {
	// VPhotonProcess Public Methods
	VPhotonProcess(u_int mp, const Point &P)
    : p(P) {
		photons = 0;
		nLookup = mp;
		nfound = 0;
	}
    void operator()(const Point &p, const Photon &photon, float dist2, float &maxDistSquared){
        if (nfound < nLookup) {
			// Add photon to unordered array of photons
			photons[nfound++] = ClosePhoton(&photon, dist2);
			if (nfound == nLookup) {
				std::make_heap(&photons[0], &photons[nLookup]);
				maxDistSquared = photons[0].distanceSquared;
			}
		}
		else {
			// Remove most distant photon from heap and add new photon
			std::pop_heap(&photons[0], &photons[nLookup]);
			photons[nLookup-1] = ClosePhoton(&photon, dist2);
			std::push_heap(&photons[0], &photons[nLookup]);
			maxDistSquared = photons[0].distanceSquared;
		}
	}
	const Point &p;
	ClosePhoton *photons;
	u_int nLookup;
	mutable u_int nfound;
};

void VolumePhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                        const Scene *scene) {
	tauSampleOffset = sample->Add1D(1);
	scatterSampleOffset = sample->Add1D(1);
}

Spectrum VolumePhotonIntegrator::Transmittance(const Scene *scene,
                       const Renderer *,
                       const RayDifferential &ray,
                       const Sample *sample,
                       RNG &rng,
                       MemoryArena &arena) const {
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}

Spectrum VolumePhotonIntegrator::Li(const Scene *scene,
            const Renderer *renderer,
            const RayDifferential &ray,
            const Sample *sample,
            RNG &rng,
            Spectrum *transmittance,
            MemoryArena &arena) const {

	// Pointers assignment
	VolumeRegion *vr = scene->volumeRegion;
	KdTree<Photon> *map = volumeMap;

	if(!map){
		Error("Volume photon map is not initialized");
		exit(1);
	}
	
    float t0, t1;
	if (!vr || !vr->IntersectP(ray, &t0, &t1)) return 0.f;
	// Do multiple scattering volume integration in _vr_
	Spectrum Lv(0.);

	// Prepare for volume integration stepping
	int N = Ceil2Int((t1-t0) / stepSize);
	float step = (t1 - t0) / N; //split into N steps
	Spectrum Tr(1.f);
	Point p = ray(t0), pPrev;
	Vector w = -ray.d; 

	if (sample)
		t0 += sample->oneD[scatterSampleOffset][0] * step;
	else
		t0 += rng.RandomFloat() * step;
	// Compute sample patterns for multiple scattering samples
	float *samp = (float *)malloc(3 * N * sizeof(float));
	
    LatinHypercube(samp, N, 3, rng);
	
    int sampOffset = 0;

	for (int i = 0; i < N; ++i, t0 += step) {
		// Advance to sample at _t0_ and update _T_
		pPrev = p;
		p = ray(t0);
		Spectrum stepTau = vr->tau(Ray(pPrev, p - pPrev, 0, 1), .5f * stepSize, rng.RandomFloat());
		Tr *= Exp(-stepTau);

		// Possibly terminate raymarching if transmittance is small
		if (Tr.y() < 1e-3) {
			const float continueProb = .5f;
			if (rng.RandomFloat() > continueProb) break;
			Tr /= continueProb;
		}
		// Compute multiple-scattering source term at _p_
		
        float maxDistS = maxDistSquared;
		Lv += Tr * vr->Lve(p, w,1e-4f);
        
		Spectrum ss(1.f); //= vr->sigma_s(p, w);
		if (!ss.IsBlack()) {
			// Collect nearby photons
			VPhotonProcess proc(nLookup, p);
			proc.photons = (ClosePhoton *)alloca(nLookup * sizeof(ClosePhoton));
			map->Lookup(p, proc, maxDistS);

			if(maxDistS > 0.f){
				float scale = 0.75f / (float(nVolumePaths) * powf(maxDistS, 1.5f)  * M_PI);
				//float scale = 0.75f / (powf(maxDistS, 1.5f)  * M_PI * vr->sigma_t(p, w));
				ClosePhoton *photons = proc.photons;
				int nFound = proc.nfound;
				Spectrum L(0.);
				for (int i = 0; i < nFound; ++i){
					if(photons[i].photon->alpha.y() < 10000)
						L += vr->p(p, photons[i].photon->wi, w,1e-4f) * photons[i].photon->alpha;
				}
				Lv += Tr * scale * L / vr->sigma_t(p, w,1e-4f);
			}

		}
		sampOffset += 3;
	}

	Lv *= step;
    *transmittance = Tr;
	free((void*)samp);
	return Lv * step;
}

// PhotonIntegrator Method Definitions
VolumePhotonIntegrator::VolumePhotonIntegrator(int nvol, int nl, int mphodepth,
                                               float mdist, int gs, float ga, float mstep, float steps)
{
    nLookup = nl;
    maxPhotonDepth = mphodepth;
    maxDistSquared = mdist * mdist;
    cosGatherAngle = cos(Radians(ga));
    gatherSamples = gs;
    lightSampleOffsets = NULL;
    bsdfSampleOffsets = NULL;
    
    nVolumePhotonsWanted = nvol;
    marchStep = mstep;
    stepSize = steps;
}

VolumePhotonIntegrator *CreatePhotonMapVolumeIntegrator(const ParamSet &params) {
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 4);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
    float maxDist = params.FindOneFloat("maxdist", .1f);
    float gatherAngle = params.FindOneFloat("gatherangle", 10.f);
    
    //new parameters
    int nVolume = params.FindOneInt("volumephotons", 100000);
    float marchStep = params.FindOneFloat("marchstep", 4.0f);
    float stepSize  = params.FindOneFloat("stepsize", 1.f);

    return new VolumePhotonIntegrator(nVolume, nUsed, maxPhotonDepth, maxDist, gatherSamples,
                                      gatherAngle, marchStep, stepSize);
}


