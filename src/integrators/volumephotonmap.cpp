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
    Photon(const Point &pp, const Spectrum &wt, const Vector &w, float n)
    : p(pp), alpha(wt), wi(w)
    {
        ri = 0.f;
    }
    Photon(const Photon&ph){
        p = ph.p;
        alpha = ph.alpha;
        wi = ph.wi;
        ri = ph.ri;
        bound = ph.bound;
    }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
    //bbh+per photon radius
    float ri;
    BBox bound;
};


// VKdTree Declarations
struct VKdNode {
    
    void init(float p, uint32_t a) {
        splitPos = p;
        splitAxis = a;
        rightChild = (1<<29)-1;
        hasLeftChild = 0;
    }
    
    void initLeaf() {
        splitAxis = 3;
        rightChild = (1<<29)-1;
        hasLeftChild = 0;
    }
    
    // VKdNode Data
    float splitPos;
    uint32_t splitAxis:2; //2 bits
    uint32_t hasLeftChild:1, rightChild:29; //30 bits
};


class VKdTree {
    friend class VBVHAccel; //can access VKdTree private
public:
    // VKdTree Public Methods
    VKdTree(const vector<Photon> &data);
    ~VKdTree() {
        FreeAligned(nodes);
        FreeAligned(photon);
    }
    template <typename LookupProc> void Lookup(const Point &p,
                                               LookupProc &process, float &maxDistSquared) const;
    
    template <typename LookupProc> void buildRadii(uint32_t nodeNum,
                                                   const Point &p, LookupProc &proc, float &maxDistSquared);
private:
    // VKdTree Private Methods
    void recursiveBuild(uint32_t nodeNum, int start, int end,
                        const Photon **buildNodes);
    
    template <typename LookupProc> void privateLookup(uint32_t nodeNum,
                                                      const Point &p, LookupProc &process, float &maxDistSquared) const;
    
    template <typename LookupProc> void radiusLookup(uint32_t nodeNum,
                                                     const Point &p, LookupProc &process, float &maxDistSquared) const;
    
    // VKdTree Private Data
    VKdNode *nodes;
    Photon *photon;
    uint32_t nNodes, nextFreeNode;
};


struct VCompareNode {
    VCompareNode(int a) { axis = a; }
    int axis;
    bool operator()(const Photon *d1, const Photon *d2) const {
        return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
        d1->p[axis] < d2->p[axis];
    }
};


// VKdTree Method Definitions
VKdTree::VKdTree(const vector<Photon> &d) {
    
    nNodes = d.size();
    nextFreeNode = 1;
    
    size_t bob = sizeof(VKdNode);
    
    nodes = AllocAligned<VKdNode>(nNodes);
    photon = VAllocAligned<Photon>(nNodes);
    
    vector<const Photon *> buildNodes(nNodes, NULL);
    
    for (uint32_t i = 0; i < nNodes; ++i)
        buildNodes[i] = &d[i];
    
    // Begin the VKdTree building process
    recursiveBuild(0, 0, nNodes, &buildNodes[0]);
}


void VKdTree::recursiveBuild(uint32_t nodeNum, int start, int end,
                             const Photon **buildNodes) {
    // Create leaf node of VKd-tree if we've reached the bottom
    if (start + 1 == end) {
        nodes[nodeNum].initLeaf();
        photon[nodeNum] = *buildNodes[start];
        return;
    }
    
    // Choose split direction and partition data
    
    // Compute bounds of data from _start_ to _end_
    BBox bound;
    for (int i = start; i < end; ++i)
        bound = Union(bound, buildNodes[i]->p);
    int splitAxis = bound.MaximumExtent();
    int splitPos = (start+end)/2;
    std::nth_element(&buildNodes[start], &buildNodes[splitPos],
                     &buildNodes[end], VCompareNode(splitAxis));
    
    // Allocate VKd-tree node and continue recursively
    nodes[nodeNum].init(buildNodes[splitPos]->p[splitAxis], splitAxis);
    photon[nodeNum] = *buildNodes[splitPos];
    if (start < splitPos) {
        nodes[nodeNum].hasLeftChild = 1;
        uint32_t childNum = nextFreeNode++;
        recursiveBuild(childNum, start, splitPos, buildNodes);
    }
    if (splitPos+1 < end) {
        nodes[nodeNum].rightChild = nextFreeNode++;
        recursiveBuild(nodes[nodeNum].rightChild, splitPos+1,
                       end, buildNodes);
    }
}


template <typename LookupProc>
void VKdTree::Lookup(const Point &p, LookupProc &proc,
                     float &maxDistSquared) const {
    privateLookup(0, p, proc, maxDistSquared);
}

//In charge of building Radii
template <typename LookupProc>
void VKdTree::buildRadii(uint32_t nodeNum, const Point &p, LookupProc &proc,
                         float &maxDistSquared){
    VKdNode *node = &nodes[nodeNum];
    float dist2 = maxDistSquared;
    int axis = node->splitAxis;
    
    //compute m at a given node, distance to the m closest photon
    privateLookup(0, photon[nodeNum].p, proc, dist2);
    photon[nodeNum].ri = .2f * max(sqrtf(dist2),sqrtf(proc.get_radius()))
    * powf((float)proc.nLookup,0.333333333333333333f); //increase size of radius
    
    Photon ndata = photon[nodeNum];
    
    //no radiuses of zero
    Assert(ndata.ri != 0.f);
    proc.clear(); //perhaps not the most efficient
    
    if (node->hasLeftChild)
        buildRadii(nodeNum+1, p, proc, maxDistSquared);
    proc.clear();
    
    if (node->rightChild < nNodes)
        buildRadii(node->rightChild, p, proc, maxDistSquared);
}


template <typename LookupProc>
void VKdTree::privateLookup(uint32_t nodeNum, const Point &p,
                            LookupProc &process, float &maxDistSquared) const {
    VKdNode *node = &nodes[nodeNum];
    // Process VKd-tree node's children
    int axis = node->splitAxis;
    
    if (axis != 3) {
        
        //distance squared from point to photon
        float dist2 = (p[axis] - node->splitPos) * (p[axis] - node->splitPos);
        
        if (p[axis] <= node->splitPos) { //below split position
            if (node->hasLeftChild)
                privateLookup(nodeNum+1, p, process, maxDistSquared);
            
            if (dist2 < maxDistSquared && node->rightChild < nNodes)
                privateLookup(node->rightChild, p, process, maxDistSquared);
        }
        else {
            if (node->rightChild < nNodes)
                privateLookup(node->rightChild, p, process, maxDistSquared);
            
            if (dist2 < maxDistSquared && node->hasLeftChild)
                privateLookup(nodeNum+1, p, process, maxDistSquared);
        }
    }
    
    
    // Hand VKd-tree node to processing function
    float dist2 = DistanceSquared(photon[nodeNum].p, p);
    
    if (dist2 < maxDistSquared && (p!=photon[nodeNum].p)) //except same
        process(p, photon[nodeNum], dist2, maxDistSquared);
}


struct VBVHBuildNode {
    // BVHBuildNode Public Methods
    VBVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
    }
    void InitInterior(uint32_t axis, VBVHBuildNode *c0, VBVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;

        if(c1 != NULL)
            photon.bound = Union(c1->photon.bound,photon.bound);

        if(c0 != NULL)
            photon.bound = Union(c0->photon.bound,photon.bound);

        splitAxis = axis;
        nPrimitives = 1;
    }
    Photon photon;
    VBVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};

//rapid construction scheme
VBVHBuildNode* VBVHAccel::recursiveBuild(MemoryArena &buildArena, VKdTree &volumemap, uint32_t nodenum, uint32_t *totalNodes){

     (*totalNodes)++;
    VBVHBuildNode *node = buildArena.Alloc<VBVHBuildNode>();
    node->photon = volumemap.photon[nodenum];
    Photon swag = volumemap.photon[nodenum];
    Assert(node->photon.ri != 0.f);
//  Bounding box given point and radius
    node->photon.bound = BBox(node->photon.p);
    node->photon.bound.Expand(node->photon.ri*3.f);

    VKdNode *kdnode = &volumemap.nodes[nodenum];
    
    VBVHBuildNode *leftchild = NULL;
    VBVHBuildNode *rightchild = NULL;
    
//  if hasleftchild
    if(kdnode->hasLeftChild){
        leftchild = recursiveBuild(buildArena, volumemap, nodenum+1,totalNodes);
    }
    
//  if has rightchild
    if(kdnode->rightChild < volumemap.nNodes){
        rightchild = recursiveBuild(buildArena, volumemap, kdnode->rightChild,totalNodes);
    }
    
    node->InitInterior(kdnode->splitAxis,leftchild, rightchild);
    return node;
}

struct LinearVBVHNode {
    BBox bounds;
    uint32_t secondChildOffset; //second child
    uint8_t axis;         // interior node: xyz
    Photon photon;
    bool hasleft;
    bool hasright;
};

uint32_t VBVHAccel::flattenVBVHTree(VBVHBuildNode *node, uint32_t *offset) {
    //no child
    if (!node)
        return *offset;    

    //    int u[64];
//    int test = 0;
//    for(int i = 0; i < 64; i++){
//        u[i]=i;
//    }
//    int a = u[test++];
//    int b = u[test--];
    
    LinearVBVHNode *linearNode = &nodes[*offset];
    
    uint32_t myOffset = (*offset)++;
    
    linearNode->bounds = node->photon.bound;
    linearNode->axis = node->splitAxis;
    linearNode->photon = node->photon;
    linearNode->photon.bound = node->photon.bound;
    
    
    if(node->children[0] != NULL){
        flattenVBVHTree(node->children[0], offset);
        linearNode->hasleft = true;
    }
    else
        linearNode->hasleft = false;
    
    //do for bottom
    if(node->children[1] != NULL){
        linearNode->secondChildOffset = flattenVBVHTree(node->children[1],offset);
        linearNode->hasright = true;
    }
    else
        linearNode->hasright = false;
    
    return myOffset;
}

VBVHAccel::VBVHAccel(VKdTree &volumemap) {

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    VBVHBuildNode *root = recursiveBuild(buildArena, volumemap, 0, &totalNodes);
    // Compute representation of depth-first traversal of BVH tree
    nodes = AllocAligned<LinearVBVHNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearVBVHNode;
    uint32_t offset = 0;
    flattenVBVHTree(root, &offset);
    Assert(offset == totalNodes);
}

VBVHAccel::~VBVHAccel() {
    FreeAligned(nodes);
}


static inline bool IntersectP(const BBox &bounds, const Ray &ray,
                              const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    
    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}

struct VIntersect{
    VIntersect(Photon p, float d):
    dist(d)
    {
        photon = Photon(p);
    }
    
    Photon photon;
    float dist;
    int a;
};

bool VBVHAccel::Intersect(const Ray &ray, vector<VIntersect> &vintersect) const {
    
    if (!nodes) return false;
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    for (int i = 0; i<64; i++){
        todo[i]=0;
    }
    
    while (true) {
        const LinearVBVHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectP(node->photon.bound, ray, invDir, dirIsNeg))
        {
            float dist = Cross(ray.d, ray.o-node->photon.p).Length();
            // Intersect ray with primitives in leaf BVH node
            if (dist < node->photon.ri)
            {
                vintersect.push_back(VIntersect(node->photon,dist)); //Bound doesn't work
//                VIntersect sample = VIntersect(node->photon,dist);
                hit = true;
                
            }

            // Put far BVH node on _todo_ stack, advance to near node
            if (dirIsNeg[node->axis])
            {
                if(node->hasleft)
                    todo[todoOffset++] = nodeNum + 1; //place
                if (node->hasright)
                    nodeNum = node->secondChildOffset;
                else{
                    if (todoOffset == 0) break;
                    nodeNum = todo[--todoOffset];
                }
            }
            else
            {
                if (node->hasright)
                    todo[todoOffset++] = node->secondChildOffset;
                if(node->hasleft)
                    nodeNum = nodeNum + 1;
                else{
                    if (todoOffset == 0) break;
                    nodeNum = todo[--todoOffset];
                }
            }
            
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    return hit;
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

struct RadiancePhoton {
    RadiancePhoton(const Point &pp, const Normal &nn)
    : p(pp), n(nn), Lo(0.f) { }
    RadiancePhoton() { }
    Point p;
    Normal n;
    Spectrum Lo;
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
        
        if (nfound < nLookup) { //still looking
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
    void clear(){
        for (int i = 0; i < nLookup; i++){
            photons[i] = NULL;
        }
        nfound = 0;
    }
    float get_radius(){
        float distance = 0.f;
        for(int i = 0; i<nfound; i++){
            if(photons[i].distanceSquared > distance){
                distance = photons[i].distanceSquared;
            }
        }
        return distance;
    }
	const Point &p;
	ClosePhoton *photons;
	u_int nLookup;
	mutable u_int nfound;
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

//float SampleScattering(const Vector &wi, float u1, float u2, Vector &wo){
//    wo = VSampleHG(wi, 0.5f, u1, u2);
//    return VPhaseHG(wi, wo, 0.5f);
//}

float SampleScattering(const Vector &wi, float u1, float u2, Vector &wo){
    float phi = u1*2*M_PI;
    float theta = u2*M_PI;
    float cosTheta = cos(theta);
    wo.x = cosTheta*cos(phi);
    wo.y = cosTheta*sin(phi);
    wo.z = sin(theta);
    return INV_TWOPI/2.f;
}

void VolumePhotonShootingTask::Run() {
    
    // Declare local variables for _VolumePhotonShootingTask_
    VolumeRegion *volume = scene->volumeRegion;
    float marchstep = integrator->marchStep;
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localVolumePhotons; //new data structure
    uint32_t totalPaths = 0;
    PermutedHalton halton(6, rng);
    //set finish conditions
    bool volumeDone = (integrator->nVolumePhotonsWanted == 0);

    while (true) {
//        const uint32_t blockSize = 4096;
        const uint32_t blockSize = 500;
        for (uint32_t i = 0; i < blockSize; ++i) {
            
            float u[6];
            halton.Sample(++totalPaths, u); //random sample
            
            // Choose light to shoot photon from
            float lightPdf; //prob distribution function
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum]; //choose light
            
            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;

            //Alpha (PBRT)
            VisibilityTester vis;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf); //initial

            if (!alpha.IsBlack()) {

                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);

                Intersection photonIsect;
                int nIntersections = 0;
                float vt0, vt1; //ray volume bounds
                Spectrum oldTransmittance(1.f);
                
                //if ray intersects with volume
                if (volume->IntersectP(photonRay, &vt0, &vt1)) {

                    ++nIntersections;
                    float inVolume = true;
                    bool firstInteraction = true;

                    // Enter volume photon mode till we leave the volumeRegion
                    float u2, u6;
                    
                    Ray volumeRay;
                    volumeRay.o = photonRay(vt0);
                    volumeRay.d = Normalize(photonRay.d); // Hard to know if it's normalized. Just do it again.
                    volumeRay.mint = 0;
                    
                    u2 = (rng.RandomFloat() + 0.5f) * marchstep; //[0,1) + 0.5 * marchstep
                    
                    if( u2 >= (vt1-vt0)/photonRay.d.Length()){ //if randomly generated step is larger than
//                        photonRay.o = volumeRay(vt1+0.00001f);
//                        alpha *= Exp(-volume->tau(volumeRay));
                        inVolume = false;
                    }
                    else{
                        volumeRay.maxt = u2; //create a new ray w/in  volume bounds
                    }                    

                    
                    while(inVolume)
                    {
                    // Ray march and decide to scattering/absorbing or move on
                        Spectrum transmittance = Exp(-volume->tau(volumeRay)) * oldTransmittance;
                        alpha *= transmittance;

                        if(1){ // Enforce interaction
                            
                            Point interactPoint = volumeRay(volumeRay.maxt);
                                                        // Store a photon
                            if(!volumeDone && !firstInteraction){
                                Photon photon(interactPoint, alpha, -volumeRay.d,integrator->blur);
                                localVolumePhotons.push_back(photon);
                                progress.Update();
                            }
                            
                            // Choose scattering or absorbing
//                            float t = volume->sigma_t(interactPoint, volumeRay.d,volumeRay.time).y();
//                            float s = volume->sigma_s(interactPoint, volumeRay.d,volumeRay.time).y();
//                            float r = rng.RandomFloat();
//                            float m = r*t;
                            
                            if((rng.RandomFloat() * volume->sigma_t(interactPoint, volumeRay.d,volumeRay.time).y()) <
                               volume->sigma_s(interactPoint, volumeRay.d,volumeRay.time).y()){

                                // Scattering, sampling a new direction
                                oldTransmittance = Spectrum(1.f);
                                Vector newDirection;
                                float spdf = SampleScattering(volumeRay.d, rng.RandomFloat(), rng.RandomFloat(), newDirection);
                                alpha *= spdf;

                                // Specify new volumeRay
                                volumeRay.o = interactPoint;
//                                volumeRay.d = newDirection; //try with same direction
                                u6 = rng.RandomFloat() + 0.5f;
                                volumeRay.maxt = marchstep * u6; //new photonray

                                firstInteraction = false; //only use scattered
                                
                                if(!volume->Inside(volumeRay(volumeRay.maxt))){
                                    break;
                                }
                                if(scene->IntersectP(volumeRay))
                                    break;
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
    
    if (volumePhotons.size() > 0){
        volumeMap = new VKdTree(volumePhotons);
    
        //traverse volumeMap, assign valid radius to Photons
        float maxdist2 = maxDistSquared * 5.f;
        VPhotonProcess proc(5,Point()); //find m closest
        proc.photons = (ClosePhoton *)alloca(5 * sizeof(ClosePhoton));
        volumeMap->buildRadii(0, Point(), proc, maxdist2);
        
        BBH = new VBVHAccel(*volumeMap);
    }
}



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
	VKdTree *map = volumeMap;
    
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

float bi_kernel(float x){
    return 3*powf(M_PI,-1)*powf(1-powf(x,2),2.f);
}

//
//Spectrum VolumePhotonIntegrator::Li(const Scene *scene,
//                                    const Renderer *renderer,
//                                    const RayDifferential &ray,
//                                    const Sample *sample,
//                                    RNG &rng,
//                                    Spectrum *transmittance,
//                                    MemoryArena &arena) const {
//    
//	// Pointers assignment
//	if(!BBH){
//		Error("Volume photon map is not initialized");
//		exit(1);
//	}
//    VolumeRegion *vr = scene->volumeRegion;
//    
//    Vector w = -ray.d;
//	vector<VIntersect> intersections;
//    intersections.reserve(100);
//    BBH->Intersect(ray, intersections);
//    Spectrum Tr(1.f);
//    
//    float nints = intersections.size();
//    Spectrum L(1.f); //accumulated radiance
//    
//    for ( auto i = intersections.begin(); i != intersections.end(); i++ ) {
//        Point p = i->photon.p;
//        float di = i->dist;
//        float ri = i->photon.ri;
//        Vector wi = i->photon.wi;
//        float ti = Dot(p-ray.o,w);
//        const Point xp = ray.o + (ti*w);
//        const Point o = ray.o;
//        
//        //Ki(x,w,s,xi,ri)
//        float Ki = powf(ri,-2.f)*bi_kernel(di/ri);
//        
//        //Transmittance
//        Spectrum stepTau = vr->tau(Ray(o, xp-o, 0, 1));
//        Tr = Exp(-stepTau);
//        
//        //Scatter
//        Spectrum ss = vr->sigma_s(xp, w,ray.time);
//        
//        //Phase Function
//        Spectrum pf = vr->p(p, wi, w,1e-4f);
//        
//        //Power
//        Spectrum weight = i->photon.alpha;
//        
//        L += Tr*Ki*ss*pf*weight;
//    }
//    if(nints>0.f)
//        L /= nints;
//    *transmittance = Tr;
//	return L;
//}


// PhotonIntegrator Method Definitions
VolumePhotonIntegrator::VolumePhotonIntegrator(int nvol, int nl, int mphodepth,
                                               float mdist, int gs, float ga, float mstep, float steps,
                                               float blr)
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
    blur = blr;
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
    float marchStep = params.FindOneFloat("marchstep", 2.0f);
    float stepSize  = params.FindOneFloat("stepsize", 4.f);
    float blur = params.FindOneFloat("blur", 4.f);

    return new VolumePhotonIntegrator(nVolume, nUsed, maxPhotonDepth, maxDist, gatherSamples,
                                      gatherAngle, marchStep, stepSize,blur);
}


