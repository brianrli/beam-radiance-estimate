#if defined(_MSC_VER)
#pragma once
#endif

#define PBRT_CORE_VKdTREE_H

// core/VKdtree.h*
#include "pbrt.h"
#include "geometry.h"


struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w, float n)
    : p(pp), alpha(wt), wi(w)
    {
        ri = 0.f;
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


template <typename Photon> struct VCompareNode {
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

    nodes = VAllocAligned<VKdNode>(nNodes);
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
                     &buildNodes[end], VCompareNode<Photon>(splitAxis));
    
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
    
    if (axis != 3) {
        //compute m at a given node, distance to the m closest photon
        privateLookup(0, photon[nodeNum].p, proc, dist2);
        photon[nodeNum].ri = max(sqrtf(dist2),sqrtf(proc.get_radius()))
            * powf((float)proc.nLookup,0.333333333333333333f);
        
        Photon ndata = photon[nodeNum];
        
        proc.clear(); //perhaps not the most efficient
        
        if (node->hasLeftChild)
            buildRadii(nodeNum+1, p, proc, maxDistSquared);

        proc.clear();
        
        if (node->rightChild < nNodes)
            buildRadii(node->rightChild, p, proc, maxDistSquared);
    }
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

