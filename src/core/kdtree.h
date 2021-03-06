#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_KDTREE_H
#define PBRT_CORE_KDTREE_H

// core/kdtree.h*
#include "pbrt.h"
#include "geometry.h"

// KdTree Declarations
struct KdNode {

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
    
    // KdNode Data
    float splitPos;
    float ri, ro;
    uint32_t splitAxis:2; //2 bits
    uint32_t hasLeftChild:1, rightChild:29; //30 bits
};


template <typename NodeData> class KdTree {
public:
    // KdTree Public Methods
    KdTree(const vector<NodeData> &data);
    ~KdTree() {
        FreeAligned(nodes);
        FreeAligned(nodeData);
    }
    template <typename LookupProc> void Lookup(const Point &p,
            LookupProc &process, float &maxDistSquared) const;
    
    template <typename LookupProc> void buildRadii(uint32_t nodeNum,
            const Point &p, LookupProc &proc, float &maxDistSquared);
private:
    // KdTree Private Methods
    void recursiveBuild(uint32_t nodeNum, int start, int end,
        const NodeData **buildNodes);
    
    template <typename LookupProc> void privateLookup(uint32_t nodeNum,
        const Point &p, LookupProc &process, float &maxDistSquared) const;

    // KdTree Private Data
    KdNode *nodes;
    NodeData *nodeData;
    uint32_t nNodes, nextFreeNode;
};


template <typename NodeData> struct CompareNode {
    CompareNode(int a) { axis = a; }
    int axis;
    bool operator()(const NodeData *d1, const NodeData *d2) const {
        return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
                                            d1->p[axis] < d2->p[axis];
    }
};


struct APhoton {
    APhoton(const Point &pp, const Spectrum &wt, const Vector &w)
    : p(pp), alpha(wt), wi(w) { }
    APhoton() { }
    Point p;
    Spectrum alpha;
    Vector wi;
    float ri;
};

// KdTree Method Definitions
template <typename NodeData>
KdTree<NodeData>::KdTree(const vector<NodeData> &d) {

    nNodes = d.size();
    nextFreeNode = 1;
    
    nodes = AllocAligned<KdNode>(nNodes);
//    nodes = (KdNode*)alloca(nNodes * sizeof(KdNode));
    nodeData = AllocAligned<NodeData>(nNodes);
    vector<const NodeData *> buildNodes(nNodes, NULL);
    for (uint32_t i = 0; i < nNodes; ++i)
        buildNodes[i] = &d[i];
    // Begin the KdTree building process
    recursiveBuild(0, 0, nNodes, &buildNodes[0]);
}


template <typename NodeData> void
KdTree<NodeData>::recursiveBuild(uint32_t nodeNum, int start, int end,
        const NodeData **buildNodes) {
    // Create leaf node of kd-tree if we've reached the bottom
    if (start + 1 == end) {
        nodes[nodeNum].initLeaf();
        nodeData[nodeNum] = *buildNodes[start];
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
                     &buildNodes[end], CompareNode<NodeData>(splitAxis));

    // Allocate kd-tree node and continue recursively
    nodes[nodeNum].init(buildNodes[splitPos]->p[splitAxis], splitAxis);
    nodeData[nodeNum] = *buildNodes[splitPos];
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


template <typename NodeData> template <typename LookupProc>
void KdTree<NodeData>::Lookup(const Point &p, LookupProc &proc,
                              float &maxDistSquared) const {
    privateLookup(0, p, proc, maxDistSquared);
}


template <typename NodeData> template <typename LookupProc>
void KdTree<NodeData>::buildRadii(uint32_t nodeNum, const Point &p, LookupProc &proc,
                              float &maxDistSquared){
    KdNode *node = &nodes[nodeNum];
    
    //compute m at a given node, distance to the m closest photon
    privateLookup(nodeNum, nodeData[nodeNum].p, proc, maxDistSquared);
    
    //node->ri = sqrtf(maxDistSquared);
    
    if (node->hasLeftChild)
        buildRadii(nodeNum+1, p, proc, maxDistSquared);

    if (node->rightChild)
        buildRadii(node->rightChild, p, proc, maxDistSquared);
}


template <typename NodeData> template <typename LookupProc>
void KdTree<NodeData>::privateLookup(uint32_t nodeNum, const Point &p,
        LookupProc &process, float &maxDistSquared) const {

    KdNode *node = &nodes[nodeNum];
    
    // Process kd-tree node's children
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

    
    // Hand kd-tree node to processing function
    float dist2 = DistanceSquared(nodeData[nodeNum].p, p);
    if (dist2 < maxDistSquared)
        process(p, nodeData[nodeNum], dist2, maxDistSquared);
}



#endif // PBRT_CORE_KDTREE_H
