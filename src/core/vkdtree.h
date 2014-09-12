#if defined(_MSC_VER)
#pragma once
#endif

#define PBRT_CORE_VKdTREE_H

// core/VKdtree.h*
#include "pbrt.h"
#include "geometry.h"
#include "integrators/volumephotonmap.h"

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


template <typename NodeData> class VKdTree {
public:
    // VKdTree Public Methods
    VKdTree(const vector<NodeData> &data);
    ~VKdTree() {
        FreeAligned(nodes); 
        FreeAligned(nodeData);
    }
    template <typename LookupProc> void Lookup(const Point &p,
                                               LookupProc &process, float &maxDistSquared) const;
    
    template <typename LookupProc> void buildRadii(uint32_t nodeNum,
                                                   const Point &p, LookupProc &proc, float &maxDistSquared);
private:
    // VKdTree Private Methods
    void recursiveBuild(uint32_t nodeNum, int start, int end,
                        const NodeData **buildNodes);
    
    template <typename LookupProc> void privateLookup(uint32_t nodeNum,
                                                      const Point &p, LookupProc &process, float &maxDistSquared) const;
    
    template <typename LookupProc> void radiusLookup(uint32_t nodeNum,
                                                      const Point &p, LookupProc &process, float &maxDistSquared) const;
    
    // VKdTree Private Data
    VKdNode *nodes;
    NodeData *nodeData;
    uint32_t nNodes, nextFreeNode;
};


template <typename NodeData> struct VCompareNode {
    VCompareNode(int a) { axis = a; }
    int axis;
    bool operator()(const NodeData *d1, const NodeData *d2) const {
        return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
        d1->p[axis] < d2->p[axis];
    }
};


// VKdTree Method Definitions
template <typename NodeData>
VKdTree<NodeData>::VKdTree(const vector<NodeData> &d) {
    
    nNodes = d.size();
    nextFreeNode = 1;

    nodes = AllocAligned<VKdNode>(nNodes);
    nodeData = AllocAligned<NodeData>(nNodes);
    vector<const NodeData *> buildNodes(nNodes, NULL);
    
    for (uint32_t i = 0; i < nNodes; ++i)
        buildNodes[i] = &d[i];
    
    // Begin the VKdTree building process
    recursiveBuild(0, 0, nNodes, &buildNodes[0]);
}


template <typename NodeData> void
VKdTree<NodeData>::recursiveBuild(uint32_t nodeNum, int start, int end,
                                 const NodeData **buildNodes) {
    // Create leaf node of VKd-tree if we've reached the bottom
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
                     &buildNodes[end], VCompareNode<NodeData>(splitAxis));
    
    // Allocate VKd-tree node and continue recursively
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
void VKdTree<NodeData>::Lookup(const Point &p, LookupProc &proc,
                              float &maxDistSquared) const {
    privateLookup(0, p, proc, maxDistSquared);
}

//In charge of building Radii
template <typename NodeData> template <typename LookupProc>
void VKdTree<NodeData>::buildRadii(uint32_t nodeNum, const Point &p, LookupProc &proc,
                                  float &maxDistSquared){
    VKdNode *node = &nodes[nodeNum];
    u_int sweg = proc.nfound;
    float dist2 = maxDistSquared;
    
    int axis = node->splitAxis;
    
    if (axis != 3) {
        //compute m at a given node, distance to the m closest photon
        privateLookup(0, nodeData[nodeNum].p, proc, dist2);
        
        if(dist2==0.f){
            for (int i = 0; i < proc.nfound; i++){
                float a = proc.photons[i].distanceSquared;
                int b = 3;
            }
        }
        
        nodeData[nodeNum].ri = max(sqrtf(dist2),sqrtf(proc.get_radius()))
            * powf((float)proc.nLookup,0.16666666666f);
        
        NodeData ndata = nodeData[nodeNum];
        
        proc.clear(); //perhaps not the most efficient
        
        if (node->hasLeftChild)
            buildRadii(nodeNum+1, p, proc, maxDistSquared);

        proc.clear();
        
        if (node->rightChild < nNodes)
            buildRadii(node->rightChild, p, proc, maxDistSquared);
    }
}


template <typename NodeData> template <typename LookupProc>
void VKdTree<NodeData>::privateLookup(uint32_t nodeNum, const Point &p,
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

    float dist2 = DistanceSquared(nodeData[nodeNum].p, p);
    if(dist2==0){
        int a = 3;
    }

    if (dist2 < maxDistSquared)
        process(p, nodeData[nodeNum], dist2, maxDistSquared);
}

