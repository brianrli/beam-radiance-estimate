#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

// accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
struct BVHBuildNode;

// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
struct LinearBVHNode;

// BVHAccel Declarations
class BVHAccel : public Aggregate {
public:
    // BVHAccel Public Methods
    BVHAccel(const vector<Reference<Primitive> > &p, uint32_t maxPrims = 1,
             const string &sm = "sah");
    BBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~BVHAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
private:
    // BVHAccel Private Methods
    BVHBuildNode *recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
        uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims);
    uint32_t flattenBVHTree(BVHBuildNode *node, uint32_t *offset);

    // BVHAccel Private Data
    uint32_t maxPrimsInNode;
    enum SplitMethod { SPLIT_MIDDLE, SPLIT_EQUAL_COUNTS, SPLIT_SAH };
    SplitMethod splitMethod;
    vector<Reference<Primitive> > primitives;
    LinearBVHNode *nodes;
};


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_BVH_H
