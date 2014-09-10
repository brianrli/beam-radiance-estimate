// volumes/exponential.cpp*
#include "stdafx.h"
#include "volumes/exponential.h"
#include "paramset.h"

// ExponentialDensity Method Definitions
ExponentialDensity *CreateExponentialVolumeRegion(const Transform &volume2world,
        const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Spectrum Le = params.FindOneSpectrum("Le", 0.);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
    float a = params.FindOneFloat("a", 1.);
    float b = params.FindOneFloat("b", 1.);
    Vector up = params.FindOneVector("updir", Vector(0,1,0));
    return new ExponentialDensity(sigma_a, sigma_s, g, Le, BBox(p0, p1),
        volume2world, a, b, up);
}


