#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

// This class represents an implementation of excercise 7.3 from the PBRT v3 textbook.
// 7.3: "Implement the improved multi-jittered sampling method introduced by Kensler (2013) as a new Sampler in pbrt."
//
// See https://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction/Exercises for details
// Paper link: http://eastfarthing.com/publications/cmj.pdf

#ifndef PBRT_SAMPLERS_MULTIJITTER_H
#define PBRT_SAMPLERS_MULTIJITTER_H

// samplers/multijitter.h*
#include "sampler.h"
#include "lowdiscrepancy.h"

namespace pbrt {

class MultiJitterSampler : public PixelSampler {
	public:
    MultiJitterSampler(
      int xPixelSamples, 
      int yPixelSamples, 
      bool jitterSamples, 
      int nSampledDimensions)
    : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions),
      xPixelSamples(xPixelSamples),
      yPixelSamples(yPixelSamples),
      jitterSamples(jitterSamples) {}
    
    void StartPixel(const Point2i &);
    std::unique_ptr<Sampler> Clone(int seed);

  private:
    const int xPixelSamples, yPixelSamples;
    const bool jitterSamples;

    // Finds reasonable n and m such that n*m >= count, for the purposes of creating strata.
    std::vector<int> DecomposeCount(int count);

    void CorrelatedMultiJitter(Point2f *samples, int m, int n, int count, RNG &rng);
};

MultiJitterSampler *CreateMultiJitterSampler(const ParamSet &params);

}
#endif 