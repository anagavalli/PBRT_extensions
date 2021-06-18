#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

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

    std::vector<int> DecomposeCount(int count);

    void CorrelatedMultiJitter(Point2f *samples, int m, int n, int count, RNG &rng);
};

MultiJitterSampler *CreateMultiJitterSampler(const ParamSet &params);

}
#endif 