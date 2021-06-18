#include "samplers/multijitter.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"
#include <math.h>

namespace pbrt {

  void MultiJitterSampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);

    // Standard shuffling, generating jittered 1D arrays then shuffling them.
    for (size_t i = 0; i < samples1D.size(); ++i) {
      StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng,
        jitterSamples);
      Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
      for (int64_t j = 0; j < samplesPerPixel; ++j) {
        int count = samples1DArraySizes[i];
        StratifiedSample1D(&sampleArray1D[i][j * count], count, rng,
          jitterSamples);
        Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
      }

    // Correlated multi jitter sampling
    for (size_t i = 0; i < samples2D.size(); ++i) {
      CorrelatedMultiJitter(
        &samples2D[i][0], 
        xPixelSamples, 
        yPixelSamples, 
        xPixelSamples * yPixelSamples, 
        rng);
    }
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
      for (int64_t j = 0; j < samplesPerPixel; ++j) {
        int count = samples2DArraySizes[i];

        std::vector<int> sizes = DecomposeCount(count);

        CorrelatedMultiJitter(
          &sampleArray2D[i][j * count], 
          sizes[0], 
          sizes[1],
          count,
          rng);
      }
    PixelSampler::StartPixel(p);
  }

  std::vector<int> MultiJitterSampler::DecomposeCount(int count) {
    std::vector<int> result(2);

    int m = floor(sqrt(count));
    if (m * m >= count) {
      result[0] = m;
      result[1] = m;
      return result;
    }
    else if ((m + 1) * m >= count) {
      result[0] = m + 1;
      result[1] = m;
      return result;
    }
    else {
      result[0] = m + 1;
      result[1] = m + 1;
      return result;
    }
  }

  void MultiJitterSampler::CorrelatedMultiJitter(
    Point2f *samples, int m, int n, int count, RNG &rng) 
  {
    std::vector<Point2f> shufflePoints = std::vector<Point2f>(n * m);

    // produce canonical arrangement
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < m; ++i) {
        shufflePoints[j * m + i].x = (i + (j + rng.UniformFloat()) / n) / m;

        // stretch coordinates along y axis 
        Float stetchedVal = ((j + (i + rng.UniformFloat()) / m) / n) * (m * n / count);
        shufflePoints[j * m + i].y = stetchedVal;
      }
    }

    // shuffle x coordinates in each column
    for (int j = 0; j < n; ++j) {
      int other_row = j + rng.UniformUInt32(n - j);
      for (int i = 0; i < m; ++i) {
        std::swap(shufflePoints[j * m + i].x, shufflePoints[other_row * m + i].x);
      }
    }

    //shuffle y coordinates in each row
    for (int i = 0; i < m; ++i) {
      int other_col = i + rng.UniformUInt32(m - i);
      for (int j = 0; j < n; ++j) {
        std::swap(shufflePoints[j * m + i].y, shufflePoints[j * m + other_col].y);
      }
    }

    // clip last mn - N samples, outside unit square
    int sampleIndex = 0;
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        Point2f point = shufflePoints[j * m + i];

        if (point.y <= 1) {
          samples[sampleIndex].x = std::min(point.x, OneMinusEpsilon);
          samples[sampleIndex].y = std::min(point.y, OneMinusEpsilon);
          sampleIndex++;
        }
      }
    }
  }

  std::unique_ptr<Sampler> MultiJitterSampler::Clone(int seed) {
    MultiJitterSampler *mjs = new MultiJitterSampler(*this);
    mjs->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(mjs);
  }

  MultiJitterSampler *CreateMultiJitterSampler(const ParamSet &params) {
    bool jitter = params.FindOneBool("jitter", true);
    int xsamp = params.FindOneInt("xsamples", 4);
    int ysamp = params.FindOneInt("ysamples", 4);
    int sd = params.FindOneInt("dimensions", 4);
    if (PbrtOptions.quickRender) xsamp = ysamp = 1;
    return new MultiJitterSampler(xsamp, ysamp, jitter, sd);
  }
}