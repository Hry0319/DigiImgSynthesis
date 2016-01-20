
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_FILTERS_NLMEAN_H
#define PBRT_FILTERS_NLMEAN_H

// filters/NLmean.h*
#include "filter.h"

// NL-Mean Filter Declarations
class NLMeanFilter : public Filter {
public:
    // Non-local Mean Filter Public Methods
    NLMeanFilter(float r, float f, float k, float a)
        : Filter(r, r), alpha(a), f(f) , k(k), r(r){

	}

    float Evaluate(float x, float y) const;

	void NLFiltering(float * rgb, int xPixelCount, int yPixelCount);

	
	void NLMeanFilter::InitSummedAreaRGBTable(float *summedArea, unsigned int width, unsigned int height, unsigned int rgb)const;

private:
    // Non-local Mean Filter Private Data
    const float alpha;
    const float f;
	const float k;
	const float r;
	//const float epslon = 1e-10;

};


NLMeanFilter *CreateNLMeanFilter(const ParamSet &ps);

#endif // PBRT_FILTERS_NLMean_H
