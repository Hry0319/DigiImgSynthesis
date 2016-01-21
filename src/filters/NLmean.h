
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

#define EPSLON 1e-10

struct Pixel {
    Pixel() {
        _weightSum = 0.f;
		_nSamplesBox = 0.f;
		for (int i = 0; i < 3; ++i) 
		{
			_Lrgb[i] = 0.f;
			_LrgbSumBox[i] = 0.f;
			_LrgbSumSqrBox[i] = 0.f;
			_varSumBox[i] = 0.f;
			_varSumSqrBox[i] = 0.f;
		}
    }
    float Lxyz[3];
    float weightSum;
    float splatXYZ[3];
    float pad;
	// These store the sample values, as well as the cumulative weights, as
	// defined by the filter.
	float _Lrgb[3];
	float _weightSum;
	// The 'box' data is used to compute the variance of samples falling within
	// the boundary of a pixel. It is not affected by the reconstruction filter.
	int _nSamplesBox;
	float _LrgbSumBox[3];
	float _LrgbSumSqrBox[3];
	// Buffer variance
	float _tmp[3];
	float _varSumBox[3];
	float _varSumSqrBox[3];
};
//typedef BlockedArray<Pixel> NLPixel;

// NL-Mean Filter Declarations
class NLMeanFilter : public Filter {
public:

    // Non-local Mean Filter Public Methods
    NLMeanFilter(float r, float f, float k, float a)
        : Filter(r, r), alpha(a), f(f) , k(k), r(r){
			epslon = EPSLON;
	}


    float Evaluate(float x, float y) const;

	float *NLFiltering(BlockedArray<Pixel> *pixelsA, BlockedArray<Pixel> *pixelsB, int xPixelCount, int yPixelCount);

	
	void NLMeanFilter::InitSummedAreaRGBTable(float *summedArea, unsigned int width, unsigned int height, unsigned int rgb)const;

	float SummedAreaValue(int x,int y,int width,int height, float *RGBSumareaTable, int rgb)const
        {
			int x0 = x - r - 1;
			int y0 = y - r - 1;
			int x1 = x + r;
			int y1 = y + r;

			if(x0 < 0) x0 = 0;
			if(y0 < 0) y0 = 0;
			if(x1 > width-1) x1 = width-1;
			if(y1 > height-1) y1 = height-1;

			float upper = RGBSumareaTable[y0 * width * 3 + x1 *3 + rgb];
			float left  = RGBSumareaTable[y1 * width * 3 + x0 *3 + rgb];
			float upper_left = RGBSumareaTable[y0 * width * 3 + x0 *3 + rgb];
			float mine = RGBSumareaTable[y1 * width * 3 + x1 *3 + rgb];
			
			return mine - upper - left + upper_left ;
        }

private:
    // Non-local Mean Filter Private Data
    const float alpha;
    const float f;
	const float k;
	const float r;
	float epslon/* = 1e-10*/;

	float *RGBSumareaTableA;
	float *RGBSumareaTableB;

};


NLMeanFilter *CreateNLMeanFilter(const ParamSet &ps);

#endif // PBRT_FILTERS_NLMean_H
