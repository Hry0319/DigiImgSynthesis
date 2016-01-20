
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


// filters/gaussian.cpp*
#include "stdafx.h"
#include "filters/NLMean.h"
#include "paramset.h"
#include "math.h"

// NLMean Filter Method Definitions
float NLMeanFilter::Evaluate(float x, float y) const {
    return 1.f;
}

void NLMeanFilter::NLFiltering(float * rgb, int xPixelCount, int yPixelCount){
	
	//int offset;
	int xMin, xMax, yMin, yMax;

	int nPix = xPixelCount * yPixelCount;

	float *RGBSumareaTable = new float[nPix*3];

	for (int index = 0; index < nPix; index++)
	{
		RGBSumareaTable[index*3]		= rgb[index*3];
		RGBSumareaTable[index*3 + 1]	= rgb[index*3 + 1];
		RGBSumareaTable[index*3 + 2]	= rgb[index*3 + 2];
	}

	//-------- sumareatable ------------------------
	InitSummedAreaRGBTable(RGBSumareaTable , xPixelCount*3, yPixelCount, 0);
	InitSummedAreaRGBTable(RGBSumareaTable , xPixelCount*3, yPixelCount, 1);
	InitSummedAreaRGBTable(RGBSumareaTable , xPixelCount*3, yPixelCount, 2);
	//----------------------------------------------

	float dis[3]  = {0.0f, 0.0f, 0.0f};
	float mean[3] = {0.0f, 0.0f, 0.0f};

	for (int y = 0; y < yPixelCount ; y++)
	{
		for (int x = 0; x < xPixelCount; x++)
		{
			for (int i = 0; i < 2*r+1 ; i++)
			{
				for (int j = 0; j < 2*r+1; j++)
				{

				}
			}
		}
	}

	for(int y = 0; y < yPixelCount ; y++){
		for(int x = 0; x < xPixelCount; x++){
			xMin = max(0, (int) (x-r));
			xMax = min(0, (int) (x+r));
			yMin = max(0, (int) (y-r));
			yMax = min(0, (int) (y+r));

			

		}
	}
}

void NLMeanFilter::InitSummedAreaRGBTable(float *summedArea, unsigned int width, unsigned int height, unsigned int rgb)const
{
    for(unsigned int i = 0; i < height; i++ )
    {
        for(unsigned int j = 0; j < width; j++)
        {
            int index = i*width + j*3 + rgb;
            if(i==0)
            {
                if(j!=0)
                    summedArea[index] += summedArea[index - 3];
            }
            else if(j==0)
            {
                if(i!=0)
                    summedArea[index] += summedArea[index - width];
            }
            else
            {
                summedArea[index] = summedArea[index] + summedArea[index-3] + summedArea[index-width] - summedArea[index-width-3];
            }
        }
    }
}

NLMeanFilter *CreateNLMeanFilter(const ParamSet &ps) {
    // Find common filter parameters
    float r = ps.FindOneFloat("rwidth", 7.f);
    float f = ps.FindOneFloat("fwidth", 3.f);
    float k = ps.FindOneFloat("kvalue", 0.45f);
	float alpha = ps.FindOneFloat("alpha", 0.5f);
    return new NLMeanFilter(r, f, k, alpha);
}


