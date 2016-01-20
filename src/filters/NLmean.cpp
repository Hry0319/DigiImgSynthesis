
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

float* NLMeanFilter::NLFiltering(float * rgb, int xPixelCount, int yPixelCount){
	
	int nPix = xPixelCount * yPixelCount;

	float *outRGB = new float[3*nPix];
   
	//int offset;
	//int xMin, xMax, yMin, yMax;

	
	RGBSumareaTable = new float[nPix*3];

	for (int index = 0; index < nPix; index++)
	{
		RGBSumareaTable[index*3]		= rgb[index*3];
		RGBSumareaTable[index*3 + 1]	= rgb[index*3 + 1];
		RGBSumareaTable[index*3 + 2]	= rgb[index*3 + 2];
	}

	//-------- sumareatable ------------------------
	InitSummedAreaRGBTable(RGBSumareaTable , xPixelCount, yPixelCount, 0);
	InitSummedAreaRGBTable(RGBSumareaTable , xPixelCount, yPixelCount, 1);
	InitSummedAreaRGBTable(RGBSumareaTable , xPixelCount, yPixelCount, 2);
	//----------------------------------------------

	//float *dis = new float[nPix*3];
	float *mean = new float[nPix*3];
	float *var  = new float[nPix*3];

	float patchArea = (2*f+1)*(2*f+1);
	

	//----------- calculate Var[p] --------------------//
	for (int y = 0; y < yPixelCount ; y++)
	{
		for (int x = 0; x < xPixelCount; x++)
		{
			int index = x + y*xPixelCount;
			
			mean[index*3     ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, 0) / patchArea;
			mean[index*3 + 1 ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, 1) / patchArea;
			mean[index*3 + 2 ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, 2) / patchArea;

			var[index*3 ] = 0;
			var[index*3 + 1] = 0;
			var[index*3 + 2] = 0;

			
			for (int py = -r; py <= r ; py++)
			{
				for (int px = -r; px <= r ; px++)
				{ 
					if( (x+px) >= 0 && (y+py) >= 0 && (x+px) < xPixelCount && (y+py) < yPixelCount)
					{
						var[index*3     ] += rgb[(x+px + (y+py)*xPixelCount)*3    ] - mean[index*3];
						var[index*3 + 1 ] += rgb[(x+px + (y+py)*xPixelCount)*3 + 1] - mean[index*3+1];
						var[index*3 + 2 ] += rgb[(x+px + (y+py)*xPixelCount)*3 + 2] - mean[index*3+2];
					}					
				}
			}
			var[index*3     ] /= patchArea;
			var[index*3 + 1 ] /= patchArea;
			var[index*3 + 2 ] /= patchArea;
			
			
		}
	}

	float weight, totalW;
	float dis;
	
	int qx, qy;

	for (int y = 0; y < yPixelCount ; y++)
	{
		for (int x = 0; x < xPixelCount; x++)
		{
			int index = x + y*xPixelCount;
			
			// filter at p
			totalW = 0;
			//outRGB = {};
			
			outRGB[index*3    ] = 0;
			outRGB[index*3 + 1] = 0;
			outRGB[index*3 + 2] = 0;

			for (int fy = -f; fy <= f ; fy++)
			{
				for (int fx = -f; fx <= f ; fx++)
				{
					qx = x + fx;
					qy = y + fy;
					
					if( (qx) >= 0 && (qy) >= 0 && (qx) < xPixelCount && (qy) < yPixelCount)
					{
						dis = 0;
						//for each patch q
						for (int py = -r; py <= r ; py++)
						{
							for (int px = -r; px <= r ; px++)
							{
								for(int i = 0; i < 3 ; i++)
								{				
									if( ( (x+px) >= 0 && (y+py) >= 0 && (x+px) < xPixelCount && (y+py) < yPixelCount) &&
										( (qx+px) >= 0 && (qy+py) >= 0 && (qx+px) < xPixelCount && (qy+py) < yPixelCount)  
									){
										int indexP = (x+px + (y+py)*xPixelCount)*3 + i;
										int indexQ = (qx+px + (qy+py)*xPixelCount)*3 + i;

										int m = min(var[indexP], var[indexQ]);
							
										dis += (pow(rgb[indexP]-rgb[indexQ] , 2) 
													- alpha*(var[indexP] - m) ) / (epslon+ k*k*(var[indexP] + var[indexQ]));													
									}
								}
							}
						}
					
						dis /= (3*(2*f+1)*(2*f+1));
						weight = exp(-1* max(0.f, dis));
						totalW += weight;

						outRGB[index*3 ] += rgb[(qx + qy*xPixelCount)*3 ] * weight;
						outRGB[index*3 + 1] += rgb[(qx + qy*xPixelCount)*3 + 1] * weight;
						outRGB[index*3 + 2] += rgb[(qx + qy*xPixelCount)*3 + 2] * weight;
					}
				}
			}
			outRGB[index*3 ] /= totalW ; 
			outRGB[index*3 + 1] /= totalW ; 
			outRGB[index*3 + 2] /= totalW ; 

		}
	}

	return outRGB ;
}

void NLMeanFilter::InitSummedAreaRGBTable(float *summedArea, unsigned int width, unsigned int height, unsigned int rgb)const
{
    for(unsigned int i = 0; i < height; i++ )
    {
        for(unsigned int j = 0; j < width; j++)
        {
			int index = (i*width + j)*3 + rgb;

			if(i==0 && j==0){
				summedArea[index] = summedArea[index];
			}else if(i==0 && j!=0){
				summedArea[index] = summedArea[index-3] + summedArea[index];
			}else if(i!=0 && j==0){
				summedArea[index] = summedArea[index-width*3] + summedArea[index];
			}else{
				summedArea[index] = summedArea[index-3] + summedArea[index-width*3] - summedArea[index-3-width*3] + summedArea[index];
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


