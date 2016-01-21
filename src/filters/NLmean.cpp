
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
#include "film/dualfilm.h"
#include "paramset.h"
#include "math.h"
#include "filters/gaussian.h"


// NLMean Filter Method Definitions
float NLMeanFilter::Evaluate(float x, float y) const {
    return 1.f;
}

float* NLMeanFilter::NLFiltering(BlockedArray<Pixel> *pixelsA, BlockedArray<Pixel> *pixelsB, int xPixelCount, int yPixelCount){
	
	int nPix = xPixelCount * yPixelCount;
	this->nPixs = nPix;
	this->_xPixelCount = xPixelCount;
	this->_yPixelCount = yPixelCount;

	float *rgbA , *rgbB;
	float *out = new float[nPix];

	rgbA = Cal(pixelsA, xPixelCount, yPixelCount);
	rgbB = Cal(pixelsB, xPixelCount, yPixelCount);

	ImgVar_A = new float[nPix*3];
	ImgVar_B = new float[nPix*3];

	 // Compute observed variance
    for (int i = 0; i < nPix*3; i++) {
		float vv = pow((rgbA[i] - rgbB[i]), 2);
        ImgVar_A[i] = 2.f * vv / (1e-3f + pow(rgbA[i], 2));
        ImgVar_B[i] = 2.f * vv / (1e-3f + pow(rgbB[i], 2));
    }
    
    // Update the pixel costs
	UpdateError(ImgVar_A, pixelsA, ImgErr_A);
	UpdateError(ImgVar_B, pixelsB, ImgErr_B);

	for(int j =0; j < yPixelCount ; j++){
		for(int i =0; i < xPixelCount; i++){
			int index = i + j*xPixelCount;
			out[index*3   ] = (rgbA[index*3   ] * (*pixelsA)(i, j)._nSamplesBox + rgbB[index*3   ] * (*pixelsB)(i, j)._nSamplesBox) / ((*pixelsA)(i, j)._nSamplesBox + (*pixelsB)(i, j)._nSamplesBox);
			out[index*3 +1] = (rgbA[index*3 +1] * (*pixelsA)(i, j)._nSamplesBox + rgbB[index*3 +1] * (*pixelsB)(i, j)._nSamplesBox) / ((*pixelsA)(i, j)._nSamplesBox + (*pixelsB)(i, j)._nSamplesBox);
			out[index*3 +2] = (rgbA[index*3 +2] * (*pixelsA)(i, j)._nSamplesBox + rgbB[index*3 +2] * (*pixelsB)(i, j)._nSamplesBox) / ((*pixelsA)(i, j)._nSamplesBox + (*pixelsB)(i, j)._nSamplesBox);
		}
	}

	return out;
}


void NLMeanFilter::UpdateError(float *ImgVar, BlockedArray<Pixel> *pixels, float *ImgErr){
	
	for (int pix = 0; pix < nPixs; pix++) {
        /*fltErr[pix]._pix = pix;
        */
        // The pixel error
        float pixel_error = (ImgVar[pix*3] + ImgVar[pix*3 +1] + ImgVar[pix*3 +2]) / 3.f;
        
        // Set the pixel error info
		ImgErr[pix] = pixel_error / (*pixels)(pix%_xPixelCount, pix/_xPixelCount)._nSamplesBox;
    }
}

void NLMeanFilter::GetSamplingMaps(int spp, int nSamples, float *mapA, float *mapB) {
    // Initialize the sampling maps using the pixel costs
    //mapA.resize(_nPix); 
	//mapB.resize(_nPix);

    for (int pix = 0; pix < nPixs; pix++) {
        mapA[pix] = ImgErr_A[pix] + ImgErr_B[pix];
    }
    
    // Blur the maps a little
    //Gauss2D gauss(.8f, KERNEL_NORM_UNIT);
    //gauss.Apply(_xPixelCount, _yPixelCount, mapA, _tmpMap, mapA);
    
    // Normalize them to sum up to nSamples/2 each
    //float sumA = accumulate(mapA.begin(), mapA.end(), 0.f);
	float sumA = 0.0f;
	for (int index = 0; index < nPixs; index++ )
	{
		sumA += mapA[index] + mapB[index];
	}

    float nSamplesA = nSamples / 2.f;
    for (int pix = 0; pix < nPixs; pix++) {
        mapA[pix] *= nSamplesA / sumA;
    }
    
    // Clamp the map to "lim" samples per pixel max. "lim" is set to spp-1, so
    // that, even with error propagation, no more than spp can be picked
    int nPixOver1, nPixOver2;
    int lim = spp;
    do {
        nPixOver1 = 0;
        for (int pix = 0; pix < nPixs; pix++) {
            if (mapA[pix] > lim) {
                mapA[pix] = lim;
                nPixOver1 += 1;
            }
        }
        // Redistribute the remaining budget over the map
        nPixOver2 = 0;
        //float distributed = accumulate(mapA.begin(), mapA.end(), 0.f);
		float distributed = 0.0f;
		for (int index = 0; index < nPixs; index++ )
		{
			distributed += mapA[index] + mapB[index];
		}

        float scale = (nSamplesA-nPixOver1*lim)/(distributed-nPixOver1*lim);
        if (scale < 0) {
            Severe("Negative scale in sample redistribution!");
        }
        for (int pix = 0; pix < nPixs; pix++) {
            if (mapA[pix] < lim)
                mapA[pix] *= scale;
            
            if (mapA[pix] > lim)
                nPixOver2 += 1;
        }
    } while(nPixOver2 > 0);
    
    //copy(mapA.begin(), mapA.end(), mapB.begin());
    /*if (PbrtOptions.verbose) {
        DumpMap(mapA, "map", DUMP_ITERATION);
    }*/
}


float* NLMeanFilter::Cal(BlockedArray<Pixel> *pixels, int xPixelCount, int yPixelCount)
{
	int nPix = xPixelCount * yPixelCount;

	float *outRGB = new float[3*nPix];
	//int offset;
	//int xMin, xMax, yMin, yMax;

	// cal mean & variance   of pixelA  pixelB
	float *mean = new float[3];
	float *var  = new float[nPix*3];
	int n, index;
	for (int y = 0; y < yPixelCount; y++)
	{
		for (int x = 0; x < xPixelCount; x++)
		{
			index = y*xPixelCount + x;
			
			n = (*pixels)(x, y)._nSamplesBox;
			mean[0] = (*pixels)(x, y)._LrgbSumBox[0] / n;
			mean[1] = (*pixels)(x, y)._LrgbSumBox[1] / n;
			mean[2] = (*pixels)(x, y)._LrgbSumBox[2] / n;

			var[index*3   ] = max( 0.f, ((*pixels)(x, y)._LrgbSumSqrBox[0] -  (*pixels)(x, y)._LrgbSumBox[0]*mean[0]) / (n+1)) ;
			var[index*3 +1] = max( 0.f, ((*pixels)(x, y)._LrgbSumSqrBox[1] -  (*pixels)(x, y)._LrgbSumBox[1]*mean[1]) / (n+1)) ;
			var[index*3 +2] = max( 0.f, ((*pixels)(x, y)._LrgbSumSqrBox[2] -  (*pixels)(x, y)._LrgbSumBox[2]*mean[2]) / (n+1)) ;
		}
	}



	//float *dis = new float[nPix*3];
	/*float *meanA = new float[nPix*3];
	float *meanB = new float[nPix*3];
	

	float patchArea = (2*f+1)*(2*f+1);*/
	

	////----------- calculate Var[p] --------------------//
	//for (int y = 0; y < yPixelCount ; y++)
	//{
	//	for (int x = 0; x < xPixelCount; x++)
	//	{
	//		int index = x + y*xPixelCount;
	//		boxVar[3*pix+0] = max(0.f, (pixel._LrgbSumSqrBox[0] - pixel._LrgbSumBox[0]*mean[0])) / (n-1);
 //           boxVar[3*pix+1] = max(0.f, (pixel._LrgbSumSqrBox[1] - pixel._LrgbSumBox[1]*mean[1])) / (n-1);
 //           boxVar[3*pix+2] = max(0.f, (pixel._LrgbSumSqrBox[2] - pixel._LrgbSumBox[2]*mean[2])) / (n-1);

	//		//meanA[index*3     ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, RGBSumareaTableA, 0) / patchArea;
	//		//meanA[index*3 + 1 ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, RGBSumareaTableA, 1) / patchArea;
	//		//meanA[index*3 + 2 ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, RGBSumareaTableA, 2) / patchArea;
	//		//
	//		//meanB[index*3     ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, RGBSumareaTableB, 0) / patchArea;
	//		//meanB[index*3 + 1 ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, RGBSumareaTableB, 1) / patchArea;
	//		//meanB[index*3 + 2 ] = SummedAreaValue(x , y, xPixelCount, yPixelCount, RGBSumareaTableB, 2) / patchArea;

	//		var[index*3 ] = 0;
	//		var[index*3 + 1] = 0;
	//		var[index*3 + 2] = 0;

	//		
	//		for (int py = -r; py <= r ; py++)
	//		{
	//			for (int px = -r; px <= r ; px++)
	//			{ 
	//				if( (x+px) >= 0 && (y+py) >= 0 && (x+px) < xPixelCount && (y+py) < yPixelCount)
	//				{
	//					var[index*3     ] += rgb[(x+px + (y+py)*xPixelCount)*3    ] - mean[index*3];
	//					var[index*3 + 1 ] += rgb[(x+px + (y+py)*xPixelCount)*3 + 1] - mean[index*3+1];
	//					var[index*3 + 2 ] += rgb[(x+px + (y+py)*xPixelCount)*3 + 2] - mean[index*3+2];
	//				}					
	//			}
	//		}
	//		var[index*3     ] /= patchArea;
	//		var[index*3 + 1 ] /= patchArea;
	//		var[index*3 + 2 ] /= patchArea;
	//		
	//		
	//	}
	//}

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
							
										dis += (pow((*pixels)(x+px, y+py)._Lrgb[i] - (*pixels)(qx+px, qy+py)._Lrgb[i] , 2) 
													- alpha*(var[indexP] - m) ) / (epslon+ k*k*(var[indexP] + var[indexQ]));													
									}
								}
							}
						}
					
						dis /= (3*(2*f+1)*(2*f+1));
						weight = exp(-1* max(0.f, dis));
						totalW += weight;

						outRGB[index*3    ] += (*pixels)(qx, qy)._Lrgb[0] * weight;
						outRGB[index*3 + 1] += (*pixels)(qx, qy)._Lrgb[1] * weight;
						outRGB[index*3 + 2] += (*pixels)(qx, qy)._Lrgb[2] * weight;
					}
				}
			}
			outRGB[index*3    ] /= totalW ; 
			outRGB[index*3 + 1] /= totalW ; 
			outRGB[index*3 + 2] /= totalW ; 

		}
	}

	delete[] mean;
	delete[] var;

	return outRGB;
}

NLMeanFilter *CreateNLMeanFilter(const ParamSet &ps) {
    // Find common filter parameters
    float r = ps.FindOneFloat("rwidth", 7.f);
    float f = ps.FindOneFloat("fwidth", 3.f);
    float k = ps.FindOneFloat("kvalue", 0.45f);
	float alpha = ps.FindOneFloat("alpha", 0.5f);
    return new NLMeanFilter(r, f, k, alpha);
}


