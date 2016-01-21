/* 
 * File:   dualfilm.cpp
 * Author: rousselle
 * 
 * Created on March 13, 2012, 11:35 AM
 */

#include "stdafx.h"
#include "dualfilm.h"
#include <iostream>
#include <fstream>
#include "spectrum.h"
#include "parallel.h"
#include "imageio.h"
#include "image.h"
#include <limits>
#include <vector>
#include <paramSet.h>
#include <string.h>


// DualFilm Method Definitions
DualFilm::DualFilm(int xres, int yres, Filter *filt, const float crop[4],
    const string &fn, bool openWindow, int wnd_rad, float k, int ptc_rad)
    : Film(xres, yres) {
    filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    filename = fn;
    // Compute film image extent
    xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
    xPixelCount = max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
    yPixelCount = max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);

    // Allocate film image storage
    /*int nPix = xPixelCount * yPixelCount;
    pixelsA.resize(nPix);   
	pixelsB.resize(nPix);*/
	pixelsA = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
	pixelsB = new BlockedArray<Pixel>(xPixelCount, yPixelCount);

    // Precompute filter weight table
#define FILTER_TABLE_SIZE 16
    filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
    float *ftp = filterTable;
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
        float fy = ((float)y + .5f) *
                   filter->yWidth / FILTER_TABLE_SIZE;
        for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
            float fx = ((float)x + .5f) *
                       filter->xWidth / FILTER_TABLE_SIZE;
            *ftp++ = filter->Evaluate(fx, fy);
        }
    }

    // Possibly open window for image display
    if (openWindow || PbrtOptions.openWindow) {
        Warning("Support for opening image display window not available in this build.");
    }
	
	NLmean = (NLMeanFilter*)filter;

    //// Allocate subpixel film image storage
    //subPixelRes = 4;
    //int nSubPix = subPixelRes * xPixelCount * subPixelRes * yPixelCount;
    //subPixelsA.resize(nSubPix);   subPixelsB.resize(nSubPix);

    //// Create the appropriate denoiser
    //_denoiser = new NlmeansDenoiser(xPixelCount, yPixelCount, filt,
    //    filename, wnd_rad, k, ptc_rad, subPixelRes);
}


void DualFilm::AddSample(const CameraSample &sample, const Spectrum &L, TargetBuffer target) {
    // Compute sample's raster extent
    float dimageX = sample.imageX - 0.5f;
    float dimageY = sample.imageY - 0.5f;
    int x0 = Ceil2Int (dimageX - filter->xWidth);
    int x1 = Floor2Int(dimageX + filter->xWidth);
    int y0 = Ceil2Int (dimageY - filter->yWidth);
    int y1 = Floor2Int(dimageY + filter->yWidth);
    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);
    if ((x1-x0) < 0 || (y1-y0) < 0)
    {
        PBRT_SAMPLE_OUTSIDE_IMAGE_EXTENT(const_cast<CameraSample *>(&sample));
        return;
    }

    // Loop over filter support and add sample to pixel arrays
    float xyz[3];
    L.ToXYZ(xyz);
	float rgb[3];
	L.ToRGB(rgb);

    // Precompute $x$ and $y$ filter table offsets
    int *ifx = ALLOCA(int, x1 - x0 + 1);
    for (int x = x0; x <= x1; ++x) {
        float fx = fabsf((x - dimageX) *
                         filter->invXWidth * FILTER_TABLE_SIZE);
        ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
    }
    int *ify = ALLOCA(int, y1 - y0 + 1);
    for (int y = y0; y <= y1; ++y) {
        float fy = fabsf((y - dimageY) *
                         filter->invYWidth * FILTER_TABLE_SIZE);
        ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
    }

    // Select the right target buffer
    //vector<NlmeansPixel> &pixels = (target == BUFFER_A) ? pixelsA : pixelsB;
    BlockedArray<Pixel> *pixels = (target == BUFFER_A) ?pixelsA : pixelsB;


    // Always use AtomicAdd since adaptive sampling might be using large kernels
    bool syncNeeded = true; // (filter->xWidth > 0.5f || filter->yWidth > 0.5f);
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            // Evaluate filter value at $(x,y)$ pixel
            int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
            float filterWt = filterTable[offset];
            
            // Update pixel values with filtered sample contribution
            /*int pix = xPixelCount * (y - yPixelStart) + x - xPixelStart;
            NlmeansPixel &pixel = pixels[pix];*/
			
			Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
			
            if (!syncNeeded) {
                pixel.Lxyz[0] += filterWt * xyz[0];
                pixel.Lxyz[1] += filterWt * xyz[1];
                pixel.Lxyz[2] += filterWt * xyz[2];
                pixel.weightSum += filterWt;
            }
            else {
                // Safely update _Lrgb_ and _weightSum_ even with concurrency
                AtomicAdd(&pixel.Lxyz[0], filterWt * xyz[0]);
                AtomicAdd(&pixel.Lxyz[1], filterWt * xyz[1]);
                AtomicAdd(&pixel.Lxyz[2], filterWt * xyz[2]);
                AtomicAdd(&pixel.weightSum, filterWt);
            }
        }
    }
    
    //// We're done is, this sample is outside the film buffer
    int x = Floor2Int(sample.imageX);
    int y = Floor2Int(sample.imageY);
    if (x < 0 || y < 0 || x >= xPixelCount || y >= yPixelCount)
        return;
    
    // Store variance information
    //int pix = xPixelCount * (y - yPixelStart) + x - xPixelStart;
    Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
    AtomicAdd(&pixel._LrgbSumBox[0], rgb[0]);
    AtomicAdd(&pixel._LrgbSumBox[1], rgb[1]);
    AtomicAdd(&pixel._LrgbSumBox[2], rgb[2]);
    AtomicAdd(&pixel._LrgbSumSqrBox[0], rgb[0]*rgb[0]);
    AtomicAdd(&pixel._LrgbSumSqrBox[1], rgb[1]*rgb[1]);
    AtomicAdd(&pixel._LrgbSumSqrBox[2], rgb[2]*rgb[2]);
    AtomicAdd((AtomicInt32*)&pixel._nSamplesBox, (int32_t)1);
    
    //// Store on sub-pixel grid
    //x = Floor2Int(subPixelRes * sample.imageX);
    //y = Floor2Int(subPixelRes * sample.imageY);
    //int pix = x + y * subPixelRes * xPixelCount;
    //Pixel subpix = (target == BUFFER_A) ? subPixelsA[pix] : subPixelsB[pix];
    //AtomicAdd(&subpix._LrgbSumBox[0], rgb[0]);
    //AtomicAdd(&subpix._LrgbSumBox[1], rgb[1]);
    //AtomicAdd(&subpix._LrgbSumBox[2], rgb[2]);
    //AtomicAdd(&subpix._nSamplesBox, 1);
}

void DualFilm::Splat(const CameraSample &sample, const Spectrum &L){
}


void DualFilm::GetSampleExtent(int *xstart, int *xend,
                                int *ystart, int *yend) const {
    *xstart = Floor2Int(xPixelStart + 0.5f - filter->xWidth);
    *xend   = Floor2Int(xPixelStart + 0.5f + xPixelCount  +
                        filter->xWidth);

    *ystart = Floor2Int(yPixelStart + 0.5f - filter->yWidth);
    *yend   = Floor2Int(yPixelStart + 0.5f + yPixelCount +
                        filter->yWidth);
}


void DualFilm::GetPixelExtent(int *xstart, int *xend,
                               int *ystart, int *yend) const {
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}


void DualFilm::WriteImage(float splatScale) {
    // Dump the noisy data by combining the two buffers
    /*ImageBuffer rgb(3 * xPixelCount * yPixelCount);
    ImageBuffer rgbA(3 * xPixelCount * yPixelCount);
    ImageBuffer rgbB(3 * xPixelCount * yPixelCount);*/
	
	int nPix = xPixelCount * yPixelCount;
    
	float *rgbA = new float[3*nPix];
	float *rgbB = new float[3*nPix];
    float *rgb = new float[3*nPix];
    int offset = 0;

    for (int y = 0, pix = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x, ++pix) {
            // Convert pixel XYZ color to RGB
            XYZToRGB((*pixelsA)(x, y).Lxyz, &rgbA[3*offset]);

            // Normalize pixel with weight sum
            Pixel &pixelA = (*pixelsA)(x, y);
            float wgtSumA = pixelA.weightSum;
            if (wgtSumA != 0.f) {
                    float invWt = 1.f / wgtSumA;
                    rgbA[3*offset  ] = max(0.f, rgbA[3*offset  ] * invWt);
                    rgbA[3*offset+1] = max(0.f, rgbA[3*offset+1] * invWt);
                    rgbA[3*offset+2] = max(0.f, rgbA[3*offset+2] * invWt);
            }
			
			XYZToRGB((*pixelsB)(x, y).Lxyz, &rgbB[3*offset]);

            // Normalize pixel with weight sum
            Pixel &pixelB = (*pixelsB)(x, y);
            float wgtSumB = pixelB.weightSum;
            if (wgtSumB != 0.f) {
                    float invWt = 1.f / wgtSumB;
                    rgbB[3*offset  ] = max(0.f, rgbB[3*offset  ] * invWt);
                    rgbB[3*offset+1] = max(0.f, rgbB[3*offset+1] * invWt);
                    rgbB[3*offset+2] = max(0.f, rgbB[ 3*offset+2] * invWt);
            }
            
            // Sum both buffers
            float wgtSum = wgtSumA + wgtSumB;
            if (wgtSumA != 0.f) {
                    rgb[3*offset  ] = (wgtSumA * rgbA[3*offset  ] + wgtSumB * rgbB[3*offset  ]) / wgtSum;
                    rgb[3*offset+1] = (wgtSumA * rgbA[3*offset+1] + wgtSumB * rgbB[3*offset+1]) / wgtSum;
                    rgb[3*offset+2] = (wgtSumA * rgbA[3*offset+2] + wgtSumB * rgbB[3*offset+2]) / wgtSum;
            }
        }
    }

	
	

    ::WriteImage(filename, &rgb[0], NULL, xPixelCount, yPixelCount, xPixelCount, yPixelCount, 0, 0);


	rgb = NLmean->NLFiltering(pixelsA, pixelsB, xPixelCount, yPixelCount);

    // Filter out noise from data and store the result
    /*
	if (!NLmean->IsReady())
			_denoiser->UpdatePixelData(pixelsA, pixelsB, subPixelsA, subPixelsB, NLM_DATA_FINAL);
	*/
    
}


void DualFilm::UpdateDisplay(int x0, int y0, int x1, int y1,
    float splatScale) {
}


DualFilm *CreateDualFilm(const ParamSet &params, Filter *filter) {
    string filename = params.FindOneString("filename", PbrtOptions.imageFile);
    if (filename == "")
#ifdef PBRT_HAS_OPENEXR
        filename = "pbrt.exr";
#else
        filename = "pbrt.tga";
#endif

    int xres = params.FindOneInt("xresolution", 640);
    int yres = params.FindOneInt("yresolution", 480);
    if (PbrtOptions.quickRender) xres = max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = max(1, yres / 4);
    bool openwin = params.FindOneBool("display", false);
    float crop[4] = { 0, 1, 0, 1 };
    int cwi;
    const float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
        crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
        crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
        crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
    }

    ///////////////////////
    // Parse denoiser flags

    // Filter width
    int wnd_rad = params.FindOneInt("wnd_rad", 10);
    
    // Damping parameter
    float k = params.FindOneFloat("k", .1f);
    
    // Patchsize
    int ptc_rad = params.FindOneInt("ptc_rad", 3);

    return new DualFilm(xres, yres, filter, crop, filename, openwin, wnd_rad,
        k, ptc_rad);
}


