
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

#ifndef PBRT_LIGHTS_MEDIAN_CUT_ENV_H
#define PBRT_LIGHTS_MEDIAN_CUT_ENV_H

// lights/infinite.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"
#include <vector>

class MedianCutRect {
    public: 
        MedianCutRect(void* _left, void* _right, int _x, int _y, int _Rect_width, int _Rect_height, RGBSpectrum _rgbspectrum)
        {
            left        = _left;
            right       = _right;
            x           = _x;
            y           = _y;
            Rect_width  = _Rect_width;
            Rect_height = _Rect_height;
            SummedRGB   = _rgbspectrum; 
            SummedLum   = SummedRGB.y();
            //rgbspectrum.ToRGB(meanRGB);
        } 
        
        ~MedianCutRect()
        {
            /*delete this->left;
            delete this->right;
            delete this;*/
        }
        
        void*       left;
        void*       right;
        int         x,y;
        int         Rect_width,Rect_height;        
        //float       meanRGB[3];
        RGBSpectrum SummedRGB;
        float       SummedLum;

        //centroid Point
        Point       LightPoint;    

        static void deleteMcr(MedianCutRect *mcr)
        {
            if (mcr->left == NULL && mcr->right == NULL)
            {
                delete mcr;
                return;
            }
            //if (mcr->left != NULL)
                deleteMcr((MedianCutRect*)mcr->left);
            //if (mcr->right != NULL)
                deleteMcr((MedianCutRect*)mcr->right);
        }
};

class MedianCutEnvironmentLight : public Light {
    public:
        MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
            const string &texmap);
        ~MedianCutEnvironmentLight();
        
        Spectrum Power(const Scene *) const;
        //Spectrum Le(const RayDifferential &r) const;
        Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls, float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
        Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2, float time, Ray *ray, Normal *Ns, float *pdf) const;
        
        float Pdf(const Point &, const Vector &) const;   
        bool IsDeltaLight() const { return false; }
        
        //void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        //                bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;

    private:
        //MIPMap<RGBSpectrum> *radianceMap;
        Distribution2D *distribution;
        RGBSpectrum     *summedArea;
        MedianCutRect   *mcr;
        int             nSamples;
        int             AreaWidth, AreaHeight;
        
        void InitSummedAreaTable(RGBSpectrum *summedArea, unsigned int width, unsigned int height)const;
        void CutCut(MedianCutRect *root, int nowTreeHeight/*,int x,int y,int width,int height*/)const;
        void CaculateLights(MedianCutRect *mcr)const;

        RGBSpectrum SummedAreaValue(int x,int y,int width,int height)const
        {
            int x1 = x + width-1;
            int y1 = y + height-1;
            if(x1 > AreaWidth)x1 = AreaWidth;
            if(y1 > AreaHeight)y1 = AreaHeight;
            if(x1 < 0)x1 = 0;
            if(y1 < 0)y1 = 0;
            RGBSpectrum upper_left   = summedArea[x  + AreaWidth*y ];
            RGBSpectrum upper_right  = summedArea[x1 + AreaWidth*y ];
            RGBSpectrum button_left  = summedArea[x  + AreaWidth*y1];
            RGBSpectrum button_right = summedArea[x1 + AreaWidth*y1];

            return button_right + upper_left - upper_right - button_left;
        }
};

MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world, const ParamSet &paramSet);

#endif // PBRT_LIGHTS_MEDIAN_CUT_ENV_H
