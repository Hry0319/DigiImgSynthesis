
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


// lights/infinite.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"
#include <math.h>


vector<MedianCutRect*> leafs;
int     nowleafs;

MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete summedArea;
}

MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &L, int ns, const string &texmap) : Light(light2world, ns)
{
    int width           = 0;
    int height          = 0;
    RGBSpectrum *texels = NULL;
    nSamples            = ns;

    // Read texel data from _texmap_ into _texels_
    if (texmap != "")
    {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    AreaWidth           = width;
    AreaHeight          = height;

    summedArea = new RGBSpectrum[width*height];
    for (int i=0;i<width*height;i++){
        summedArea[i]=0;
    }
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            summedArea[u + v*width] = texels[u + v*width]/*.y()*/ * sinTheta;
        }
    }
    InitSummedAreaTable(summedArea,width,height);

    MedianCutRect *mcr = new MedianCutRect(NULL, NULL, 0, 0, width, height, summedArea[width*height-1]);
    CutCut(mcr, 0);    
    delete mcr;

    for(unsigned int index = 0; index < leafs.size(); index++)
    {
        CaculateLights(leafs[index]);
    }
}

void MedianCutEnvironmentLight::InitSummedAreaTable(RGBSpectrum *summedArea, unsigned int width, unsigned int height)const
{
    for(unsigned int i = 0; i < height; i++ )
    {
        for(unsigned int j = 0; j < width; j++)
        {
            int _i_j = i*width+j;
            if(i==0)
            {
                if(j!=0)
                    summedArea[_i_j] += summedArea[_i_j - 1];
            }
            else if(j==0)
            {
                if(i!=0)
                    summedArea[_i_j] += summedArea[_i_j - width];
            }
            else
            {
                summedArea[_i_j] = summedArea[_i_j] + summedArea[_i_j-1] + summedArea[_i_j-width] - summedArea[_i_j-width-1];
            }
        }
    }
}

void MedianCutEnvironmentLight::CutCut( MedianCutRect *root, int nowTreeHeight)const
{
    int LimitedHeight = Log2(nSamples);
    if (nowTreeHeight > LimitedHeight)
    {
        return;
    }
    else
    {
        if(root->Rect_width > root->Rect_height)
        {
            int         medianCC = 0;
            RGBSpectrum leftRect, RightRect;

            for(int indexX = 1 ; indexX < root->Rect_width - 1; indexX++)
            {
                leftRect  = SummedAreaValue(
                                root->x,
                                root->y,
                                indexX,
                                root->Rect_height
                                );
                if (leftRect.y() >= root->SummedLum/2)
                {
                    medianCC = indexX;

                    RightRect = SummedAreaValue(
                                    root->x + indexX,
                                    root->y,
                                    root->Rect_width - indexX -1,
                                    root->Rect_height
                                    );
                    break;
                }

            }
            if(nowTreeHeight + 1 > LimitedHeight)
            {
                leafs.push_back((MedianCutRect*)root);
            }
            else
            {
                root->left  = new MedianCutRect(NULL, NULL, root->x, root->y, medianCC, root->Rect_height, leftRect);
                root->right = new MedianCutRect(NULL, NULL, root->x + medianCC, root->y, root->Rect_width - medianCC -1, root->Rect_height, RightRect);

                CutCut((MedianCutRect*)root->left, nowTreeHeight + 1);
                CutCut((MedianCutRect*)root->right, nowTreeHeight + 1);
            }

        }
        else
        {
            int         medianCC = 0;
            RGBSpectrum UpRect, ButtonRect;
            for(int indexY = 1; indexY < root->Rect_height - 1; indexY++)
            {
                UpRect  = SummedAreaValue(
                                root->x,
                                root->y,
                                root->Rect_width,
                                indexY
                                );
                if (UpRect.y() >= root->SummedLum/2)
                {
                    medianCC = indexY;

                    ButtonRect = SummedAreaValue(
                                        root->x,
                                        root->y + indexY,
                                        root->Rect_width,
                                        root->Rect_height - indexY - 1
                                        );
                    break;
                }
            }
            if(nowTreeHeight + 1 > LimitedHeight)
            {
                leafs.push_back((MedianCutRect*)root);
            }
            else
            {
                root->left  = new MedianCutRect(NULL, NULL, root->x, root->y, root->Rect_width, medianCC, UpRect);
                root->right = new MedianCutRect(NULL, NULL, root->x, root->y + medianCC, root->Rect_width, root->Rect_height - medianCC - 1, ButtonRect);

                CutCut((MedianCutRect*)root->left, nowTreeHeight + 1);
                CutCut((MedianCutRect*)root->right, nowTreeHeight + 1);
            }
        }
    }
}


void MedianCutEnvironmentLight::CaculateLights(MedianCutRect *mcr)const
{
    // Light point
    const float phi = (mcr->Rect_width * 0.5f + mcr->x) / float(AreaWidth) * 2.f * M_PI;
    const float theta = (mcr->Rect_height * 0.5f + mcr->y) / float(AreaHeight) * M_PI;
    const float costheta = cosf(theta), sintheta = sinf(theta);
    const float cosphi = cosf(phi), sinphi = sinf(phi);
    float px = sintheta * cosphi, py = sintheta * sinphi, pz = costheta;
    mcr->LightPoint = Point(px,py,pz);
}

Spectrum
MedianCutEnvironmentLight::Sample_L(
                            const Point &p,
                            float pEpsilon,
                            const LightSample &ls,
                            float time,
                            Vector *wi,
                            float *pdf,
                            VisibilityTester *visibility
                            ) const
{
    int randomLightSampleNum = rand()%nSamples;

    *wi = Normalize( LightToWorld(Vector(leafs[randomLightSampleNum]->LightPoint)) );
    visibility->SetRay(p, pEpsilon, *wi, time);
    *pdf = 1.f/*/nSamples*/;
    Spectrum Ls = Spectrum(leafs[randomLightSampleNum]->SummedRGB, SPECTRUM_ILLUMINANT);
    return Ls;
}

Spectrum
MedianCutEnvironmentLight::Sample_L(
                            const Scene *scene,
                            const LightSample &ls,
                            float u1,
                            float u2,
                            float time,
                            Ray *ray,
                            Normal *Ns,
                            float *pdf
                            ) const
{
    int randomLightSampleNum = rand()%nSamples;
    nowleafs = randomLightSampleNum;

    Point lightPos = LightToWorld(leafs[randomLightSampleNum]->LightPoint);
    *ray = Ray(lightPos, - Normalize(Vector(lightPos)), 0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    //*pdf = UniformSpherePdf();
    *pdf = 1.f/nSamples;
    Spectrum Ls = Spectrum(leafs[randomLightSampleNum]->SummedRGB, SPECTRUM_ILLUMINANT);
    return Ls;
}

Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const
{
    return Spectrum(SummedAreaValue(0, 0,AreaWidth, AreaHeight), SPECTRUM_ILLUMINANT);
}

MedianCutEnvironmentLight
*CreateMedianCutEnvironmentLight(
    const Transform &light2world,
    const ParamSet &paramSet
    )
{
    Spectrum L          = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc         = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string   texmap     = paramSet.FindOneFilename("mapname", "");
    int      nSamples   = paramSet.FindOneInt("nsamples", 1);

    if (PbrtOptions.quickRender)
        nSamples = max(1, nSamples / 4);

    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}

float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const
{
    //return 1.f/nSamples;
    return 0.0f;
}

//Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const
//{
//    Vector wh = Normalize(WorldToLight(r.d));
//    float s = SphericalPhi(wh) * INV_TWOPI;
//    float t = SphericalTheta(wh) * INV_PI;
//
//    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
//}

//void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
//        int lmax, const Scene *scene, bool computeLightVis,
//        float time, RNG &rng, Spectrum *coeffs) const
//{
//  for (int i = 0; i < SHTerms(lmax); ++i)
//    coeffs[i] = 0.f;
//
//  Point worldCenter;
//  float worldRadius;
//  scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
//
//  for (int light_id = 0; light_id < leafs.size(); ++light_id) {
//      Point lightPos = LightToWorld(leafs[light_id]->LightPoint);
//    if (computeLightVis &&
//        scene->IntersectP(Ray(p, Normalize(Vector(lightPos)), pEpsilon,
//                              INFINITY, time)))
//      continue;
//    // Project point light source to SH
//    float *Ylm = ALLOCA(float, SHTerms(lmax));
//    Vector wi = Normalize(Vector(lightPos));
//    SHEvaluate(wi, lmax, Ylm);
//    Spectrum Li = Spectrum(leafs[light_id]->SummedRGB * nSamples / 128, SPECTRUM_ILLUMINANT);
//    for (int i = 0; i < SHTerms(lmax); ++i)
//      coeffs[i] += Li * Ylm[i];
//  }
//}
