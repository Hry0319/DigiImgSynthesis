
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


vector<MedianCutRect*> leafs;

// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    //delete distribution;
    //delete radianceMap;
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

    summedArea = (float*)malloc(sizeof(float)*width*height);
    for (int i=0;i<width*height;i++){
        summedArea[i]=0;
    }
    //float rgb[3];
    for(int index = 0; index < width*height; index++ ){
        //texels[index].ToRGB(rgb);
        //summedArea[index] = rgb[0]+rgb[1]+rgb[2];
        
        float vp = (float)(index/width) / (float)height;
        float sinTheta = sinf(M_PI * float((index/width)+.5f)/float(height));
        summedArea[index] = texels[index].y() * sinTheta;
    }

    AllocateSummedAreaTable(summedArea,width,height);

    MedianCutRect *mcr = new MedianCutRect(NULL,NULL,0,0,width,height,summedArea[width*height]); 

    CutCut(mcr, 0);
    
    for(unsigned int index = 0; index < leafs.size(); index++)
    {
        int area = (leafs[index]->Rect_width - leafs[index]->x) * (leafs[index]->Rect_height - leafs[index]->y);
        leafs[index]->rgbspectrum = SummedAreaValue(leafs[index]->x, leafs[index]->y, leafs[index]->Rect_width, leafs[index]->Rect_height);
        leafs[index]->SummedValue = leafs[index]->rgbspectrum.y();
        /*leafs[index]->rgbspectrum.ToRGB(rgb);
        rgb[0] = rgb[0]/area;
        rgb[1] = rgb[1]/area;
        rgb[2] = rgb[2]/area;*/
        CaculateLights(leafs[index], texels);
    }
}

void MedianCutEnvironmentLight::CaculateLights(MedianCutRect *mcr, RGBSpectrum *texels)const
{
    int area = (mcr->Rect_width-mcr->x)*(mcr->Rect_height-mcr->y);
    for (int y = mcr->y; y < mcr->Rect_height;y++)
    {
        for (int x = mcr->x; x < mcr->Rect_width;x++)
        {
            float rgb[3];
            texels[x+y*AreaWidth].ToRGB(rgb);
            if(rgb[0] > 0)
                mcr->meanRGB[0] += rgb[0];
            if(rgb[1] > 0)
                mcr->meanRGB[1] += rgb[1];
            if(rgb[2] > 0)
                mcr->meanRGB[2] += rgb[2];
        }
    }
    if(mcr->meanRGB[0] > 0)
        mcr->meanRGB[0] /= area;
    if(mcr->meanRGB[1] > 0)
        mcr->meanRGB[1] /= area;
    if(mcr->meanRGB[2] > 0)
        mcr->meanRGB[2] /= area;    
        
    const float phi = (mcr->Rect_width * 0.5f + mcr->x) / float(AreaWidth) * 2.f * M_PI;
    const float theta = (mcr->Rect_height * 0.5f + mcr->y) / float(AreaHeight) * M_PI;
    // Add Light
    const float costheta = cosf(theta), sintheta = sinf(theta);
    const float cosphi = cosf(phi), sinphi = sinf(phi);
    float px = sintheta * cosphi, py = sintheta * sinphi, pz = costheta;
    mcr->LightPoint.x = px;
    mcr->LightPoint.y = py;
    mcr->LightPoint.z = pz;
    
}

//RGBSpectrum SummedAreaValue(int x,int y,int width,int height)
void MedianCutEnvironmentLight::CutCut( MedianCutRect *root, int nowTreeHeight/*,int x,int y,int width,int height*/ )const
{
    if (nowTreeHeight > Log2(nSamples))
    {
        return;
    }
    else
    {
        if(root->Rect_width > root->Rect_height)
        {
            float variance = SummedAreaValue(root->x, root->y, root->Rect_width, root->Rect_height).y();
            int   medianCC = 0;
            float leftRect, RightRect;
            //RGBSpectrum  minL,minR;
            for(int indexX = 1 ; indexX < root->Rect_width - 1; indexX++)
            {
                leftRect  = SummedAreaValue(
                                root->x, 
                                root->y, 
                                indexX /*root->Rect_width*/, 
                                root->Rect_height
                                ).y();
                RightRect = SummedAreaValue(
                                root->x + indexX,
                                root->y, 
                                root->Rect_width - indexX -1,
                                root->Rect_height
                                ).y();
                if (variance > abs(leftRect-RightRect))
                {
                    medianCC = indexX;
                    variance = abs(leftRect-RightRect);

                }
            }
            if(nowTreeHeight + 1 > Log2(nSamples))
            {
                leafs.push_back((MedianCutRect*)root);
            }
            else
            {                
                root->left  = new MedianCutRect(NULL, NULL, root->x, root->y, medianCC, root->Rect_height, 0);
                CutCut((MedianCutRect*)root->left, nowTreeHeight + 1);
                root->right = new MedianCutRect(NULL, NULL, root->x + medianCC, root->y, root->Rect_width - medianCC -1, root->Rect_height, 0);
                CutCut((MedianCutRect*)root->right, nowTreeHeight + 1);
            }

        }
        else
        {
            float variance = SummedAreaValue(root->x, root->y, root->Rect_width, root->Rect_height).y();
            int   medianCC = 0;
            float UpRect, ButtonRect;
            for(int indexY = 1; indexY < root->Rect_height - 1; indexY++)
            {
                UpRect  = SummedAreaValue(
                                root->x, 
                                root->y, 
                                root->Rect_width, 
                                indexY
                                ).y();
                ButtonRect = SummedAreaValue(
                                root->x,
                                root->y + indexY, 
                                root->Rect_width,
                                root->Rect_height - indexY - 1
                                ).y();
                if (variance > abs(UpRect-ButtonRect))
                {
                    medianCC = indexY;
                    variance = abs(UpRect-ButtonRect);
                }
            }
            if(nowTreeHeight + 1 > Log2(nSamples))
            {
                leafs.push_back((MedianCutRect*)root);
            }
            else
            {
                root->left  = new MedianCutRect(NULL, NULL, root->x, root->y, root->Rect_width, medianCC, 0);
                CutCut((MedianCutRect*)root->left, nowTreeHeight + 1);
                root->right = new MedianCutRect(NULL, NULL, root->x, root->y + medianCC, root->Rect_width, root->Rect_height - medianCC - 1, 0);
                CutCut((MedianCutRect*)root->right, nowTreeHeight + 1);
            }
        }
    }
}

void MedianCutEnvironmentLight::AllocateSummedAreaTable(float *summedArea, unsigned int width, unsigned int height)const
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
                summedArea[_i_j] = summedArea[_i_j] + summedArea[_i_j-1] + summedArea[_i_j-width] -summedArea[_i_j-width-1]; 
            }
        }
    }
}

Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon, const LightSample &ls, float time, Vector *wi, float *pdf, VisibilityTester *visibility) const 
{
    //PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    //// Find $(u,v)$ sample coordinates in infinite light texture
    //float uv[2], mapPdf;
    ////distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    ////if (mapPdf == 0.f) return 0.f;

    //// Convert infinite light sample point to direction
    //float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    //float costheta = cosf(theta), sintheta = sinf(theta);
    //float sinphi = sinf(phi), cosphi = cosf(phi);
    //*wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
    //                          costheta));

    //// Compute PDF for sampled infinite light direction
    ////*pdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    //if (sintheta == 0.f) *pdf = 0.f;

    //// Return radiance value for infinite light direction
    //visibility->SetRay(p, pEpsilon, *wi, time);
    //Spectrum Ls = Spectrum(radianceMap->Lookup(uv[0], uv[1]),
    //                       SPECTRUM_ILLUMINANT);
    //PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    //return Ls;

    int randomLightSampleNum = rand()%nSamples;

    *wi = LightToWorld(Vector(leafs[randomLightSampleNum]->LightPoint));
    visibility->SetRay(p, pEpsilon, *wi, time); 
    *pdf = 1.0f;
    Spectrum Ls = Spectrum(leafs[randomLightSampleNum]->rgbspectrum);
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
                            float *pdf) const 
{
    //PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    //// Compute direction for infinite light sample ray

    //// Find $(u,v)$ sample coordinates in infinite light texture
    //float uv[2], mapPdf;
    //distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    //if (mapPdf == 0.f) return Spectrum(0.f);

    //float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    //float costheta = cosf(theta), sintheta = sinf(theta);
    //float sinphi = sinf(phi), cosphi = cosf(phi);
    //Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
    //                                costheta));
    //*Ns = (Normal)d;

    //// Compute origin for infinite light sample ray
    //Point worldCenter;
    //float worldRadius;
    //scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    //Vector v1, v2;
    //CoordinateSystem(-d, &v1, &v2);
    //float d1, d2;
    //ConcentricSampleDisk(u1, u2, &d1, &d2);
    //Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    //*ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    //// Compute _MedianCutEnvironmentLight_ ray PDF
    //float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    //float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    //*pdf = directionPdf * areaPdf;
    //if (sintheta == 0.f) *pdf = 0.f;
    //Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    //PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();

    //return Ls;

    int randomLightSampleNum = rand()%nSamples;
    
    Point lightPos = LightToWorld(leafs[randomLightSampleNum]->LightPoint);
    *ray = Ray(lightPos, - Normalize(Vector(lightPos)), 0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    *pdf = UniformSpherePdf();
    Spectrum Ls = Spectrum(leafs[randomLightSampleNum]->rgbspectrum);
    return Ls;
}

float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const 
{
    /*PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;*/
    return 1.f/nSamples;
    //return 0.0f;
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

Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const 
{
    //Point worldCenter;
    //float worldRadius;
    //scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);

    //return M_PI * worldRadius * worldRadius *
    //    Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
    
    return Spectrum(SummedAreaValue(0, 0,AreaWidth, AreaHeight), SPECTRUM_ILLUMINANT);
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
//        float time, RNG &rng, Spectrum *coeffs) const {
//    // Project _MedianCutEnvironmentLight_ to SH using Monte Carlo if visibility needed
//    if (computeLightVis) {
//        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
//                         time, rng, coeffs);
//        return;
//    }
//    for (int i = 0; i < SHTerms(lmax); ++i)
//        coeffs[i] = 0.f;
//    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
//    if (min(ntheta, nphi) > 50) {
//        // Project _MedianCutEnvironmentLight_ to SH from lat-long representation
//
//        // Precompute $\theta$ and $\phi$ values for lat-long map projection
//        float *buf = new float[2*ntheta + 2*nphi];
//        float *bufp = buf;
//        float *sintheta = bufp;  bufp += ntheta;
//        float *costheta = bufp;  bufp += ntheta;
//        float *sinphi = bufp;    bufp += nphi;
//        float *cosphi = bufp;
//        for (int theta = 0; theta < ntheta; ++theta) {
//            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
//            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
//        }
//        for (int phi = 0; phi < nphi; ++phi) {
//            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
//            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
//        }
//        float *Ylm = ALLOCA(float, SHTerms(lmax));
//        for (int theta = 0; theta < ntheta; ++theta) {
//            for (int phi = 0; phi < nphi; ++phi) {
//                // Add _MedianCutEnvironmentLight_ texel's contribution to SH coefficients
//                Vector w = Vector(sintheta[theta] * cosphi[phi],
//                                  sintheta[theta] * sinphi[phi],
//                                  costheta[theta]);
//                w = Normalize(LightToWorld(w));
//                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
//                                       SPECTRUM_ILLUMINANT);
//                SHEvaluate(w, lmax, Ylm);
//                for (int i = 0; i < SHTerms(lmax); ++i)
//                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
//                        (M_PI / ntheta) * (2.f * M_PI / nphi);
//            }
//        }
//
//        // Free memory used for lat-long theta and phi values
//        delete[] buf;
//    }
//    else {
//        // Project _MedianCutEnvironmentLight_ to SH from cube map sampling
//        SHProjectCube(InfiniteAreaCube(this, scene, time, computeLightVis,
//                                       pEpsilon),
//                      p, 200, lmax, coeffs);
//    }
//}