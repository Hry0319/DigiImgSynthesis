
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_TWOSTAGESSAMPLERRENDERER_H
#define PBRT_RENDERERS_TWOSTAGESSAMPLERRENDERER_H

// renderers/twostages.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
#include "film.h"
#include "integrator.h"

// TwoStagesSamplerRenderer Declarations
class TwoStagesSamplerRenderer : public Renderer {
public:
    // TwoStagesSamplerRenderer Public Methods
    TwoStagesSamplerRenderer(Sampler *s, Camera *c, SurfaceIntegrator *si,
                    VolumeIntegrator *vi, bool visIds);
    ~TwoStagesSamplerRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // TwoStagesSamplerRenderer Private Data
    bool visualizeObjectIds;
    Sampler *sampler;
    Camera *camera;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
};



// TwoStagesSamplerRendererTask Declarations
class TwoStagesSamplerRendererTask : public Task {
public:
    // TwoStagesSamplerRendererTask Public Methods
    TwoStagesSamplerRendererTask(const Scene *sc, TwoStagesSamplerRenderer *ren, Camera *c,
                        ProgressReporter &pr, Sampler *ms, Sample *sam,
                        bool visIds, int tn, int tc, bool isDual = false)
      : reporter(pr)
    {
        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc;
        dualSampler = isDual;
    }
    void Run();
private:
    // TwoStagesSamplerRendererTask Private Data
    bool dualSampler;
    const Scene *scene;
    const TwoStagesSamplerRenderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
};



#endif // PBRT_RENDERERS_TWOSTAGESSAMPLERRENDERER_H
