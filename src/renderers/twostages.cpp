
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


// renderers/twostages.cpp*
#include "stdafx.h"
#include "renderers/twostages.h"
#include "scene.h"
#include "film/dualfilm.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "samplers/dualsampler.h"

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t   hash, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

// TwoStagesSamplerRendererTask Definitions
void TwoStagesSamplerRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _TwoStagesSamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    //
    DualFilm *dualfilm = dynamic_cast<DualFilm*> (camera->film);
    bool singleBuffered = (dualfilm == NULL);
    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (visualizeObjectIds) {
                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                    // random shading based on shape id...
                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash((char *)ids, sizeof(ids));
                    float rgb[3] = { (h & 0xff), (h >> 8) & 0xff, (h >> 16) & 0xff };
                    Ls[i] = Spectrum::FromRGB(rgb);
                    Ls[i] /= 255.f;
                }
                else
                    Ls[i] = 0.f;
            }
            else {
            if (rayWeight > 0.f) {
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i]);
            }
            else {
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }

            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned"
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned"
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            if (singleBuffered) {
                for (int i = 0; i < sampleCount; ++i)
                {
                    PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                    camera->film->AddSample(samples[i], Ls[i]);
                    PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
                }
            }
            else {
                TargetBuffer target = (taskNum < taskCount/2) ? BUFFER_A : BUFFER_B;
                for (int i = 0; i < sampleCount; ++i)
                {
                    PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                    dualfilm->AddSample(samples[i], Ls[i], target);
                    PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
                }
            }
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _TwoStagesSamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// TwoStagesSamplerRenderer Method Definitions
TwoStagesSamplerRenderer::TwoStagesSamplerRenderer(Sampler *s, Camera *c,
    SurfaceIntegrator *si, VolumeIntegrator *vi, bool visIds) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


TwoStagesSamplerRenderer::~TwoStagesSamplerRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}

//
// call by api.c => pbrtWorldEnd()
//
void TwoStagesSamplerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

    // Create and launch _TwoStagesSamplerRendererTask_s for rendering image

    DualSampler *dualSampler = dynamic_cast<DualSampler *> (sampler);

    if (dualSampler != NULL) {
        // Initialization phase using uniform sampling
        // Compute number of _TwoStagesSamplerRendererTask_s to create for rendering
        int nPixels = camera->film->xResolution * camera->film->yResolution;
        int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
        nTasks = RoundUpPow2(nTasks);
        // Distribute the workload
        ProgressReporter reporter(nTasks, "Rendering");
        vector<Task *> renderTasks;
        for (int i = 0; i < nTasks; ++i)
            renderTasks.push_back(new TwoStagesSamplerRendererTask(scene, this,
                camera, reporter, sampler, sample, visualizeObjectIds, nTasks-1-i, nTasks, true));
        EnqueueTasks(renderTasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < renderTasks.size(); ++i)
            delete renderTasks[i];
        renderTasks.clear();
        if (dualSampler->PixelsToSampleTotal() == 0)
            dualSampler->Finalize();
        reporter.Done();

        // Adaptive phase. If the there is no pixels left to sample, the user
        // requested uniform sampling and we skip the adaptive phase.
        if (dualSampler->PixelsToSampleTotal() > 0) {
            dualSampler->SetAdaptiveMode();
            // Compute number of _TwoStagesSamplerRendererTask_s to create for rendering
            int nIterations = dualSampler->GetIterationCount();
            int nPixelsPerIteration = Ceil2Int(float(dualSampler->PixelsToSampleTotal()) / nIterations);
            nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
            nTasks = RoundUpPow2(nTasks);
            int nTasksTotal = nIterations * nTasks;
            // Set the progress reporter
            ProgressReporter reporterAdapt(nTasksTotal, "Adaptive Rendering");
            for (int iter = 0; iter < nIterations; iter++) {
                dualSampler->GetSamplingMaps(nPixelsPerIteration);

                // Generate tasks
                for (int i = 0; i < nTasks; ++i)
                    renderTasks.push_back(new TwoStagesSamplerRendererTask(
                        scene, this, camera, reporterAdapt, sampler, sample,
                        visualizeObjectIds, nTasks-1-i, nTasks, true));

                // Do the work
                EnqueueTasks(renderTasks);
                WaitForAllTasks();

                // Clean up
                for (uint32_t i = 0; i < renderTasks.size(); ++i)
                    delete renderTasks[i];
                renderTasks.clear();
            }
            dualSampler->Finalize();
            reporterAdapt.Done();
        }
    }
    else {
        // Compute number of _TwoStagesSamplerRendererTask_s to create for rendering
        int nPixels = camera->film->xResolution * camera->film->yResolution;
        int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
        nTasks = RoundUpPow2(nTasks);
        ProgressReporter reporter(nTasks, "Rendering");
        vector<Task *> renderTasks;
        for (int i = 0; i < nTasks; ++i)
            renderTasks.push_back(new TwoStagesSamplerRendererTask(scene, this,
                camera, reporter, sampler, sample, visualizeObjectIds,
                nTasks-1-i, nTasks));
        EnqueueTasks(renderTasks);
        WaitForAllTasks();
        for (uint32_t i = 0; i < renderTasks.size(); ++i)
            delete renderTasks[i];
        reporter.Done();
    }

    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();
}


Spectrum TwoStagesSamplerRenderer::Li(const Scene *scene,
    const RayDifferential &ray, const Sample *sample, RNG &rng,
    MemoryArena &arena, Intersection *isect,
    Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)) {
        Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena);
    }
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    // Restore the old isect, in case it was NULL
    return *T * Li + Lvi;
}


Spectrum TwoStagesSamplerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


