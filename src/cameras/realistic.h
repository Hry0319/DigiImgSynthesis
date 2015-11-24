
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "paramset.h"
#include "film.h"

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Public Methods
	struct Lens{
		// Realistic camera-specific parameters
		string	specfile;	
		float	filmdistance;
		float	aperture_diameter;
		float	filmdiag;
		// Extract common camera parameters from \use{ParamSet}
		float	hither;		
		float	yon;			
		float   shutteropen;	
		float	shutterclose;	
	} lens;
    float lens_radius;
    float z_axis_intercept;
    float index_of_refraction;
    float aperture;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H