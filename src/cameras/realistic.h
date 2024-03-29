
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "paramset.h"
#include "film.h"
#include <vector>

#define IN
#define OUT


struct SceneCamera {
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
};


// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	struct Lens {
	public:
		float radius;
		float thickness;
		float axpos;
		float N;
		float aperture;
	};
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	~RealisticCamera();
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
    //Transform       RasterToCamera;
    SceneCamera     scenecam;
    vector<Lens>    lens;
	float			SumofThick;
	float			Xres;
	float			Yres;
	float			ScaleRate;
    float			filmpos;
	float			RasterDiag;

	// RealisticCamera Public Methods
    void  ParseLens(const string& filename);
	void  RasterToScreen(IN const Point Praster, OUT Point *P_camera) const;
	float RasterToCamera(float in, int dim) const;
    bool  SnellsLaw(Vector s1, Vector N, float n1, float n2, Vector *s2) const;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H