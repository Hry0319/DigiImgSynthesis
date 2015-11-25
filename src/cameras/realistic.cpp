
#include "stdafx.h"
#include "cameras/realistic.h"
#include "math.h"
#include <stdio.h>
#include <iostream> 
#include <fstream>

using namespace std; 


RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
   // YOUR CODE HERE -- build and store datastructures representing the given lens
   // and film placement.
	scenecam.specfile			= specfile;	
	scenecam.filmdistance		= filmdistance;
	scenecam.aperture_diameter	= aperture_diameter;
	scenecam.filmdiag			= filmdiag;
	scenecam.hither				= hither;		
	scenecam.yon				= yon;		
	scenecam.shutteropen		= sopen;
	scenecam.shutterclose		= sclose;
	
	Xres				= f->xResolution;
	Yres				= f->yResolution;
	double RasterDiag	= sqrt( Xres*Xres + Yres*Yres );
	ScaleRate			= abs(filmdiag/RasterDiag);
	SumofThick          = 0;
    
   //current_path();
	if (specfile.compare("") != 0)  {        
        ParseLens("C:\\RENDERING\\HW3\\project2_template\\"+specfile);
    }
 //   for (int i=0;i<lens.size();i++){
 //       printf("%d : [ %5.5f %5.5f %5.5f %5.5f ] \n", i, lens[i].lens_radius, lens[i].axpos, lens[i].N, lens[i].aperture);
 //   }
	//printf("SumofThick = %f \n", SumofThick);    

    // Assume the origin of the camera system is at the front lens
    RasterToCamera = Scale(-1,1,1) * // fix the issue of horizontal flipping
                   Translate(Vector(-1*Xres*ScaleRate/2, -1*Yres*ScaleRate/2, -1*filmdistance)) * 
                   Scale(ScaleRate, ScaleRate, 1.0f);
}



float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
	
	Point P_ras(sample.imageX, sample.imageY, 0);
	Point P_screen;
	//RasterToCamera(Pras, &Pcamera);
	RasterToScreen(P_ras, &P_screen);
	P_screen.x *= -1;
	//P_screen.y *= -1;
	P_screen.z  = -1*(SumofThick+scenecam.filmdistance);

	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);	
    lensU *= lens.back().aperture/2;
	lensV *= lens.back().aperture/2;    

    double _d = sqrt( lens.back().lens_radius*lens.back().lens_radius - sqrt( lensU*lensU + lensV*lensV ) );
    float initRayZ = 0;
    if (lens.back().lens_radius > 0)
    {
        initRayZ = -1*SumofThick - lens.back().lens_radius - _d;
    } else {
        initRayZ = -1*SumofThick + lens.back().lens_radius - _d;
    }
	Point FirstSamplePoint (lensU, lensV, initRayZ);
    *ray = Ray(P_screen, Normalize(Vector(FirstSamplePoint - P_screen)), 0.f, INFINITY);
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);

	
	//double thickLen = -1*SumofThick;
	//for ( int index = lens.size()-1; index >= 0; index--){
	//	//// Modify ray for depth of field
	//	if (lens[index].lens_radius > 0.) {
	//		// Sample point on lens
	//		float lensU, lensV;
	//		ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	//		lensU *= (lens[index].aperture/2);
	//		lensV *= (lens[index].aperture/2);

	//		double _d = sqrt( lens[index].lens_radius*lens[index].lens_radius - sqrt( lensU*lensU + lensV*lensV ) );
	//		Point Rayway(lensU, lensV, thickLen-lens[index].lens_radius+_d);

	//		// Compute point on plane of focus
	//		//float ft = focaldistance / ray->d.z;
	//		//Point Pfocus = (*ray)(ft);

	//		// Update ray for effect of lens
	//		ray->o = Point(lensU, lensV, 0.f);
	//		ray->d = Normalize(Rayway - ray->o);
	//		//ray->o = Point(lensU, lensV, 0.f);
	//		
	//		thickLen += lens[index].axpos;
	//	}
	//}

    CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);

	return 1.f;
}

void RealisticCamera::RasterToScreen(IN const Point Praster, OUT Point *P_screen) const {
	P_screen->x = (Praster.x  - (Xres/2) ) * ScaleRate;
	P_screen->y = (Praster.y  - (Yres/2) ) * ScaleRate;

	return;
}

void RealisticCamera::ParseLens(const string& filename)  {
    ifstream specfile(filename.c_str());
    if (!specfile) {
        fprintf(stderr, "Cannot open file %s\n", filename.c_str());
        exit (-1);
    }

    char    line[512];
    int     index = 0;
    while (!specfile.eof()) {
        specfile.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
        {
		    lens.resize(lens.size()+1);
		    Lens& len = lens[lens.size()-1];
            sscanf(line, "%f %f %f %f\n", &len.lens_radius, &len.axpos, &len.N, &len.aperture);
			SumofThick+=len.axpos;
        }
    }

    printf("Read in %Iu lens from %s\n", lens.size(), filename.c_str());
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {

	// Realistic camera-specific parameters
	string specfile		= params.FindOneString("specfile", "");
	float filmdistance	= params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop			= params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag		= params.FindOneFloat("filmdiag", 35.0);

	// Extract common camera parameters from \use{ParamSet}
	float hither		= params.FindOneFloat("hither", -1);
	float yon			= params.FindOneFloat("yon", -1);
	float shutteropen	= params.FindOneFloat("shutteropen", -1);
	float shutterclose	= params.FindOneFloat("shutterclose", -1);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
