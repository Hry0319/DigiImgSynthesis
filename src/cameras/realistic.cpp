
#include "stdafx.h"
#include "cameras/realistic.h"
#include <stdio.h>
#include <iostream> 
#include <fstream>
#include <algorithm>
//#include <string>
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
	//scenecam.specfile			= specfile;	
	//scenecam.filmdistance		= filmdistance;
	//scenecam.aperture_diameter	= aperture_diameter;
	//scenecam.filmdiag			= filmdiag;
	//scenecam.hither				= hither;		
	//scenecam.yon				= yon;		
	//scenecam.shutteropen		= sopen;
	//scenecam.shutterclose		= sclose;
    
   
    printf("Read pbrt ___________________________ %s\n", specfile.c_str());
	if (specfile.compare("") != 0)  {
        ParseLens(specfile);
    }
    

}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
	
	//Point Pras(sample.imageX, sample.imageY, 0);
	//Point Pcamera;
	//RasterToCamera(Pras, &Pcamera);


	return 0;
}

void RealisticCamera::ParseLens(const string& filename)  {
    ifstream specfile(filename.c_str());
    if (!specfile) {
        fprintf(stderr, "Cannot open file %s\n", filename.c_str());
        exit (-1);
    }

    char line[512];
    //int index;
    while (!specfile.eof()) {
        specfile.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
        {
		    lens.resize(lens.size()+1);
		    Lens& len = lens[lens.size()-1];
            sscanf(line, "%f %f %f %f\n", len.aperture, len.index_of_refraction, len.lens_radius, len.z_axis_intercept);
        }
    }

    //printf("Read in %zu lens from %s\n", lens.size(), filename.c_str());
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
