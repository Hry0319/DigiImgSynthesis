
#include "stdafx.h"
#include "cameras/realistic.h"
#include "math.h"
#include <stdio.h>
#include <iostream> 
#include <fstream>
#include "shapes/sphere.h"

using namespace std; 


RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
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
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	Point P_ras(sample.imageX, sample.imageY, 0);
	Point P_screen;
	RasterToScreen(P_ras, &P_screen);
	P_screen.x *= -1;
	P_screen.z  = -1*(SumofThick+scenecam.filmdistance);

	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);	
    lensU *= lens.back().aperture/2;
	lensV *= lens.back().aperture/2;    

    double _d = sqrt( lens.back().radius*lens.back().radius - sqrt( lensU*lensU + lensV*lensV ) );
    float initRayZ = 0;
    if (lens.back().radius > 0)
    {
        initRayZ = -1*SumofThick - lens.back().radius - _d;
    } else {
        initRayZ = -1*SumofThick + lens.back().radius - _d;
    }
	Point FirstSamplePoint (lensU, lensV, initRayZ);
    *ray = Ray(P_screen, Normalize(Vector(FirstSamplePoint - P_screen)), 0.f, INFINITY);
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	//Ray r = *ray;
	
	Point p_hit;
	//double thickLen = -1*SumofThick;
	for ( int index = lens.size()-1; index >= 0; index--){
		
		if(lens[index].N == 0 || lens[index].radius == 0)
		{
			float scale = fabs((lens[index].axpos - ray->o.z)/ray->d.z);
			p_hit = ray->o + scale*ray->d; 
		}
		else
		{
			float tHit = 0;
			Vector		w2oV(0.f, 0.f, lens[index].radius - lens[index].axpos);
			Transform	o2w = Translate(-1*w2oV), 
						w2o = Translate(1*w2oV);

			Sphere sphere(&o2w, &w2o, false, fabs(lens[index].radius), -1*lens[index].radius, lens[index].radius, 360);

			float rayEpsilon;
			DifferentialGeometry dg;
			if (!sphere.Intersect(*ray, &tHit, &rayEpsilon, &dg)) return 0.f;
			if (tHit > ray->maxt || tHit < ray->mint) return 0.f;

			p_hit = (*ray)(tHit);
		}		
		if ( (p_hit.x * p_hit.x + p_hit.y * p_hit.y) >= (lens[index].aperture*lens[index].aperture)/4) 
			return 0.f; 

		Vector Normal;
		// returned normal always points from the sphere surface outward
		Normal = Normalize(p_hit - Point(0, 0, lens[index].axpos - lens[index].radius));

		/*if (lens[index].radius == 0) {
			assert(lens[indexi].N == 1.f);
			assert(lens[index-1].N == 1.f);
		}*/

		
		// update ray direction
		Vector newD;
		//float n1;
		float n2 = (index==0)? 1 : lens[index-1].N;
		if (lens[index].radius * ray->d.z > 0) 
			Normal = Normal*-1;
		
		if (!SnellsLaw(ray->d, Normal, lens[index].N, n2, &newD)) {
			//cout << "total internal reflection" << endl;
			return 0.f;
		} 
		*ray = Ray(p_hit, Normalize(newD), 0.f, INFINITY, ray->time);
	}
	
	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);

	return 1.f;
}

void RealisticCamera::RasterToScreen(IN const Point Praster, OUT Point *P_screen) const {
	P_screen->x = (Praster.x  - (Xres/2) ) * ScaleRate;
	P_screen->y = (Praster.y  - (Yres/2) ) * ScaleRate;

	return;
}

bool RealisticCamera::SnellsLaw(Vector s1, Vector N, float n1, float n2, Vector *s2) const {

	float mu = n1/n2;
	
	float costheta = Dot(-1*s1, N); // note normal is pointing in opposite direction of s1
	float toSqrt = 1 - mu*mu* ( 1- costheta*costheta);
	
	if (toSqrt < 0) return false; // total internal reflection

	int sign = 1;
	if (costheta < 0) sign = -1;

	float gamma = mu*costheta - sign*sqrtf(toSqrt); 

	*s2 = mu*s1 + gamma*N;
	return true;
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
            sscanf(line, "%f %f %f %f\n", &len.radius, &len.thickness, &len.N, &len.aperture);
			SumofThick+=len.thickness;
			len.axpos = -1*SumofThick;
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
