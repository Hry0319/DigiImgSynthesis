
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
    
	
    // Trace the ray through all the lens surfaces
    Point P_hit;
	double thickLen = -1*SumofThick;
	for ( int index = lens.size()-1; index >= 0; index--){    
        Ray r = *ray;
        thickLen+=lens[index].axpos;
        //
        // intersection on a sphere
        //
        Vector IntersectNormal;
        if (lens[index].lens_radius != 0) {
		    float thit = 0;	
		    Vector w2oV(0.f, 0.f, lens[index].lens_radius - thickLen);
		    Transform   o2w = Translate(-1*w2oV), 
			            w2o = Translate( 1*w2oV);
            
		    Sphere sphere(&o2w, &w2o, false, fabs(lens[index].lens_radius), -1*lens[index].lens_radius, lens[index].lens_radius, 360);

		    float rayEpsilon;
		    DifferentialGeometry dg;
            CameraToWorld(r, &r);
		    if (!sphere.Intersect(r, &thit, &rayEpsilon, &dg)) return 0.f;
		    //if (thit > r.maxt || thit < r.mint) return 0.f;
            //w2o(r, &r);
		    P_hit = (*ray)(thit);
	    }
        else
        {            
		    // using absolute value because ray can come from either side
		    float scale = fabs((thickLen - r.o.z)/r.d.z);
		    P_hit = r.o + scale*r.d; 
        }       
	    // aperture
	    if ( (P_hit.x * P_hit.x + P_hit.y * P_hit.y) >= (lens[index].aperture*lens[index].aperture)/4) 
        {
            return 0.f;
        }
        IntersectNormal = Normalize(P_hit - Point(0, 0, thickLen - lens[index].lens_radius));
        

        // assert the aperture is in air
		if (lens[index].lens_radius == 0) {
			assert(lens[index].N == 1.f);
			assert(lens[index-1].N == 1.f);
		}


        // update ray direction
		Vector newD;
		float n2 = (index==0)? 1 : lens[index-1].N;
		// flipping normal direction if radius is positive
		// so that the normal will always be on the same side as the incident ray
		if (lens[index].lens_radius * ray->d.z > 0) 
			IntersectNormal *= -1;
		
		if (!SnellsLaw(ray->d, IntersectNormal, lens[index].N, n2, &newD)) {
			//cout << "total internal reflection" << endl;
			return 0.f;
		} 
		
		*ray = Ray(P_hit, Normalize(newD), 0.f, INFINITY, ray->time);
	}
    
	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);




 //   // GenerateRay() should return the weight of the generated ray
	//// E = A*cos^4(theta)/Z^2 
	//double A = M_PI* pow(lens.back().aperture/2, 2);
 //   double Z = lens.back().axpos;
 //   double costheta = Dot(ray->d, Vector(0,0,1));    
 //   double E = A*pow(costheta,4) /(Z*Z);
    
    // BONUS: varying pixel weight
    //float weight = Dot(Normalize(ray->d),Vector(0,0,1));
    //weight = weight * weight / fabsf(P_screen.z);
    //weight = weight * weight * (lens.back().aperture*lens.back().aperture * M_PI);
  
    return 1.f;
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

void RealisticCamera::RasterToScreen(IN const Point Praster, OUT Point *P_screen) const {
	P_screen->x = (Praster.x  - (Xres/2) ) * ScaleRate;
	P_screen->y = (Praster.y  - (Yres/2) ) * ScaleRate;
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
