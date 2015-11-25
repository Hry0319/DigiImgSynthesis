
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

    // Assume the origin of the camera system is at the front lens
    RasterToCamera = Scale(-1,1,1) * // fix the issue of horizontal flipping
                   Translate(Vector(-1*Xres*ScaleRate/2, -1*Yres*ScaleRate/2, -1*(filmdistance+SumofThick) )) * Scale(ScaleRate, ScaleRate, 1.0f);
}



float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
	
	Point P_ras(sample.imageX, sample.imageY, 0);
	Point P_screen;
	//RasterToCamera(Pras, &Pcamera);
	//RasterToScreen(P_ras, &P_screen);
    RasterToCamera(P_ras, &P_screen);
	//P_screen.x *= -1;
	//P_screen.y *= -1;
	//P_screen.z  = -1*(SumofThick+scenecam.filmdistance);

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
    //*ray = Ray(P_screen, Normalize(Vector(FirstSamplePoint - P_screen)), 0.f, INFINITY);
    ray->d = Normalize(Vector(FirstSamplePoint - P_screen));
    ray->o = P_screen;
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);

    
	
    // Trace the ray through all the lens surfaces
    Point P_hit;
	double thickLen = -1*SumofThick;
	for ( int index = lens.size()-1; index >= 0; index--){  
        if (lens[index].N == 0) {
          // Find intersection (ray-disk)
          float t = (lens[index].lensCenterZ - ray->o.z) / ray->d.z;
          float x = ray->o.x + t * ray->d.x;
          float y = ray->o.y + t * ray->d.y;
          if (x * x + y * y > lens[index].aperture*lens[index].aperture) return 0.0f;
        } else {
          // Find the intersection point, a.k.a. P_hit
          {
            Vector co = ray->o - Point(0.0f, 0.0f, lens[index].lensCenterZ);
            float c = co.LengthSquared();
            float d = Dot(co, ray->d);
            float e = d * d - c + lens[index].radius*lens[index].radius;
            if (e < 0.0f) return 0.0f;
            float t = ((lens[index].radius > 0.0f) ? (-d + sqrtf(e)) : (-d - sqrtf(e)));
            P_hit = ray->o + t * ray->d;
            if (P_hit.x*P_hit.x + P_hit.y*P_hit.y > lens[index].aperture*lens[index].radius)
              return 0.0f;
          }
      
          // get normal vector on the intersection point
          Vector normal;
          if (lens[index].radius > 0.0f)
            normal = Normalize(Point(0.0f, 0.0f, lens[index].lensCenterZ) - P_hit);
          else
            normal = Normalize(P_hit - Point(0.0f, 0.0f, lens[index].lensCenterZ));
      
          Vector input = ray->d, output;
      
          // Heckber's method
          {
            float c1, c2;
            c1 = -Dot(input, normal);
            c2 = 1.0f - lens[index].n_ratio2 * (1.0f - c1 * c1);
            if (c2 < 0.0f) return 0.0f;
            c2 = sqrtf(c2);
            output = lens[index].n_ratio * input + (lens[index].n_ratio * c1 - c2) * normal;
          }
      
          // next!
          ray->o = P_hit;
          ray->d = Normalize(output);
        }
	}
    
    ray->mint = scenecam.hither;
    ray->maxt = (scenecam.yon - scenecam.hither) / ray->d.z;
    CameraToWorld(*ray, ray);
	//ray->d = Normalize(ray->d);

    
    // BONUS: varying pixel weight
    float weight = Dot(Normalize(P_hit-P_screen),Vector(0,0,1));
    weight = weight * weight / fabsf(SumofThick + scenecam.filmdistance);
    weight = weight * weight * (lens.back().radius*lens.back().radius * M_PI);
  
    return weight;
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
    while (!specfile.eof()) {
        specfile.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
        {
		    lens.resize(lens.size()+1);
		    Lens& len = lens[lens.size()-1];
            sscanf(line, "%f %f %f %f\n", &len.radius, &len.axpos, &len.N, &len.aperture);
			SumofThick+=len.axpos;
            len.lensZoffset-=len.axpos;
            len.lensCenterZ = len.lensZoffset+len.radius;
        }
    }

    //printf("Read in %Iu lens from %s\n", lens.size(), filename.c_str());

    
	for ( int index = lens.size()-1; index >= 0; index--)
    {
        if (lens[index].N==0)
            lens[index].n_ratio = 1.0f;
        else {
            Lens len2;
            if (index-1 >= 0)
                Lens len2 = lens[index-1];
            else
                len2.N = 0.0f;
            if (index-1 == 0 || len2.N == 0.0f)
                lens[index].n_ratio = lens[index].N;
            else
                lens[index].n_ratio = lens[index].N / len2.N;
        }
        lens[index].n_ratio2 = lens[index].n_ratio * lens[index].n_ratio;
    }
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
