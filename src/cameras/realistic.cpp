
#include "stdafx.h"
#include "cameras/realistic.h"

#include <iostream> 
#include <fstream>
#include "shapes/sphere.h"
#include "intersection.h"

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

	SumofThick					= 0;

	if (specfile.compare("") != 0)  {        
        ParseLens("C:\\RENDERING\\HW3\\project2_template\\"+specfile);
    } 
	lens[lens.size()-1].thickness = filmdistance;

	Xres				= f->xResolution;
	Yres				= f->yResolution;
	RasterDiag			= sqrtf( Xres*Xres + Yres*Yres );
	ScaleRate			= abs(filmdiag/RasterDiag);
    filmpos             = lens.back().axpos - lens.back().thickness;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	//Point P_ras(sample.imageX, sample.imageY, filmpos);
	Point P_ras(
		RasterToCamera(sample.imageX, 0), 
		RasterToCamera(sample.imageY, 1), 
		filmpos);
	/*
	Point P_camera;
	RasterToScreen(P_ras, &P_camera);*/
	
	float _lensU, _lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &_lensU, &_lensV);	
    _lensU *= lens.back().aperture/2;
	_lensV *= lens.back().aperture/2;    

    /*
	double _d = sqrtf( lens.back().radius*lens.back().radius - sqrtf( _lensU*_lensU + _lensV*_lensV ) );
    float initRayZ = 0;
    if (lens.back().radius > 0)
    {
        initRayZ = -1*SumofThick - lens.back().radius - _d;
    } else {
        initRayZ = -1*SumofThick + lens.back().radius - _d;
    }
	*/

	Point FirstSamplePoint (_lensU, _lensV, /*initRayZ*/lens.back().axpos);
	Ray r(P_ras, Normalize(FirstSamplePoint - P_ras), 0.f, INFINITY);
	
	float A = M_PI* pow(lens.back().aperture/2., 2);
    float Z = lens.back().thickness;
    float costheta = Dot(r.d, Vector(0,0,1));    
    float E = A*pow(costheta,4) /(Z*Z);

	for ( int index = lens.size()-1; index >= 0; index--){		
		Point p_hit;
		Vector Normal;
		
		if(lens[index].radius == 0)
		{
			float scale = fabs((lens[index].axpos - r.o.z)/r.d.z);
			p_hit = r.o + scale*r.d; 
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
			if (!sphere.Intersect(r, &tHit, &rayEpsilon, &dg)) return 0.f;
			if (tHit > r.maxt || tHit < r.mint) return 0.f;

			p_hit = r(tHit);
		}		
		if ( (p_hit.x * p_hit.x + p_hit.y * p_hit.y) >= (lens[index].aperture*lens[index].aperture)/4) 
			return 0.f; 

		// returned normal always points from the sphere surface outward
		Normal = Normalize(p_hit - Point(0, 0, lens[index].axpos - lens[index].radius));

		
		// update ray direction
		Vector newD;
		float n2 = (index==0)? 1 : lens[index-1].N;
		if (lens[index].radius * r.d.z > 0) 
			Normal = Normal*-1;
		
		if (!SnellsLaw(r.d, Normal, lens[index].N, n2, &newD)) {
			return 0.f;
		} 
		r = Ray(p_hit, Normalize(newD), 0.f, INFINITY);
	}
	
	r.time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(r, ray);
	ray->d = Normalize(ray->d);

	return E;
}

//void RealisticCamera::RasterToScreen(IN const Point Praster, OUT Point *P_camera) const {
//	P_camera->x = ( (Xres*0.5) - Praster.x ) * ScaleRate;
//	P_camera->y = ( (Yres*0.5) - Praster.y ) * ScaleRate;
//
//	//P_camera->x *= -1;
//	////P_camera->z  = -1*(SumofThick+scenecam.filmdistance);
//
//
//	return;
//}
float RealisticCamera::RasterToCamera(float in, int dim) const {

	if (dim == 0) {
		float res = film->xResolution;
		return (res/2.f - in) * scenecam.filmdiag/RasterDiag;
	} else {
		float res = film->yResolution;
		return (in - res/2.f) * scenecam.filmdiag/RasterDiag;
	}
	
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

    while (!specfile.eof()) {
        specfile.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
        {
		    //lens.resize(lens.size()+1);
		    Lens len;
            sscanf(line, "%f %f %f %f\n", &len.radius, &len.thickness, &len.N, &len.aperture);

			if (lens.size() == 0) 
				len.axpos = 0;
			else 
				len.axpos = lens.back().axpos - lens.back().thickness; // note this is negative

			// real aperture diameter is given by parameter rather than by data file
			// (apeture has radius of 0)
			if (len.radius == 0) {
				len.aperture = scenecam.aperture_diameter; 
				len.N = 1.f;
			}

			//len.axpos = -1*SumofThick;
			//SumofThick+=len.thickness;
			lens.push_back(len);
        }
    }
    //printf("Read in %Iu lens from %s\n", lens.size(), filename.c_str());
}


RealisticCamera::~RealisticCamera()
{

}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {	
	// Extract common camera parameters from \use{ParamSet}
	float hither		= params.FindOneFloat("hither", -1);
	float yon			= params.FindOneFloat("yon", -1);
	float shutteropen	= params.FindOneFloat("shutteropen", -1);
	float shutterclose	= params.FindOneFloat("shutterclose", -1);
	// Realistic camera-specific parameters
	string specfile		= params.FindOneString("specfile", "");
	float filmdistance	= params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop			= params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag		= params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
