// A simple set of realistic camera lens classes
// Nolan Goodnight <ngoodnight@cs.virginia.edu>
#ifndef _LENSDATA_H_
#define _LENSDATA_H_

#include <vector>
#include <geometry.h>
#include "lensview.h"
#include "geometry.h"
#include "transform.h"

enum LENSMODE {
	PRECISE,
	THICK, 
};
#define 	RAY_EPSILON   1e-3f
class CRay {
public:
	// CRay Public Methods
	CRay(): mint(RAY_EPSILON), maxt(INFINITY), time(0.f) {}
	CRay(const Point &origin, const Vector &direction,
		float start = RAY_EPSILON, float end = INFINITY, float t = 0.f)
		: o(origin), d(direction), mint(start), maxt(end), time(t) { }
	Point operator()(float t) const { return o + d * t; }
	CRay& operator = ( const CRay& r) { 
		this -> d = r.d; this->o= r.o; this->maxt = r.maxt; this->mint = r.mint; this->time = r.time;
		return (*this); 
	}
	void print() { 
		printf("\no.x %f, o.y %f, o.z %f, d.x %f, d.y %f, d.z %f, mint %f, maxt %f, time %f",o.x, o.y, o.z, d.x, d.y, d.z, mint, maxt, time);
	}	
	// CRay Public Data
	Point o;
	Vector d;
	mutable float mint, maxt;
	float time;
};

class CRayChain {
public:
	// CRayChain Public Methods
	CRayChain() {};
	void draw();
	void print();
	//
	vector<CRay> rayChain; 
};


// A simple class for the camera lens interface elements
class lensElement {
public:
    lensElement(void): 
	    isStop(true), rad(0), zpos(0), ndr(0), aper(0), isActive(false) {
		    edge[X] = edge[Y] = 0.0f;
		}
	lensElement(float _rad, float _xpos, float _ndr, float _aper) {
		this->Init(_rad, _xpos, _ndr, _aper);
		edge[X] = edge[Y] = 0.0f;
	}
   
	// Method to initialize a lens interface object
	void Init(float _rad, float _xpos, float _ndr, float _aper) {
		this->isStop = false;
		rad = _rad; zpos = _xpos; 
		ndr = _ndr; aper = _aper;
		this->isActive = true;
		if (rad == 0.0f && ndr == 0.0f)
			this->isStop = true;
	}
	// Draw a single lens interfaces
	void Draw(void);
	bool isAir(void) { return (ndr == 1.0f)?true:false; }
	bool lensElement::refractOnce(CRay * outRay, const CRay * inRay, const float input_miu );
	bool lensElement::checkStop(CRay * outRay, const CRay * inRay);


    bool isStop;   // Flag to identify the aperture stop
    float rad;     // Radius for the lens interface
	float zpos;    // Lens interface position (1D);
	float ndr;     // Index of refraction
	float aper;    // Aperture (diameter) of the interface
	bool isActive; // Flag to set interface active
	float edge[2]; // Top edge point of each interface  //edge[X] = z edge[Y] = y, are used in drawing
};

// A simple class for the camera lens systems: lens interfaces
class lensSystem {
public:
	lensSystem(void): f1(0), f2(0), p1(0), p2(0), fstop(1) {
		this->oPlane = INFINITY; this->iPlane = f2; m_fleffective_aperture = -2;}
	lensSystem(const char *file): f1(0), f2(0), p1(0), p2(0), fstop(1) {
		this->oPlane = INFINITY; this->iPlane = f2;
		this->Load(file); m_fleffective_aperture = -2;
	}
	// Method to load a lens specification file and build
	// a vector of lens interfaces. 
	bool Load(const char *file);
	// Get the number of lens interfaces
	int numLenses(void) { return int(lenses.size()); }
	// Draw the entire lens system (all interfaces)
	void Draw(void);
	// Get the maximum aperture in the lens system
	float maxAperture(void);
	// Get the minimum aperture in the lens system
	float minAperture(void);
	// Get the aperture interface
	lensElement *getAperture(void);
	void lensSystem::printLenses(void) ;
	
	
	//Given all the setting, the OUTPUT the f1 and f2
	void lensSystem::calF1F2P1P2(void);	  	  
	//Given all the setting, the OUTPUT the p1 and p2
	float lensSystem::cross2RaysAtZ(const CRay * ray1, const CRay * ray2) ;
	
	void lensSystem::refractExact_ImgToObj(CRayChain * pRayChain) ;
	void lensSystem::refractExact_ObjToImg(CRayChain * pRayChain);

	float lensSystem::calAperture(float input_fstop) ;
	
	void lensSystem::setMaxFSMinFS(void);
	Point lensSystem::calExitPupil() ;
	void lensSystem::initlsystem() ;
	void lensSystem::printConfig(void);
	float lensSystem::caliPlane(float input_oPlane);
	float lensSystem::caloPlane(float input_iPlane); 
	void DrawRaysFromiPlane(float lviPlaneY, float lviPlaneZ,LENSMODE mode );
	void lensSystem::refractThick_ImgToObj(CRayChain * pRayChain); 

    vector<lensElement> lenses; // Vector of lens elements
	float f1, f2;               // Focal points for the lens
	float p1, p2;               // Locations of principle planes
	Point pupil;			    // Exit pupil for the lens system
	
	float fstop;		        // Current fstop for the camera
	float maxFS;                // Maximum fstop for the camera
	float minFS;                // Minimum fstop for the camera

	float iPlane;				// Image plane
	float oPlane;				// Object focal plane

private: 
	float m_fleffective_aperture;
};


#endif

