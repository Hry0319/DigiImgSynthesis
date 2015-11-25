// A simple set of realistic camera lens classes
// Nolan Goodnight <ngoodnight@cs.virginia.edu>
#ifndef _LENSDATA_H_
#define _LENSDATA_H_

#include <vector>
#include <geometry.h>
#include "lensview.h"

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

    bool isStop;   // Flag to identify the aperture stop
    float rad;     // Radius for the lens interface
	float zpos;    // Lens interface position (1D);
	float ndr;     // Index of refraction
	float aper;    // Aperture (diameter) of the interface
	bool isActive; // Flag to set interface active
	float edge[2]; // Top edge point of each interface
};

// A simple class for the camera lens systems: lens interfaces
class lensSystem {
public:
	lensSystem(void): f1(0), f2(0), p1(0), p2(0), fstop(1) {
		this->oPlane = INFINITY; this->iPlane = f2; }
	lensSystem(const char *file): f1(0), f2(0), p1(0), p2(0), fstop(1) {
		this->oPlane = INFINITY; this->iPlane = f2;
		this->Load(file);
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

    vector<lensElement> lenses; // Vector of lens elements
	float f1, f2;               // Focal points for the lens
	float p1, p2;               // Locations of principle planes
	Point pupil;			    // Exit pupil for the lens system
	
	float fstop;		        // Current fstop for the camera
	float maxFS;                // Maximum fstop for the camera
	float minFS;                // Minimum fstop for the camera

	float iPlane;				// Image plane
	float oPlane;				// Object focal plane
};

#endif

