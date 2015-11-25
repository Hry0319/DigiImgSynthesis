// A simple set of realistic camera lens classes
// Nolan Goodnight <ngoodnight@cs.virginia.edu>
#include <fstream>
#include "lensdata.h"
using namespace std;

// Method to load a lens specification file and build
// a vector of lens interfaces. 
bool lensSystem::Load(const char *file) {
	ifstream lensfile(file);
	if (!lensfile) return false;
     
	this->lenses.clear();
	// Iterate over the lens element specs and build
	// the lens system
	float dist = 0.0f, zpos = 0.0f;
    while (!lensfile.eof()) {
        lensElement lens; 
		string data;

		lensfile >> data;
		// Remove comments and empty lines
		if (!data.compare("#") || !data.compare("")) {
		    getline(lensfile, data);
			continue;
		}
		float rad  = (float)atof(data.c_str());
		lensfile >> data;
		float axpos = dist - zpos; 
		dist -= zpos;
		zpos = (float)atof(data.c_str());
		
		lensfile >> data;
		float ndr  = (float)atof(data.c_str());
		lensfile >> data;
		float aper = (float)atof(data.c_str());
		if (!aper && !rad && !zpos && !ndr) 
			return false;

		// Initialize the lens interface and store
		lens.Init(rad, axpos, ndr, aper);		
		this->lenses.push_back(lens);
	}
	if (this->numLenses() < 5) 
		return false;
	return true;
}

// Draw the entire lens system (all interfaces and aperture)
void lensSystem::Draw(void) {

	glColor3f(0.75f, 0.75f, 0.75f); glLineWidth(1);
	vector<lensElement>::iterator Li, Ln;
	
	for (Li = lenses.begin(); Li != lenses.end(); Li++) 
		if (Li->isActive && !Li->isStop) 
			Li->Draw();

	for (Li = lenses.begin(); Li != lenses.end() - 1; Li++) {
		Ln = Li + 1; 
		if ((Ln->isStop || Li->isStop) || Li->isAir() || !Li->isActive) 
			continue;

		float lip[2] = {Li->edge[X], Li->edge[Y]};
		float lnp[2] = {Ln->edge[X], Ln->edge[Y]};

		if (Li->aper < Ln->aper) {
			drawLine(lip[X],  lip[Y], lip[X],  lnp[Y]);
	    	drawLine(lnp[X],  lnp[Y], lip[X],  lnp[Y]);	
			drawLine(lip[X], -lip[Y], lip[X], -lnp[Y]);
			drawLine(lnp[X], -lnp[Y], lip[X], -lnp[Y]);	
		} else {
			drawLine(lip[X],  lip[Y], lnp[X],  lip[Y]);
			drawLine(lnp[X],  lnp[Y], lnp[X],  lip[Y]);
			drawLine(lip[X], -lip[Y], lnp[X], -lip[Y]);
			drawLine(lnp[X], -lnp[Y], lnp[X], -lip[Y]);
		}
	}

	// Draw the aperture stop interface
	lensElement *stop = this->getAperture(); 
	float maxp = this->maxAperture();
	float haxp = stop->aper / 2.0f;
	float taxp = haxp + maxp / 2.0f;
	glColor3f(1, 1, 1);
	glLineWidth(2);
	drawLine(stop->zpos, haxp, stop->zpos, taxp);
	drawLine(stop->zpos, -haxp, stop->zpos, -taxp);
	glLineWidth(1);
}

// Get the maximum aperture in the lens system
float lensSystem::maxAperture(void) {
	float dist = 0.0f;
	vector<lensElement>::iterator Li;
	for (Li = lenses.begin(); Li != lenses.end(); Li++) 
		if (Li->aper > dist) 
			dist = Li->aper;
	return dist;
}

// Get the maximum aperture in the lens system
float lensSystem::minAperture(void) {
	float dist = INFINITY;
	vector<lensElement>::iterator Li;
	for (Li = lenses.begin(); Li != lenses.end(); Li++) 
		if (Li->aper < dist) 
			dist = Li->aper;
	return dist;
}

// Method to get the aperture interface
lensElement* lensSystem::getAperture(void) {
	int i;
	for (i = 0; i < numLenses(); i++) 
		if (lenses[i].isStop)
			return (&lenses[i]);
	return NULL;
}

// Draw a single lens interfaces
void lensElement::Draw(void) {
	float nextX, nextY, thisX, thisY;
	thisY = -this->aper / 2.0f;
	float rs = this->rad / fabs(this->rad);
	int i;
	
	float step = this->aper / LVSD;		
	float posx = sqrt(rad*rad - thisY*thisY);
	thisX = rs * posx - rad + zpos;
	
	for (i = 1; i <= LVSD; i++) {
	  nextY = thisY + step;
	  posx = sqrt(rad*rad - nextY*nextY);
	  nextX = rs * posx - rad + zpos;
	  drawLine(thisX, thisY, nextX, nextY);
	  thisX = nextX; thisY = nextY;
	}
	this->edge[X] = thisX;
	this->edge[Y] = thisY;
}
	  	  
