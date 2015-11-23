// A simple set of realistic camera lens classes
// Nolan Goodnight <ngoodnight@cs.virginia.edu>
// Define lensSystem  lensElement
#include <fstream>
#include "lensdata.h"
using namespace std;

#define PLANE_HEIGHT  30 
#define SQUARE_WIDTH  4
#define OFFSET  10
// Method to load a lens specification file and build
// a vector of lens interfaces. 
bool lensSystem::Load(const char *file) {
	ifstream lensfile(file);
	if (!lensfile) return false;
     
	this->lenses.clear();
	// Iterate over the lens element specs and build
	// the lens system
	float dist = 0.0f;
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
		float zpos = (float)atof(data.c_str());
		float axpos = dist - zpos; 
		lensfile >> data;
		float ndr  = (float)atof(data.c_str());
		lensfile >> data;
		float aper = (float)atof(data.c_str());
		if (!aper && !rad && !zpos && !ndr) 
			return false;

		// Initialize the lens interface and store
		lens.Init(rad, axpos, ndr, aper);		
		this->lenses.push_back(lens);
		dist -= zpos;
	}
	if (this->numLenses() < 5) 
		return false;
	return true;
}

// Draw the entire lens system (all interfaces and aperture)
void lensSystem::Draw(void) {

	glColor3f(0.75f, 0.75f, 0.75f); glLineWidth(1);
	vector<lensElement>::iterator Li, Ln;
	
	//The curve for the Lens
	for (Li = lenses.begin(); Li != lenses.end(); Li++) 
		if (Li->isActive && !Li->isStop) 
			Li->Draw();
	//The Horizontal lines for the Lens
	for (Li = lenses.begin(); Li != lenses.end() - 1; Li++) {
		Ln = Li + 1; //!!!
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
	if (m_fleffective_aperture < 0) {
		drawLine(stop->zpos, haxp, stop->zpos, taxp);
		drawLine(stop->zpos, -haxp, stop->zpos, -taxp);
	} else {
		haxp = m_fleffective_aperture / 2.0f; 
		drawLine(stop->zpos, haxp, stop->zpos, taxp);
		drawLine(stop->zpos, -haxp, stop->zpos, -taxp);
	}
	glLineWidth(1);
	
	//Draw the exit pupil
	Point ep = this ->calExitPupil();
	float hep = ep.y; 
	float tep = hep + maxp/2.0f; 
	float zep = ep.z; 
	glColor3f(0, 1, 0);
	glLineWidth(2);
    drawLine(zep, hep, zep, tep);
	drawLine(zep, -hep, zep, -tep);			
	glLineWidth(1);

	//a. P and P' 
	float p1z = this->p1; 
	float p2z = this->p2;
	glColor3f(1, 0, 0);
	glLineWidth(2);
	drawLine(p1z, -PLANE_HEIGHT , p1z, PLANE_HEIGHT);	
	drawLine(p2z, -PLANE_HEIGHT , p2z, PLANE_HEIGHT);	
	drawText(p1z , -PLANE_HEIGHT - OFFSET, SBFONT, "P"); 
	drawText(p2z , -PLANE_HEIGHT - OFFSET,  SBFONT, "P'"); 

	glLineWidth(1);
	
	//b. F and F'    //by lsystem.draw()
	float f1z = this->f1; 
	float f2z = this->f2; 	
	drawSquare(1,1,1, f1z - SQUARE_WIDTH/2, - SQUARE_WIDTH/2  , f1z + SQUARE_WIDTH/2 ,  SQUARE_WIDTH/2);	
	drawSquare(1,1,1, f2z - SQUARE_WIDTH/2, - SQUARE_WIDTH/2  , f2z + SQUARE_WIDTH/2 ,  SQUARE_WIDTH/2);	
	drawText(f1z , - OFFSET,  SBFONT, "F"); 
	drawText(f2z , - OFFSET,  SBFONT, "F'"); 
	
	//c. I and O
	float iPlanez = this->iPlane; 
	float oPlanez = this->oPlane;
	glColor3f(1, 0, 1);
	glLineWidth(2);
	drawLine(iPlanez, -PLANE_HEIGHT , iPlanez, PLANE_HEIGHT);	
	drawLine(oPlanez, -PLANE_HEIGHT , oPlanez, PLANE_HEIGHT);	
	drawText(iPlanez, -PLANE_HEIGHT - OFFSET , SBFONT, "I");	
	drawText(oPlanez, -PLANE_HEIGHT - OFFSET , SBFONT, "O");		
	glLineWidth(1);
}

void lensSystem::DrawRaysFromiPlane(float iPlaneY, float iPlaneZ, LENSMODE mode ){
	Vector upperDirection = Normalize( Vector(0, pupil.y - iPlaneY, pupil.z - iPlaneZ) );
	Vector lowerDirection = Normalize( Vector(0, -pupil.y - iPlaneY, pupil.z - iPlaneZ) );
	int numRays = 20; 
	for (int i = 0 ; i <= numRays; i++) {
		CRayChain raychains; 
		float w1 = i; 
		float w2 = numRays - i ;
		Vector curDirection = Normalize( upperDirection * w1 + lowerDirection * w2 );
		CRay firstRay(Point(0, iPlaneY, iPlaneZ), curDirection);
		raychains.rayChain.push_back(firstRay);
		if (mode == PRECISE)
			refractExact_ImgToObj(&raychains); 	
		else 
			refractThick_ImgToObj(&raychains); 
		raychains.draw();
	}
	
}

void lensSystem::printConfig(void) {		
	printf("f1 %1.2f \t f2 %1.2f \t fstop %1.2f \t iPlane %1.2f \t	m_fleffective_aperture %1.2f \t maxFS %1.2f \t minFS %1.2f \t oPlane %1.2f \t p1 %1.2f \t p2 %1.2f \t pupil(%1.2f %1.2f %1.2f) \t \n",
		f1, 
		f2, 
		fstop, 
		iPlane,
		m_fleffective_aperture,
		maxFS, 
		minFS,
		oPlane,
		p1, 
		p2,
		pupil.x, pupil.y, pupil.z 
	);	
}

void lensSystem::printLenses(void) {
	vector<lensElement>::iterator Li;
	
	for (Li = lenses.begin(); Li != lenses.end(); Li++) {
		
			printf( "aper %f;\t edge %f,%f;\t isActive %d; \t isAir() %d; \t isStop %d;\t ndr %f;\t rad %f; \t zpos %f;\t \n", 
			Li->aper,
			Li->edge[0],Li->edge[1],		//edge[X] = z  , edge[Y] = y, are used in drawing

			Li->isActive,					//???
			Li->isAir(),

			Li->isStop,
			Li->ndr,
			Li->rad,
			Li->zpos );
		
	}
	//The data for the lensSystem:
}
// Get the maximum aperture in the lens system; mainly for display()
float lensSystem::maxAperture(void) {
	float dist = 0.0f;
	vector<lensElement>::iterator Li;
	for (Li = lenses.begin(); Li != lenses.end(); Li++) 
		if (Li->aper > dist) 
			dist = Li->aper;
	return dist;
}

// Get the maximum aperture in the lens system;  mainly for display()
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
	this->edge[X] = thisX;   //edge[X] = z edge[Y] = y, are used in drawing
	this->edge[Y] = thisY;
}

		
inline bool myQuadratic(float A, float B, float C, float *t0,
		float *t1) {
	// Find quadratic discriminant
	float discrim = B * B - 4.f * A * C;
	if (discrim < 0.) return false;
	float rootDiscrim = sqrtf(discrim);
	// Compute quadratic _t_ values
	float q1, q2;
	//if (B < 0) q = -.5f * (B - rootDiscrim);
	//else       q = -.5f * (B + rootDiscrim);
	q1 = .5f * (-B - rootDiscrim);
	q2 = .5f * (-B + rootDiscrim);
	*t0 = q1 / A;
	*t1 = q2 / A;
	if (*t0 > *t1) swap(*t0, *t1);
	return true;
}		
				
inline void invalidateRay(CRay * outRay) {
	outRay->maxt = 0;
	outRay->mint = 1;
	return;
}

//Given an input ray, calculate the intersection, and do the refraction once, then OUTPUT the out ray. 
//OUTPUT: (*outRay).o and (*outRay).d;
//		   if the ray cannot penetrate, invalidate the ray, and return false
bool lensElement::refractOnce(CRay * outRay, const CRay * inRay, const float input_miu ){
	
	//The inRay intersect with the interface (SPHERE), then REFRACT to produce the outRay. 	
	Point c (0,0,(zpos - rad)); //center of the sphere
	float r = rad;				//radius of the sphere
	Point o = inRay->o;			//Input Ray. origin
	Vector d = inRay->d;		//Input Ray. direction
	float A = d.Length() ; 
	float B = 2* Dot((o-c), d); 
	float C = Dot((o-c),(o-c)) - r*r;	
	
	float t0, t1;
	if (!myQuadratic(A, B, C, &t0, &t1)) {
		//set the ray to be invalid	
		invalidateRay(outRay);
		return false ;
	}
	if (t0 >  (*inRay).maxt || t1 <  (*inRay).mint){
		//set the ray to be invalid		
		invalidateRay(outRay);				
		return false;
	}
	//All we set is outRay->o , ->d; the other parts we use inRay
	(*outRay) = (*inRay);
	//Now we already got t0 and t1. So which one is the intersection depends on which one is close to the lense itself. 
	Point point_inRay_t0 = (*inRay)(t0); 
	Point point_inRay_t1 = (*inRay)(t1); 
	Point lens_pos (0, 0, zpos);
	float dist_t0 = (point_inRay_t0 - lens_pos).Length();
	float dist_t1 = (point_inRay_t1 - lens_pos).Length();
	if (dist_t0 < dist_t1 ) {
		(*outRay).o = point_inRay_t0;
	} else {
		(*outRay).o = point_inRay_t1; 
	}
	//If the intersection point is outside of the lense itself, then return false;
	float intersectionPy = (*outRay).o.y;
	if ( abs(intersectionPy) > (aper/2) ) {
		invalidateRay(outRay);
		return false;
	}



	//Here we should cal the refraction now	
	Point intersectionPoint = (*outRay).o;
	Vector I = Normalize( (*inRay).d );
	Vector N = Normalize( intersectionPoint - c ) ;
	float miu = input_miu;
	float DotIN = Dot(I, N);
	float gamma_first = -miu * DotIN ;
	float Delta = 1 - miu* miu * (1- DotIN*DotIN);
	if (Delta < 0) {
		//set the ray to be invalid		
		invalidateRay(outRay);				
		return false ;
	}
	float gamma_second = sqrtf( Delta ) ; 
	float gamma01 = gamma_first - gamma_second;
	float gamma02 = gamma_first + gamma_second;
	Vector possibleT01 = miu * I + gamma01 * N;
	Vector possibleT02 = miu * I + gamma02 * N;
	if (Dot(possibleT01, I) >= 0 ) 
		(*outRay).d =Normalize( possibleT01);
	else {
		assert( Dot(possibleT02, I) >= 0 );
		(*outRay).d =Normalize( possibleT02);
	}
	//Here do one more assert
#ifdef _DEBUG 
	Vector T = (*outRay).d ; 
	Vector NxT = Cross(N, T);
	Vector NxI = Cross(N, I);
	Vector difference = (NxT - miu * NxI);
	assert ( difference.Length() < RAY_EPSILON );
#endif 

	return true; 
}


bool lensElement::checkStop(CRay * outRay, const CRay * inRay ){	
	//Totally, we want to cal the O.z + t * D.z = zpos
	float stopZ = zpos;	

	//All we set is outRay->o , ->d; the other parts we use inRay
	(*outRay) = (*inRay);

	if ( abs(inRay->d.z) < RAY_EPSILON) {
		//set the ray to be invalid		
		invalidateRay(outRay);				
		return false;
	}
	float t = (stopZ - inRay->o.z ) / inRay->d.z;
	
	if (t >  (*inRay).maxt || t <  (*inRay).mint){
		//set the ray to be invalid		
		invalidateRay(outRay);				
		return false;
	}
	
	//If the intersection point is outside of the lense itself, then return false;
	float intersectionPy = (*inRay)(t).y;
	if ( abs(intersectionPy) > (aper/2) ) {
		invalidateRay(outRay);
		return false;
	}
	(*outRay).d = (*inRay).d;
	(*outRay).o = (*inRay)(t);


	return true; 
}

//This is from Image to Obejct.
void lensSystem::refractExact_ImgToObj(CRayChain * pRayChain) {
	vector<lensElement>::iterator Li, Ln;		
	//Refract for the Lens
	int index = 0; 
	CRay inRay, outRay; 
	inRay = pRayChain->rayChain[0];
	float miu = -1;  
	//As to the miu = ni/nt, we need to add the air to the leftmost side!!! also take care of the Stop???
	//This is from Image to Obejct. (lenses.end()-1 -> lenses.begin)
	Li = lenses.end();
	for (Li = lenses.end() - 1; Li >= lenses.begin(); Li--)  {
		if (Li->isActive && !Li->isStop) {
			Ln = Li - 1; 
			while (Ln->isStop) Ln = Ln - 1; 
			if (Ln < lenses.begin() )  {
				//we reach the object side air now
				miu =  Li->ndr;  // Li->ndr/1.0
			} else {
				miu = Li->ndr / Ln->ndr;
			}
			bool b_penetrated = Li->refractOnce(&outRay, &inRay, miu);
			pRayChain->rayChain.push_back(outRay);						//Now we also push_back the invalidated ray!
			index ++;			
			inRay = outRay;
			if (!b_penetrated) return;			
		} else if (Li->isActive && Li->isStop ) {
			// Check the stop. If penetrated, then nothing happen. If not, push back an invalidated Ray.
			bool b_penetrated = Li-> checkStop(&outRay, &inRay);
			pRayChain->rayChain.push_back(outRay);						//Now we also push_back the invalidated ray!
			index ++;			
			inRay = outRay;
			if (!b_penetrated) {
				//invalidateRay(outRay); 
				return;
			}
		}
	}
}

//This is from Image to Obejct.
void lensSystem::refractThick_ImgToObj(CRayChain * pRayChain) {
	
	//Suppose we already got the first ray, then we need to take care of all these stuffs. 	
	//0. got the input Point
	vector<CRay>::iterator Ri = pRayChain->rayChain.begin(); 
	Point inputP, outputP; 
	inputP = Ri->o ; 
	//drawSquare(1,0,0, inputP.z, inputP.y, inputP.z + SQUARE_WIDTH, inputP.y + SQUARE_WIDTH); 
	//1. Matrix
	float matrix_t = p1 - p2; 
	float matrix_f = f1 - p1; 
	Matrix4x4 T(1, 0, 0, 0,    0, 1, 0, 0,     0, 0, 1+ matrix_t/matrix_f, matrix_t,     0, 0 ,1/matrix_f, 1);
	
	
	float inputW = 1; 
	float outputW = -1; 	
	
	inputP.z -= p2;
	outputP.x = T.m[0][0] * inputP.x + T.m[0][1] * inputP.y + T.m[0][2] * inputP.z + T.m[0][3] * inputW;
	outputP.y = T.m[1][0] * inputP.x + T.m[1][1] * inputP.y + T.m[1][2] * inputP.z + T.m[1][3] * inputW;
	outputP.z = T.m[2][0] * inputP.x + T.m[2][1] * inputP.y + T.m[2][2] * inputP.z + T.m[2][3] * inputW;
	outputW = T.m[3][0] * inputP.x + T.m[3][1] * inputP.y + T.m[3][2] * inputP.z + T.m[3][3] * inputW;
	
	assert(outputW != 0); 
	outputP /= outputW; 
	outputP.z += p2; 

	//2. Map first Ray's Original to Object space. it should be the same of the precise mode
	//Show the Mapped Point; 
	//drawSquare(0,1,0, outputP.z, outputP.y, outputP.z + SQUARE_WIDTH, outputP.y + SQUARE_WIDTH); 
	glColor3f(1,1,1);
	
	//3. Calculate all the raychain, at most 3 ray. 
	//3.a: Get the min of p1, p2 AND ACTUAL STOP'S. 
	inputP = Ri->o;
	float intersectNearZ = p1, intersectFarZ = p2, actualStopZ = getAperture()->zpos;
	if ( abs(intersectNearZ - inputP.z) > abs(intersectFarZ - inputP.z) ){
		swap(intersectFarZ, intersectNearZ); 
		assert( abs(intersectNearZ - inputP.z) <= abs(intersectFarZ - inputP.z) );
	}
	int numberActualStop = 0; 
	if (abs(actualStopZ - inputP.z) > abs(intersectNearZ - inputP.z)) numberActualStop = 1; 
	if (abs(actualStopZ - inputP.z) > abs(intersectFarZ - inputP.z)) numberActualStop = 2; 	
	char str[20]; 
	sprintf(str, "%d", numberActualStop); 
	//drawText(10,10, SBFONT, str);	
	
	//4. Whenever we reach the actual stop, If can't penetrate the stop, then checkStop will invalidated it, and we just return; 
	CRay inRay = pRayChain->rayChain[0]; 
	CRay outRay; 
	bool b_penetrated = false; 
	if (numberActualStop == 0) {
		b_penetrated = getAperture()->checkStop(&outRay, &inRay);
		pRayChain->rayChain.push_back(outRay);
		inRay = outRay;
		if (! b_penetrated) return; 
	}
	//intersecting inRay with intersectNearZ;
	assert(inRay.d.z != 0); 
	float firstHitT = (intersectNearZ - inRay.o.z) / inRay.d.z; 
	outRay.o = inRay(firstHitT); 
	outRay.d = Vector(0,0,1); 
	pRayChain->rayChain.push_back(outRay);
	inRay = outRay;	

	if (numberActualStop == 1) {
		b_penetrated = getAperture()->checkStop(&outRay, &inRay);
		pRayChain->rayChain.push_back(outRay);
		inRay = outRay;
		if (! b_penetrated) return; 
	}
	///intersectingZ(intersectFarZ);
	assert(inRay.d.z != 0); 
	float secondHitT = (intersectFarZ - inRay.o.z) / inRay.d.z; 
	outRay.o = inRay(secondHitT); 
	outRay.d = Normalize( outputP - outRay.o);
	pRayChain->rayChain.push_back(outRay);
	inRay = outRay;	

	if (numberActualStop == 2) {
		b_penetrated = getAperture()->checkStop(&outRay, &inRay);
		pRayChain->rayChain.push_back(outRay);
		inRay = outRay;
		if (! b_penetrated) return; 
	}
	
	return; 
}

//This is from Obejct To Image
void lensSystem::refractExact_ObjToImg(CRayChain * pRayChain) {
	vector<lensElement>::iterator Li, Lp;		 //Lp means previous one
	//Refract for the Lens
	int index = 0; 
	CRay inRay, outRay; 
	inRay = pRayChain->rayChain[0];
	float miu = -1;  
	//As to the miu = ni/nt, we need to add the air to the leftmost side!!! also take care of the Stop???
	//This is from Image to Obejct. (lenses.end()-1 -> lenses.begin)	
	for (Li = lenses.begin(); Li < lenses.end() ; Li++)  {
		if (Li->isActive && !Li->isStop) {
			Lp = Li - 1; 
			while (Lp->isStop) Lp = Lp - 1; 
			if (Lp < lenses.begin() )  {
				//we reach the object side air now
				miu =  1.0 / Li->ndr;  // Li->ndr/1.0
			} else {
				miu = Lp->ndr / Li->ndr;
			}
			bool b_penetrated = Li->refractOnce(&outRay, &inRay, miu);
			pRayChain->rayChain.push_back(outRay);						//Now we also push_back the invalidated ray!
			index ++;			
			inRay = outRay;
			if (!b_penetrated) return;			
		} else if (Li->isActive && Li->isStop ) {
			// Check the stop. If penetrated, then nothing happen. If not, push back an invalidated Ray.
			bool b_penetrated = Li-> checkStop(&outRay, &inRay);
			pRayChain->rayChain.push_back(outRay);						//Now we also push_back the invalidated ray!
			index ++;			
			inRay = outRay;
			if (!b_penetrated) {
				//invalidateRay(outRay); 
				return;
			}
		}
	}
}

//Given all the setting, the OUTPUT the f1 and f2
void lensSystem::calF1F2P1P2(void) {
	//Set two parallel rays, then call the refractExact() -> We got two outRays, then we got the F1, F2
	int lensNumber = numLenses() - 1; //Since there is always one stop.
	float tempImageZ = -5000; 
	float tempObjectZ = 5000; 
	CRayChain rayChain1; // = new CRay [lensNumber + 1];	//Here + 1 means, n lenses will generate n+1 rays
	CRayChain rayChain2; // = new CRay [lensNumber + 1];
	float rayHeight = minAperture()/8;

	//10 Setup for tracing
	CRay ray1(Point(0, rayHeight, tempImageZ ), Vector(0,0, 1) );
	CRay ray2(Point(0, rayHeight, tempObjectZ), Vector(0,0, -1) );
	rayChain1.rayChain.push_back( ray1 );
	rayChain2.rayChain.push_back( ray2 );
	
	//20 Trace
	refractExact_ImgToObj(&rayChain1);
	refractExact_ObjToImg(&rayChain2);	

	//rayChain1.draw();
	//rayChain2.draw();
	
	//30. Now cal the F1, 
	//vector<CRay>::iterator Rif1, Rif2; 
	CRay  Rif1,  Rif2; 
	float oyf1, dyf1, tf1, oyf2, dyf2, tf2; 

	//Rif1 = rayChain1.rayChain.end() - 1 ;	
	Rif1 = rayChain1.rayChain[rayChain1.rayChain.size()-1];
	oyf1 = Rif1.o.y;
	dyf1 = Rif1.d.y;
	tf1 = - oyf1 / dyf1; 
	this->f1 = (Rif1)(tf1).z;
	
	//Rif2 = rayChain2.rayChain.end() - 1 ;	
	Rif2 = rayChain2.rayChain[rayChain2.rayChain.size()-1];
	oyf2 = Rif2.o.y;
	dyf2 = Rif2.d.y;
	tf2 = - oyf2 / dyf2; 
	this->f2 = (Rif2)(tf2).z;

	
	//

	//40 Cal the P1, P2. 
	//Rif1 <-> ray1 ==> P1
	//Rif2 <-> ray2 ==> P2		
	float z1 = cross2RaysAtZ( &Rif1, &ray1);
	float z2 = cross2RaysAtZ( &Rif2, &ray2);
	this->p1 = z1; 
	this->p2 = z2; 

}	  	  

float lensSystem::cross2RaysAtZ(const CRay * ray1, const CRay * ray2) {
	float resultZ = ray1->o.z;
	// O1+ t1* D1 = O2 + t2 * D2; on .y, .z, since .x is trivial
	// oy1+ t1* dy1 = oy2 + t2 * dy2;
	// oz1+ t1* dz1 = oz2 + t2 * dz2;
	float oy1 = ray1->o.y, dy1 = ray1->d.y, oy2 = ray2->o.y, dy2 = ray2->d.y;
	float oz1 = ray1->o.z, dz1 = ray1->d.z, oz2 = ray2->o.z, dz2 = ray2->d.z;
	if ( (dy1) != 0) {
		float leftside = oz1 + (oy2-oy1)/dy1*dz1 -oz2; 
		float rightside = dz2 - dz1*dy2/dy1;
		float t2  = leftside/rightside;
		resultZ = (*ray2)(t2).z;
	} else {
		//dy1 == 0, since all the ray are initialized to be correct, so dy1 cannot be 0!!!
		assert(false);
		//t1*dz1 = oz2 + (oy2-oy1)/oy2*dz2 - oz1;
		if ( dy2 == 0) {
			if (oz1 == oz2) ; //The same line
			else ;  //never intersect
		}else if (dz1 ==0) {
			//trivial, |ray1.d| = 0; 
		}
	}
	return resultZ; 
}

void lensSystem::initlsystem() {
	///This is for the rear part. And after we get the exit pupil, we want to 	
	////This is only for the whole system
	//this->printLenses();
	
	this->calF1F2P1P2();	//f1, f2, p1, p2	
	this->setMaxFSMinFS();  //fstop and m_flefective_aperture, minFS, maxFS will be set
	this->calExitPupil();
	this->caliPlane(2000);
	this->printConfig();
	//lvFStop = lsystem.minFS;
}

Point lensSystem::calExitPupil() {
	//Save a lot of parameter, since this function will change all of them 
	float pref1 = this->f1, pref2 = this->f2, prep1 = this->p1, prep2 = this->p2; 
	float prefstop = fstop, pre_m_fleffective_aperture = m_fleffective_aperture;
	float premaxFS = maxFS, preminFS = minFS; 
	
	vector<lensElement>::iterator  Li, Ln, Lend;

	//15.Set the active part as the rear part.
	Li = lenses.begin();	
	while (! Li->isStop && Li < lenses.end() )  {
		Li->isActive = false; 
		Li ++;
	}
	
	//20. cal all the f1, f2, ....
	this->calF1F2P1P2();	//f1, f2, p1, p2	

	//30 cal the Image Point of the stopPoint for the only rear part
	//Using p1, p2, f2, m_fleffective_aperture [changed by calAperture()] getAperture()->zpos[Fixed]
	//a. cal the matrix
	//t = P' - P
	float t = p2-p1;
	 
	Point result; 
	float ix, iy, iz, iw, ex, ey, ez, ew;  // (i-p1) = T (e-p1) ; i = T(e-p1) + p1;
	
	float ff = f2 - p2; 
	//assert ((abs(ff) - abs(f1-p1)) <RAY_EPSILON); 

	Matrix4x4 T(1, 0, 0, 0,       0, 1, 0, 0,     0, 0, 1+ t/ff, t,     0, 0 ,1/ff, 1);
	ex = 0; ey = m_fleffective_aperture/2; ez = getAperture()->zpos; ew = 1;
	ez = ez - p1; 
	ix = T.m[0][0] * ex + T.m[0][1] * ey + T.m[0][2] * ez + T.m[0][3] * ew;
	iy = T.m[1][0] * ex + T.m[1][1] * ey + T.m[1][2] * ez + T.m[1][3] * ew;
	iz = T.m[2][0] * ex + T.m[2][1] * ey + T.m[2][2] * ez + T.m[2][3] * ew;
	iw = T.m[3][0] * ex + T.m[3][1] * ey + T.m[3][2] * ez + T.m[3][3] * ew;
	
	assert(iw != 0);
	
	result = Point (ix/iw, iy/iw, iz/iw);
	result.z = result.z + p1;
	//20 with stopPoint, I can find the exit pupil.	
	this->pupil = result; 
	
	//85.Set the rear part as active.
	Li = lenses.begin();	
	while (! Li->isStop && Li < lenses.end() )  {
		Li->isActive = true; 
		Li ++;
	}
	//this->printConfig();
	//95. Recover all the parameters.
	f1 = pref1; f2 = pref2; p1 = prep1; p2 =prep2; 
	assert (this->fstop == prefstop); 
	assert (this->m_fleffective_aperture ==  pre_m_fleffective_aperture); 
	assert (this->minFS == preminFS);
	assert (this->maxFS == premaxFS);
	//this->printConfig();
	return result;

}

float lensSystem::caloPlane(float input_iPlane){
	// 1 / z' - 1 / z = 1 / f'
	// z' = z2 = iPlane - P'[p2]  ; iPlane = z2 + P'
	// z  = z1 = oPlane - P [p1]  ; oPlane = z1 + P
	// f' = ff = F'-P' = f2-p2
	this->iPlane = input_iPlane;
	if (this->iPlane >= this->f2) this->iPlane = this->f2 - 1;

	float z2 = iPlane - this->p2; 
	float ff = this->f2 - this->p2; 
	assert(ff != 0); 
	assert(z2 != 0); 
	float reciprocalff = 1 / ff; 
	float reciprocalz2 = 1 / z2; 
	float reciprocalz1 = - reciprocalff + reciprocalz2; 
	assert(reciprocalz1 != 0);
	float tempz1 = 1 / reciprocalz1; 
	this->oPlane = tempz1 + this->p1;
	if(this->oPlane <= this->f1) this->oPlane = this->f1 + 1;

	return this->oPlane;
}

float lensSystem::caliPlane(float input_oPlane) {
	// 1 / z' - 1 / z = 1 / f'
	// z' = z2 = iPlane - P'[p2]  ; iPlane = z2 + P'
	// z  = z1 = oPlane - P [p1]
	// f' = ff = F'-P' = f2-p2
	this->oPlane = input_oPlane;
	if(this->oPlane <= this->f1) this->oPlane = this->f1 -1 ;

	float z1 = oPlane - this->p1; 
	float ff = this->f2 - this->p2; 
	assert(ff != 0); 
	assert(z1 != 0); 
	float reciprocalff = 1 / ff; 
	float reciprocalz1 = 1 / z1; 
	float reciprocalz2 = reciprocalff + reciprocalz1; 
	assert(reciprocalz2 != 0);
	float tempz2 = 1 / reciprocalz2; 
	this->iPlane = tempz2 + this->p2;

	if (this->iPlane >= this->f2) this->iPlane = this->f2 +1;
	return this->iPlane;
}


//Set the fstop and m_fleffective_aperture
float lensSystem::calAperture(float input_fstop) {
	vector<lensElement>::iterator  Li, Ln, Lend;
	vector<CRay>::iterator  Ri; 
	CRayChain rayChain4pupil;
	//float effective_aperture;	//The result of the calculation. if not penetrate, then -1
	//10 Get the fstop zpos, determine the aperture size. 
	fstop = input_fstop;	

	//15.Set the active part as the rear part.
	Li = lenses.begin();	
	while (! Li->isStop && Li < lenses.end() )  {
		Li->isActive = false; 
		Li ++;
	}

	//20 Set the stopPoint, 
	//a. get the fd, 
	Lend = lenses.end() - 1; 
	float fd = Lend ->zpos - f2;
	//b. get the sin(theta).
	float sin_theta = 0.5 / fstop; 
	//c. trace a ray at theta angel.
	CRay rayAtTheta (Point(0,0,f2),Vector(0, sin_theta,sqrtf(1-sin_theta*sin_theta) ) );
	rayChain4pupil.rayChain.push_back(rayAtTheta);
	refractExact_ImgToObj(&rayChain4pupil);
	
	//rayChain4pupil.draw();	
	Ri = rayChain4pupil.rayChain.end() - 1;	
	//d. intersect this ray with the stop, now get the effective aperture. 	
	if (Ri->maxt < Ri->mint) 
		m_fleffective_aperture = -1; 
	else 
		m_fleffective_aperture = Ri ->o.y *2;  //since I always send to the y+, so it is always the positive.
	
	//95 Set back the active part as usual
	Li = lenses.begin();	
	while (! Li->isStop && Li < lenses.end() )  {
		Li->isActive = true; 
		Li ++;
	}	
	return m_fleffective_aperture;
}

void lensSystem::setMaxFSMinFS(void) {
	this->minFS = 0.6; 
	this->maxFS = 30; 
	while (calAperture(minFS) < 0 )	this->minFS += 0.1;
}



////////////////////////////////////////////////////////////////////
//Interface
void applyTransforms(void);   //from lensviw.cpp
void CRayChain::draw()
{
	applyTransforms();
	vector<CRay>::iterator Ri, Rn;
	//If the last ray is invalid, then abort;
	Ri = this->rayChain.end() - 1; 
	if (Ri->maxt < Ri->mint ) return; 

	//draw all the rays except last one
	for (Ri= this->rayChain.begin(); Ri < this->rayChain.end()-1; Ri ++) {
		//glColor3f(0.15,1 - deltaColor, deltaColor);
		//deltaColor += stepColor;
		
		glColor3f(  rand()/(float)RAND_MAX,rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);		
		Rn = Ri + 1;
		drawLine( Ri->o.z, Ri->o.y, Rn->o.z, Rn->o.y);				
	}
	//draw the last ray.
	float extentedLen = 1000; 
	if (Ri->maxt >= Ri->mint) {
		glColor3f(  rand()/(float)RAND_MAX,rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);		
		drawLine( Ri->o.z, Ri->o.y, (*Ri)(extentedLen).z, (*Ri)(extentedLen).y );				
	} else {
		//assert(false);
		glColor3f(  rand()/(float)RAND_MAX,rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);		
		glLineWidth(2);
		drawLine( Ri->o.z, Ri->o.y, (*Ri)(extentedLen).z, (*Ri)(extentedLen).y );				
		glLineWidth(1);

	}
}

void CRayChain::print()
{	
	vector<CRay>::iterator Ri;	
	for (Ri= this->rayChain.begin(); Ri < this->rayChain.end(); Ri ++) {				
		Ri->print();
	}
}
