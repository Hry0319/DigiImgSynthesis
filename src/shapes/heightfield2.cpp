
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"
#include "error.h"
#include <stdio.h>

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(
				const Transform *o2w,
				const Transform *w2o,
        		bool ro,
        		int x,
        		int y,
        		const float *zs
        		) : Shape(o2w, w2o, ro)
{
    nx = x;
    ny = y;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    min_z = max_z = z[0];
    for (int index = 0; index < nx*ny; index++){
		if (min_z > z[index]) min_z = z[index];
		if (max_z < z[index]) max_z = z[index];
	}
	width  = nx;
	height = ny;
	depth  = 1;

//	BBox bounds = BBox(Point(0,0,min_z), Point(1,1,max_z));
//	Vector _delta = bounds.pMax - bounds.pMin;
//	voxelwidth = Vector(_delta[0]/width, _delta[1]/height, 0.0f);
	voxelwidth = Vector(1.f/width, 1.f/height, 0.0f);

    nVoxels[0] = nx - 1;
    nVoxels[1] = ny - 1;
	nVoxels[2] = 1;

	InitVertexNormals();
	InitTriangles();
}


Heightfield2::~Heightfield2() {
    delete[] z;
}


//
// Heightfield2::InitTriangles()
//
void Heightfield2::InitTriangles() {

}
//
// Heightfield2::InitNormals()
//
/**
-*
 *	    0  1  2
 *      3  4  5
 *      6  7  8
-*
**/
void Heightfield2::InitVertexNormals() {
	Point       p[9];
	int 		voxel2posX = width-1;
	int 		voxel2posY = height-1;
  	int 		voxel2posZ = 1;

	vertexNormals = new Normal[nx*ny];
	points        = new Point[nx*ny];

	for ( int j = 0; j < ny; j++)
	{
		for ( int i = 0; i < nx; i++)
		{
			p[4] = Point(i/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);
			points[i+j*nx] = p[4];
			//printf("/// %f %f %f  \n", points[i+j*nx].x, points[i+j*nx].y, points[i+j*nx].z);

			//
			// inside middle
			//
			if (i > 0 && j > 0 && i < width && j < height)
			{
				p[0] = Point( (i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]);
				p[1] = Point( (i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
				p[2] = Point( (i+1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i+1)+(j-1)*nx]);
				p[3] = Point( (i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]);
				p[5] = Point( (i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]);
				p[6] = Point( (i-1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i-1)+(j+1)*nx]);
				p[7] = Point( (i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]);
				p[8] = Point( (i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]);
			}
			else
			{
				//
				// left top
				//
				if (i == 0 &&  j == 0)
				{
					p[0] = p[4];
					p[1] = p[4];
					p[2] = p[4];
					p[3] = p[4];
					p[6] = p[4];

					p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]);
					p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]);
					p[8] = Point((i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]);
				}
				else if (i == 0 && j > 0 )
				{
					p[0] = Point(i/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
					p[3] = Point(i/(float)voxel2posX, (j  )/(float)voxel2posY, z[i+j*nx]);
					p[6] = Point(i/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]);

					p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
					p[2] = Point((i+1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i+1)+(j-1)*nx]);
					p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]);
					p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]);
					p[8] = Point((i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]);
				}
				else if (i > 0  && j == 0)
				{
					p[0] = Point((i-1)/(float)voxel2posX, j/(float)voxel2posY, z[(i-1)+j*nx]);
					p[1] = Point((i  )/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]);
					p[2] = Point((i+1)/(float)voxel2posX, j/(float)voxel2posY, z[(i+1)+j*nx]);

					p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]);
					p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]);
					p[6] = Point((i-1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i-1)+(j+1)*nx]);
					p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]);
					p[8] = Point((i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]);
				}
				//
				// right button
				//
				if (i == width && j == height)
				{
					p[2] = p[4];
					p[5] = p[4];
					p[6] = p[4];
					p[7] = p[4];
					p[8] = p[4];

					p[0] = Point((i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]);
					p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
					p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]);
				}
				else if (i == width && j <  height)
				{
					p[2] = Point(i/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
					p[5] = Point(i/(float)voxel2posX, (j  )/(float)voxel2posY, z[i+j*nx]);
					p[8] = Point(i/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j-1)*nx]);

					p[0] = Point((i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]);
					p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
					p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]);
					p[6] = Point((i-1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i-1)+(j+1)*nx]);
					p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]);
				}
				else if (i <  width && j == height)
				{
					p[6] = Point((i-1)/(float)voxel2posX, j/(float)voxel2posY, z[(i-1)+j*nx]);
					p[7] = Point((i  )/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]);
					p[8] = Point((i+1)/(float)voxel2posX, j/(float)voxel2posY, z[(i+1)+j*nx]);

					p[0] = Point((i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]);
					p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]);
					p[2] = Point((i+1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i+1)+(j-1)*nx]);
					p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]);
					p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]);
				}
			}
			vertexNormals[i+j*nx] = Normal(	Normalize(
													Cross( (p[0]-p[4]) ,(p[1]-p[4]) ) + \
													Cross( (p[1]-p[4]) ,(p[2]-p[4]) ) + \
													Cross( (p[2]-p[4]) ,(p[5]-p[4]) ) + \
													Cross( (p[5]-p[4]) ,(p[8]-p[4]) ) + \
													Cross( (p[8]-p[4]) ,(p[7]-p[4]) ) + \
													Cross( (p[7]-p[4]) ,(p[6]-p[4]) ) + \
													Cross( (p[6]-p[4]) ,(p[3]-p[4]) ) + \
													Cross( (p[3]-p[4]) ,(p[0]-p[4]) )
													)
												)*(-1);

		}
	}
}

bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const{
	Ray ray;
//	(*WorldToObject)(ray, &rayW2O);
	(*WorldToObject)(r, &ray);

	float rayT;
	BBox bounds = ObjectBound();
    if (bounds.Inside(ray(ray.mint)))
    {
        rayT = ray.mint;
    }
    else if (!bounds.IntersectP(ray, &rayT))
    {
        return false;
	}

    // FirstHitPoint
    // ray's member : Point operator()(float t) const { return o + d * t; }
    Point gridIntersect = ray(rayT);
    //
	// Set up 3D DDA for ray (ref to GridAccel)
	//
    float NextCrossingT[3], DeltaT[3];
    int Step[3], Out[3], Pos[3];
    for (int axis = 0; axis < 3; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = pos2Voxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
        // + Handle ray with "positive" direction for voxel stepping
            NextCrossingT[axis] = rayT + (voxel2Pos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = 	(voxelwidth[axis] / ray.d[axis] );
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else
        {
        // - Handle ray with "negative" direction for voxel stepping
            NextCrossingT[axis] = rayT + (voxel2Pos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = -(voxelwidth[axis] / ray.d[axis] );
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

	// Walk grid for shadow ray
    bool hitSomething = false;
	Intersection isect;
    for (;;) {
		int i = Pos[0], j = Pos[1];

		/**
		 *      4 points
		 *      2 triangles
         *
         *
		 *    (i,j)    (i+1,j)
		 *     .0 _______.1
		 *       |\      |
         *       |  \    |
		 *       |    \  |
	 	 *       |______\|
	  	 *     .2        .3
		 *   (i,j+1)   (i+1,j+1)
		 *
		**/
		Point *triangle = new Point[4];

		triangle[0] = Point(voxel2Pos(i,0)		, voxel2Pos(j,1)	,z[i+ j*nx]);
		triangle[1] = Point(voxel2Pos(i+1,0)	, voxel2Pos(j,1)	,z[(i+1) + j*nx]);
		triangle[2] = Point(voxel2Pos(i,0)		, voxel2Pos(j+1,1)	,z[i + (j+1)*nx]);
		triangle[3] = Point(voxel2Pos(i+1,0)	, voxel2Pos(j+1,1)	,z[(i+1) + (j+1)*nx]);

		BBox bounds(Union(BBox(triangle[0], triangle[1]), BBox(triangle[2], triangle[3])));
		if (bounds.IntersectP(ray))
		{
			int vptr[6] = {	0,1,3,0,2,3 };

			Normal normals[4] = {
						vertexNormals[j*nx + i],
						vertexNormals[j*nx + i+1],
						vertexNormals[(j+1)*nx + i],
						vertexNormals[(j+1)*nx + i+1]};

			TriangleMesh *triMesh = new TriangleMesh(ObjectToWorld, WorldToObject, ReverseOrientation,
					2, 4, vptr, triangle, normals, NULL, NULL, NULL);
			Triangle *triangle1 = new Triangle(ObjectToWorld, WorldToObject, ReverseOrientation, triMesh, 0);
			Triangle *triangle2 = new Triangle(ObjectToWorld, WorldToObject, ReverseOrientation, triMesh, 1);

			float tHit1,tHit2;
			Intersection in1, in2;
			bool tri1 = triangle1->Intersect((*ObjectToWorld)(ray), &(tHit1), &(in1.rayEpsilon), &(in1.dg));
			bool tri2 = triangle2->Intersect((*ObjectToWorld)(ray), &(tHit2), &(in2.rayEpsilon), &(in2.dg));

			if (!tri1 && !tri2) {
				hitSomething = false;
			}
			else
			{
				if (tri1 && tri2) {
					if (tHit1 < tHit2) {
						isect = in1;
						tri2 = false;
						*tHit = tHit1;
					} else {
						isect = in2;
						tri1 = false;
						*tHit = tHit2;
					}
				} else if (tri1) {
					isect = in1;
					*tHit = tHit1;
				} else {
					isect = in2;
					*tHit = tHit2;
				}
				hitSomething = true;
			}
		}

		// in heightfields, there will be no overlap
		if (hitSomething) break;

        // Advance to next voxel
        // Find _stepAxis_ for stepping to next voxel
        int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
                   ((NextCrossingT[0] < NextCrossingT[2]) << 1) +
                   ((NextCrossingT[1] < NextCrossingT[2]));
		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
        int stepAxis = cmpToAxis[bits];
		if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }

	*dg = isect.dg;
	*rayEpsilon = isect.rayEpsilon;

	return hitSomething;
}

bool Heightfield2::IntersectP(const Ray &ray) const {
	//float tHit, rayEpsilon;
	//DifferentialGeometry dg;
	//return Intersect(ray, &tHit, &rayEpsilon, &dg);
	return false;
}
bool Heightfield2::IntersectP(const Ray &ray, float *hit0, float *hit1) const{
	float t0 = ray.mint, t1 = ray.maxt;
	for(int axis = 0; axis < 3; ++axis)
	{
		float invRayDir = 1.f / ray.d[axis];
		float tNear = (nVoxels[axis] - ray.o[axis]) * invRayDir;
		float tFar  = (nVoxels[axis] - ray.o[axis]) * invRayDir;

		if (tNear > tFar) swap(tNear, tFar);

		t0 = tNear > t0 ? tNear : t0;
		t1 = tFar < t1 ? tFar : t1;
		if (t0 > t1) return false;
	}
	if (hit0) *hit0 = t0;
	if (hit1) *hit1 = t1;
	return false;
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const {
	dg.shape->GetShadingGeometry(obj2world,dg,dgShading);
//	* dgShading = dg;
}


BBox Heightfield2::ObjectBound() const {
//    BBox objectBounds;
//    for (int i = 0; i < nx*ny; i++)
//        objectBounds = Union(objectBounds, (*WorldToObject)(points[i]));
//	return objectBounds;
	return BBox(Point(0,0,min_z), Point(1,1,max_z));
}


bool Heightfield2::CanIntersect() const {
    return true;
}

Heightfield2 *CreateHeightfield2Shape (const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}

