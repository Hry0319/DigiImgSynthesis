
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
    z = new float[2*nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    min_z = max_z = z[0];
    for (int index = 0; index < nx*ny; index++){
		if (min_z > z[index]) min_z = z[index];
		if (max_z < z[index]) max_z = z[index];
	}

//	BBox bounds = BBox(Point(0,0,min_z), Point(1,1,max_z));
//	Vector _delta = bounds.pMax - bounds.pMin;
//	voxelwidth = Vector(_delta[0]/width, _delta[1]/height, 0.0f);
	
	voxelwidth = Vector(1.f/nx, 1.f/ny, 0.0f);

    nVoxels[0] = nx -1;
    nVoxels[1] = ny -1;
	nVoxels[2] = 1;


	float cubeRoot = 3.f * powf(1.f, 1.f/3.f);
    float voxelsPerUnitDist = cubeRoot;
    for (int axis = 0; axis < 2; ++axis) {
        nVoxels[axis] = Clamp(nVoxels[axis], 1, 64);
    }

	InitVertexNormals();
}


Heightfield2::~Heightfield2() {
    delete[] z;
    delete[] vertexNormals;
    delete[] points;
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
	int 		voxel2posX = nx-1;
	int 		voxel2posY = ny-1;
  	int 		voxel2posZ = 1;

	vertexNormals = new Normal[nx*ny];
	points        = new Point[nx*ny];

	for ( int j = 0; j < ny; j++)
	{
		for ( int i = 0; i < nx; i++)
		{
			p[4] = Point(i/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);
			points[i+j*nx] = p[4];

			//
			// inside middle
			//
			if (i > 0 && j > 0 && i < nx && j < ny)
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
				if (i == nx && j == ny)
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
				else if (i == nx && j <  ny)
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
				else if (i <  nx && j == ny)
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
			printf("z:%2.5f    ",vertexNormals[i+j*nx].x,vertexNormals[i+j*nx].y,vertexNormals[i+j*nx].z);

		}
	}
}

bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const{
	Ray ray;
	(*WorldToObject)(r, &ray);

	float rayT;
	BBox bounds = ObjectBound();
    if (bounds.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if (!bounds.IntersectP(ray, &rayT))
        return false;

    // FirstHitPoint
    // ray's member : Point operator()(float t) const { return o + d * t; }
    Point gridIntersect = ray(rayT);
    //
	// Set up 2D DDA for ray (ref to GridAccel)
	//
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
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

	//
	// Walk grid for shadow ray
	//
    bool hitSomething = false;
    for (;;) {
		int i = Pos[0], j = Pos[1];

		/**
		 *        4 points
		 *        2 triangles
         *
         *
		 *        (i,j)    (i+1,j)
		 *         .0 _______.1
		 *           |\      |
         *           |  \    |
		 *           |    \  |
	 	 *           |______\|
	  	 *         .2        .3
		 *       (i,j+1)   (i+1,j+1)
		 *
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
			Normal normals[4] = {
						vertexNormals[j*nx + i],
						vertexNormals[j*nx + i+1],
						vertexNormals[(j+1)*nx + i],
						vertexNormals[(j+1)*nx + i+1]};


			Point  IntersectTriangles[3] = {triangle[0], triangle[1], triangle[3]};
			Normal IntersectNormals[3]   = {normals[0] , normals[1] , normals[3] };

			Intersection in1, in2;
			float tHit1,tHit2;
			bool  isHit1,isHit2;
			isHit1 = TriangleIntersect(ray, &in1, IntersectTriangles, IntersectNormals, &tHit1);
			IntersectTriangles[1] = triangle[2];
			IntersectNormals[1]   = normals[2];
			isHit2 = TriangleIntersect(ray, &in2, IntersectTriangles, IntersectNormals, &tHit2);

			if (!isHit1 && !isHit2) {
				hitSomething = false;
			}
			else
			{
				if (isHit1 && isHit2) {
					if (tHit1 < tHit2) {
						*dg 		= in1.dg;
						*rayEpsilon = in1.rayEpsilon;
						isHit2 = false;
						*tHit = tHit1;
					} else {
						*dg 		= in2.dg;
						*rayEpsilon = in2.rayEpsilon;
						isHit1 = false;
						*tHit = tHit2;
					}
				} else if (isHit1) {
					*dg 		= in1.dg;
					*rayEpsilon = in1.rayEpsilon;
					*tHit = tHit1;
				} else {
					*dg 		= in2.dg;
					*rayEpsilon = in2.rayEpsilon;
					*tHit = tHit2;
				}
				hitSomething = true;
			}

		}


		// in heightfields, there will be no overlap
		if (hitSomething) break;

        // Advance to next voxel
        // Find stepAxis for stepping to next voxel
		int stepAxis = NextCrossingT[0] < NextCrossingT[1] ? 0 : 1;
		if (ray.maxt < NextCrossingT[stepAxis])
			break;
		Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
			break;
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }

	return hitSomething;
}

bool Heightfield2::IntersectP(const Ray &r) const {
//	float t0 = r.mint, t1 = r.maxt;
//	for(int axis = 0; axis < 3; ++axis)
//	{
//		float invRayDir = 1.f / r.d[axis];
//		float tNear = (nVoxels[axis] - r.o[axis]) * invRayDir;
//		float tFar  = (nVoxels[axis] - r.o[axis]) * invRayDir;
//
//		if (tNear > tFar) swap(tNear, tFar);
//
//		t0 = tNear > t0 ? tNear : t0;
//		t1 = tFar < t1 ? tFar : t1;
//		if (t0 > t1) return false;
//	}
//	return false;
	Ray ray;
	(*WorldToObject)(r, &ray);

	float rayT;
	BBox bounds = ObjectBound();
    if (bounds.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if (!bounds.IntersectP(ray, &rayT))
        return false;

    // FirstHitPoint
    // ray's member : Point operator()(float t) const { return o + d * t; }
    Point gridIntersect = ray(rayT);
    //
	// Set up 2D DDA for ray (ref to GridAccel)
	//
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = pos2Voxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
        // + Handle ray with "positive" direction for voxel stepping
            NextCrossingT[axis] = rayT + (voxel2Pos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = 	(float(voxelwidth[axis]) / ray.d[axis] );
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else
        {
        // - Handle ray with "negative" direction for voxel stepping
            NextCrossingT[axis] = rayT + (voxel2Pos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = -(float(voxelwidth[axis]) / ray.d[axis] );
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

	//
	// Walk grid for shadow ray
	//
    bool hitSomething = false;
    for (;;) {
		int i = Pos[0], j = Pos[1];

		Point *triangle = new Point[4];

		triangle[0] = Point(voxel2Pos(i,0)		, voxel2Pos(j,1)	,z[i+ j*nx]);
		triangle[1] = Point(voxel2Pos(i+1,0)	, voxel2Pos(j,1)	,z[(i+1) + j*nx]);
		triangle[2] = Point(voxel2Pos(i,0)		, voxel2Pos(j+1,1)	,z[i + (j+1)*nx]);
		triangle[3] = Point(voxel2Pos(i+1,0)	, voxel2Pos(j+1,1)	,z[(i+1) + (j+1)*nx]);

		BBox bounds(Union(BBox(triangle[0], triangle[1]), BBox(triangle[2], triangle[3])));

		if (bounds.IntersectP(ray))
		{
			Point  IntersectTriangles[3] = {triangle[0], triangle[1], triangle[3]};
			bool  isHit1,isHit2;

			isHit1 = TriangleIntersectP(ray, IntersectTriangles);
			IntersectTriangles[1] = triangle[2];
			isHit2 = TriangleIntersectP(ray, IntersectTriangles);

			if (!isHit1 && !isHit2) {
				hitSomething = false;
			}
			else
			{
				hitSomething = true;
			}
		}

		// in heightfields, there will be no overlap
		if (hitSomething) break;

        // Advance to next voxel
        // Find _stepAxis_ for stepping to next voxel
		int stepAxis = (NextCrossingT[0] < NextCrossingT[1]) ? 0 : 1;
		if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
	return hitSomething;

}


void Heightfield2::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const {
	//dg.shape->GetShadingGeometry(obj2world,dg,dgShading);
	* dgShading = dg;
}

BBox Heightfield2::ObjectBound() const {
	return BBox(Point(0,0,min_z), Point(1,1,max_z));
}

bool Heightfield2::CanIntersect() const {
    return true;
}

bool Heightfield2::TriangleIntersect(const Ray &r, Intersection *instect, const Point *triangle, const Normal *normals, float *tHit) const{
	Ray ray;
	(*ObjectToWorld)(r, &ray);

    const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);

	float divisor = Dot(s1, e1);

    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
	float uvs[6] =	{
					p1.x, p1.y,
					p2.x, p2.y,
					p3.x, p3.y
					};

	// Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*uvs[0] + b1*uvs[2] + b2*uvs[4];
    float tv = b0*uvs[1] + b1*uvs[3] + b2*uvs[5];





	// Compute deltas for triangle partial derivatives
	/*float du1 = uvs[0][0] - uvs[2][0];
      float du2 = uvs[1][0] - uvs[2][0];
      float dv1 = uvs[0][1] - uvs[2][1];
      float dv2 = uvs[1][1] - uvs[2][1];*/
    float du1 = uvs[0] - uvs[4];
    float du2 = uvs[2] - uvs[4];
    float dv1 = uvs[1] - uvs[5];
    float dv2 = uvs[3] - uvs[5];
    Vector dp1 = p1 - p3, dp2 = p2 - p3;

    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

	// manipulate the normal
	/*
	Vector normal00= Vector(normals[0]);
	Vector normal01= Vector(normals[1]);
	Vector normal02= Vector(normals[2]);
	Vector avg_normal = normal00 * b0 + normal01 * b1 + normal02 * b2;
	Vector temp_tangent(1,0,0);
	Vector temp_surface(1,0,0);

	avg_normal = Normalize(avg_normal);
	temp_surface = Normalize(Cross(avg_normal, temp_tangent) ); //now surface is |_ avg_normal
	temp_tangent = Normalize(Cross(temp_surface, avg_normal) );

	dpdu = temp_tangent;
	dpdv = temp_surface;
	*/

    // Fill in _DifferentialGeometry_ from triangle hit
    instect->dg = DifferentialGeometry( ray(t),
										dpdu,
										dpdv,
		                                Normal(0,0,0),
										Normal(0,0,0),
		                                tu, tv,
										this);
    *tHit = t;
    instect->rayEpsilon = 1e-3f * *tHit;

	return true;
}


bool Heightfield2::TriangleIntersectP(const Ray &r, const Point *triangle) const {
	Ray ray;
	(*ObjectToWorld)(r, &ray);

    const Point &p1 = (*ObjectToWorld)(triangle[0]);
    const Point &p2 = (*ObjectToWorld)(triangle[1]);
    const Point &p3 = (*ObjectToWorld)(triangle[2]);
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);

	float divisor = Dot(s1, e1);

    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

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


