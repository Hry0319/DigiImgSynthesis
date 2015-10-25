
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

float tHit1,tHit2;

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


    nVoxels[0] = nx - 1;
    nVoxels[1] = ny - 1;
	nVoxels[2] = 1;
//	for (int axis = 0; axis < 3; ++axis) {
//        nVoxels[axis] = Clamp(nVoxels[axis], 1, 64);
//    }

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
void Heightfield2::InitVertexNormals() {
	Point          p[9];
//	Vector         vector[8];


	int voxel2posX = width-1;
	int voxel2posY = height-1;
  	int voxel2posZ = 1;

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
			// left top
			//
			if (i == 0 &&  j == 0)
			{
				p[0] = p[4];
				p[1] = p[4];
				p[2] = p[4];
				p[3] = p[4];
				p[6] = p[4];

				p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);
				p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]/(float)voxel2posZ);
				p[8] = Point((i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]/(float)voxel2posZ);
			}
			else if (i == 0 && j > 0 )
			{
				p[0] = Point(i/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[3] = Point(i/(float)voxel2posX, (j  )/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);
				p[6] = Point(i/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]/(float)voxel2posZ);

				p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[2] = Point((i+1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i+1)+(j-1)*nx]/(float)voxel2posZ);
				p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);
				p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]/(float)voxel2posZ);
				p[8] = Point((i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]/(float)voxel2posZ);
			}
			else if (i > 0  && j == 0)
			{
				p[0] = Point((i-1)/(float)voxel2posX, j/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
				p[1] = Point((i  )/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);
				p[2] = Point((i+1)/(float)voxel2posX, j/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);

				p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
				p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);
				p[6] = Point((i-1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i-1)+(j+1)*nx]/(float)voxel2posZ);
				p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]/(float)voxel2posZ);
				p[8] = Point((i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]/(float)voxel2posZ);
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

				p[0] = Point((i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]/(float)voxel2posZ);
				p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
			}
			else if (i == width && j <  height)
			{
				p[2] = Point(i/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[5] = Point(i/(float)voxel2posX, (j  )/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);
				p[8] = Point(i/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);

				p[0] = Point((i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]/(float)voxel2posZ);
				p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
				p[6] = Point((i-1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i-1)+(j+1)*nx]/(float)voxel2posZ);
				p[7] = Point((i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]/(float)voxel2posZ);
			}
			else if (i <  width && j == height)
			{
				p[6] = Point((i-1)/(float)voxel2posX, j/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
				p[7] = Point((i  )/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);
				p[8] = Point((i+1)/(float)voxel2posX, j/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);

				p[0] = Point((i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]/(float)voxel2posZ);
				p[1] = Point((i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[2] = Point((i+1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i+1)+(j-1)*nx]/(float)voxel2posZ);
				p[3] = Point((i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
				p[5] = Point((i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);
			}

			//
			// inside middle
			//
			if (i > 0 && j > 0 && i < width && j < height)
			{
				p[0] = Point( (i-1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i-1)+(j-1)*nx]/(float)voxel2posZ);
				p[1] = Point( (i  )/(float)voxel2posX, (j-1)/(float)voxel2posY, z[i+(j-1)*nx]/(float)voxel2posZ);
				p[2] = Point( (i+1)/(float)voxel2posX, (j-1)/(float)voxel2posY, z[(i+1)+(j-1)*nx]/(float)voxel2posZ);
				p[3] = Point( (i-1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i-1)+j*nx]/(float)voxel2posZ);
				p[5] = Point( (i+1)/(float)voxel2posX, (j  )/(float)voxel2posY, z[(i+1)+j*nx]/(float)voxel2posZ);
				p[6] = Point( (i-1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i-1)+(j+1)*nx]/(float)voxel2posZ);
				p[7] = Point( (i  )/(float)voxel2posX, (j+1)/(float)voxel2posY, z[i+(j+1)*nx]/(float)voxel2posZ);
				p[8] = Point( (i+1)/(float)voxel2posX, (j+1)/(float)voxel2posY, z[(i+1)+(j+1)*nx]/(float)voxel2posZ);
			}
			vertexNormals[i+j*nx] = Normal( Normalize(	Cross( (p[0]-p[4]) ,(p[1]-p[4]) ) + \
												Cross( (p[1]-p[4]) ,(p[2]-p[4]) ) + \
												Cross( (p[2]-p[4]) ,(p[5]-p[4]) ) + \
												Cross( (p[5]-p[4]) ,(p[8]-p[4]) ) + \
												Cross( (p[8]-p[4]) ,(p[7]-p[4]) ) + \
												Cross( (p[7]-p[4]) ,(p[6]-p[4]) ) + \
												Cross( (p[6]-p[4]) ,(p[3]-p[4]) ) + \
												Cross( (p[3]-p[4]) ,(p[0]-p[4]) )
												))*(-1);

		}
	}
}

bool Heightfield2::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const{
	Ray rayW2O;
	(*WorldToObject)(ray, &rayW2O);

	float rayT;
	BBox bounds = ObjectBound();

    if (bounds.Inside(rayW2O(rayW2O.mint)))
    {
        rayT = rayW2O.mint;
    }
    else if (!bounds.IntersectP(rayW2O, &rayT))
    {
        return false;
	}


    Point gridIntersect = rayW2O(rayT);
	//
	// ray position ??
	//

	// Set up 3D DDA for ray
    float NextCrossingT[3], DeltaT[3];
    int Step[3], Out[3], Pos[3];
    for (int axis = 0; axis < 3; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = pos2Voxel(gridIntersect, axis);
        if (rayW2O.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxel2Pos(Pos[axis]+1, axis) - gridIntersect[axis]) / rayW2O.d[axis];
//            DeltaT[axis] = widthV[axis] / ray.d[axis];
            DeltaT[axis] = 1.f / (nVoxels[axis] * rayW2O.d[axis]);
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxel2Pos(Pos[axis], axis) - gridIntersect[axis]) / rayW2O.d[axis];
//            DeltaT[axis] = -widthV[axis] / ray.d[axis];
            DeltaT[axis] = -1.f / (nVoxels[axis] * rayW2O.d[axis]);
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

	// Walk grid for shadow ray
    bool hitSomething = false;
	Intersection isect;
    for (;;) {
		int i = Pos[0], j = Pos[1];
		Point TL(voxel2Pos(i,0),   voxel2Pos(j,1),   getZ(i,   j));
		Point TR(voxel2Pos(i+1,0), voxel2Pos(j,1),   getZ(i+1, j));
		Point BR(voxel2Pos(i+1,0), voxel2Pos(j+1,1), getZ(i+1, j+1));
		Point BL(voxel2Pos(i,0),   voxel2Pos(j+1,1), getZ(i,   j+1));
		Point pts[4] = {
				(TL), (TR),
				(BR), (BL)};

		hitSomething = IntersectHelper(rayW2O, pts, i, j, &isect);

		if (!bounds.IntersectP(rayW2O)) return false;

		// in heightfields, there will be no overlap
		if (hitSomething) break;

        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
                   ((NextCrossingT[0] < NextCrossingT[2]) << 1) +
                   ((NextCrossingT[1] < NextCrossingT[2]));
		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
        int stepAxis = cmpToAxis[bits];
		if (rayW2O.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }

	if (!hitSomething) return false;

	*tHit = tHit1;
	*dg = isect.dg;
	*rayEpsilon = isect.rayEpsilon;


	return true;
}

bool Heightfield2::IntersectHelper(const Ray &ray, const Point *pts, int i, int j, Intersection *in) const {


	BBox bounds(Union(BBox(pts[0], pts[1]), BBox(pts[2], pts[3])));
	if (!bounds.IntersectP(ray)) return false;


	int vptr[6] = {0,1,2,0,2,3};

	Normal normals[4] = {
		vertexNormals[j*nx + i],
		vertexNormals[j*nx + i + 1],
		vertexNormals[(j+1)*nx + i + 1],
		vertexNormals[(j+1)*nx + i]};
	float uvs[8] = {
		pts[0].x, pts[0].y,
		pts[1].x, pts[1].y,
		pts[2].x, pts[2].y,
		pts[3].x, pts[3].y};

	TriangleMesh *triMesh = new TriangleMesh(ObjectToWorld, WorldToObject, ReverseOrientation,
		2, 4, vptr, pts, normals, NULL, uvs, NULL);
	Triangle *triangle = new Triangle(ObjectToWorld, WorldToObject, ReverseOrientation, triMesh, 0);

	Triangle *triangle2 = new Triangle(ObjectToWorld, WorldToObject, ReverseOrientation, triMesh, 1);


	Intersection in1, in2;
	bool tri1 = triangle->Intersect((*ObjectToWorld)(ray), &(tHit1), &(in1.rayEpsilon), &(in1.dg));
	bool tri2 = triangle2->Intersect((*ObjectToWorld)(ray), &(tHit2), &(in2.rayEpsilon), &(in2.dg));

	if (!tri1 && !tri2) return false;

	if (tri1 && tri2) {
		if (tHit1 < tHit2) {
			tHit2 = tHit1;
			*in = in1;
			tri2 = false;
		} else {
			tHit1 = tHit2;
			*in = in2;
			tri1 = false;
		}
	} else if (tri1) {
		*in = in1;
		tHit2 = tHit1;
	} else {
		*in = in2;
		tHit1 = tHit2;
	}

	return true;
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


//void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
//    int ntris = 2*(nx-1)*(ny-1);
//    refined.reserve(ntris);
//    int *verts = new int[3*ntris];
//    Point *P = new Point[nx*ny];
//    float *uvs = new float[2*nx*ny];
//    int nverts = nx*ny;
//    int x, y;
//    // Compute heightfield2 vertex positions
//    int pos = 0;
//    for (y = 0; y < ny; ++y) {
//        for (x = 0; x < nx; ++x) {
//            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
//            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
//            P[pos].z = z[pos];
//            ++pos;
//        }
//    }
//
//    // Fill in heightfield2 vertex offset array
//    int *vp = verts;
//    for (y = 0; y < ny-1; ++y) {
//        for (x = 0; x < nx-1; ++x) {
//#define VERT(x,y) ((x)+(y)*nx)
//            *vp++ = VERT(x, y);
//            *vp++ = VERT(x+1, y);
//            *vp++ = VERT(x+1, y+1);
//
//            *vp++ = VERT(x, y);
//            *vp++ = VERT(x+1, y+1);
//            *vp++ = VERT(x, y+1);
//        }
//#undef VERT
//    }
//    ParamSet paramSet;
//    paramSet.AddInt("indices", verts, 3*ntris);
//    paramSet.AddFloat("uv", uvs, 2 * nverts);
//    paramSet.AddPoint("P", P, nverts);
//    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
//    delete[] P;
//    delete[] uvs;
//    delete[] verts;
//}


