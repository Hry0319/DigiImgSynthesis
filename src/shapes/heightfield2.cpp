
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

	InitNormals();

}


Heightfield2::~Heightfield2() {
    delete[] z;
}


//
// Heightfield2::InitNormals()
//
void Heightfield2::InitNormals() {
	Point          p[9];
//	Vector         vector[8];


	int voxel2posX = width-1;
	int voxel2posY = height-1;
  	int voxel2posZ = 1;

	normal = new Normal[nx*ny];

	for ( int j = 0; j < ny; j++)
	{
		for ( int i = 0; i < nx; i++)
		{
			p[4] = Point(i/(float)voxel2posX, j/(float)voxel2posY, z[i+j*nx]/(float)voxel2posZ);

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



			normal[i+j*nx] = Normal( Normalize(	Cross( (p[0]-p[4]) ,(p[1]-p[4]) ) + \
												Cross( (p[1]-p[4]) ,(p[2]-p[4]) ) + \
												Cross( (p[2]-p[4]) ,(p[5]-p[4]) ) + \
												Cross( (p[5]-p[4]) ,(p[8]-p[4]) ) + \
												Cross( (p[8]-p[4]) ,(p[7]-p[4]) ) + \
												Cross( (p[7]-p[4]) ,(p[6]-p[4]) ) + \
												Cross( (p[6]-p[4]) ,(p[3]-p[4]) ) + \
												Cross( (p[3]-p[4]) ,(p[0]-p[4]) )
												))*(-1);

//			normal[i+j*nx] = Normal( Normalize(	Cross( (p[0]-p[4]) ,(p[1]-p[4]) ) + \
//												Cross( (p[1]-p[4]) ,(p[5]-p[4]) ) + \
//												Cross( (p[5]-p[4]) ,(p[8]-p[4]) ) + \
//												Cross( (p[8]-p[4]) ,(p[7]-p[4]) ) + \
//												Cross( (p[7]-p[4]) ,(p[3]-p[4]) ) + \
//												Cross( (p[3]-p[4]) ,(p[0]-p[4]) )
//												))*(-1);
			printf("%e %e %e  \n" , normal[i+j*nx].x, normal[i+j*nx].y, normal[i+j*nx].z);

		}
	}
}

bool Heightfield2::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const{
	float rayT;
	BBox bounds = ObjectBound();
    if (bounds.Inside(ray(ray.mint)))
        return false;
    else if (!bounds.IntersectP(ray, &rayT))
        return false;

	return false;
}
bool Heightfield2::IntersectP(const Ray &ray) const{

	return false;
}

BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
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


