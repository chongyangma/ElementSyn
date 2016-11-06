
#ifndef	DELAUNAYTRI_H
#define DELAUNAYTRI_H

// Implement of Delaunay triangulation

#include "main.h"
#include "vec.h"

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif


// TYPES //////////////////////////////////////////////////
typedef struct VERTEX2D_TYP
{
	REAL x;
	REAL y;

} VERTEX2D, *VERTEX2D_PTR;

typedef struct TRIANGLE_TYP
{
	int i1; // vertex index
	int i2; 
	int i3; 

	TRIANGLE_TYP* pNext;
	TRIANGLE_TYP* pPrev;

} TRIANGLE, *TRIANGLE_PTR;

typedef struct MESH_TYP
{
	int vertex_num;
	int triangle_num;

	VERTEX2D_PTR pVerArr; // point to outer vertices array
	TRIANGLE_PTR pTriArr; // point to outer triangles array

} MESH, *MESH_PTR;

class CDelaunayTri
{
public:
	CDelaunayTri();

	~CDelaunayTri();

	void InitMesh(vector<Vec2f> vecPoint);

	REAL CounterClockWise(VERTEX2D_PTR pa, VERTEX2D_PTR pb, VERTEX2D_PTR pc);

	REAL CounterClockWise(VERTEX2D va, VERTEX2D vb, VERTEX2D vc);

	TRIANGLE_PTR AddTriangleNode(TRIANGLE_PTR pPrevTri, int i1, int i2, int i3);

	void AddBoundingBox();

	void RemoveBoundingBox();

	void RemoveTriangleNode(TRIANGLE_PTR pTri);

	REAL InCircle(VERTEX2D_PTR pa, VERTEX2D_PTR pb, VERTEX2D_PTR pp, VERTEX2D_PTR  pd);

	REAL InTriangle(VERTEX2D_PTR pVer, TRIANGLE_PTR pTri);

	bool FlipTest(TRIANGLE_PTR pTestTri);

	void InsertOnEdge(TRIANGLE_PTR pTargetTri, int ver_index);

	void InsertInTriangle(TRIANGLE_PTR pTargetTri, int ver_index);

	void Insert(int ver_index);

	void IncrementalDelaunay();

	MESH_PTR GetMeshPtr() { return m_pMesh; }

private:
	MESH_PTR m_pMesh;

};

#endif DELAUNAYTRI_H
