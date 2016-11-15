
#ifndef BOUNDARYCONSTRAINT_H
#define BOUNDARYCONSTRAINT_H

#include "main.h"
#include "vec.h"
#include "GridND.h"
#include "DelaunayTri.h"

typedef struct BoundarySample
{
	Vec3f m_samplePos;
	Vec3f m_sampleNorm;
} BoundarySample;

class CBoundaryConstraint
{
public:
	CBoundaryConstraint();

	bool ApplyBoundaryConstraint(Vec3f& posOld, Vec3f& posNew);

	bool ApplyBoundaryConstraintNew(Vec3f& posOld, Vec3f& posNew);

	bool InsideBoundary(Vec3f& pos);

	bool InsideBoundaryNew(Vec3f& pos);

	bool InsideBoundary(int px, int py);

	void LoadBoundaryConstraintFromImage(string fileName);

	CGrid2D3f LoadGrid2D3fFromCImage(const char* fileName);

	void SaveGrid2D3fAsCImage(const CGrid2D3f& grid, const char* fileName);

	static float m_colorThresh;

	inline void SetPosMin(Vec3f posMin) { m_bcPosMin = posMin; }
	inline void SetPosMax(Vec3f posMax) { m_bcPosMax = posMax; }
	inline Vec3f GetPosMin() { return m_bcPosMin; }
	inline Vec3f GetPosMax() { return m_bcPosMax; }

	inline void SetFrameIdx(int idx) { m_bcFrameIdx = idx; }
	inline int GetFrameIdx() { return m_bcFrameIdx; }

	inline void SetInputIdx(int idx) { m_bcInputIdx = idx; }
	inline int GetInputIdx() { return m_bcInputIdx; }

private:
	CGrid2D3f IdentifyGridContour(const CGrid2D3f& gridInput);

	void Triangulate();

	Vec3f ComputeSampleNorm(int px, int py);

	CGrid2D3f m_grid;
	CGrid2D3f m_gridContour;
	vector<BoundarySample> m_vecBoundarySample;
	Vec3f m_bcPosMin;
	Vec3f m_bcPosMax;
	int m_bcFrameIdx;
	int m_bcInputIdx;
};

#endif // BOUNDARYCONSTRAINT_H
