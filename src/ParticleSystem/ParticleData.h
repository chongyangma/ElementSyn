
#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

#include "ParticleSystemConfig.h"
#include "../DynamicElement/MathFunc.h"
using namespace machy_math;

class CNeighboringSample
{
public:
	CNeighboringSample(Vec3f pr, int sampleIdx, int neighBodyIdx, int neighSampleIdx)
		: m_neighPr(pr), m_sampleIdx(sampleIdx), m_neighBodyIdx(neighBodyIdx), m_neighSampleIdx(neighSampleIdx)
	{}

	Vec3f m_neighPr;
	int m_sampleIdx;
	int m_neighBodyIdx;
	int m_neighSampleIdx;
};

class CParticleData
{
public:
	CParticleData();

	inline void SetPos(Vec3f& pos) { m_pos = pos; }
	inline Vec3f& GetPos() { return m_pos; }

	inline void SetVel(Vec3f& vel) { m_vel = vel; }
	inline Vec3f& GetVel() { return m_vel; }

	inline void SetQuat(Vec4f& quat) { m_quat = quat; }
	inline Vec4f& GetQuat() { return m_quat; }

	inline void SetVecSamplePos(vector<Vec3f> vecSamplePos) { m_vecSamplePos = vecSamplePos; }
	inline vector<Vec3f>& GetVecSamplePos() { return m_vecSamplePos; }

	inline void SetFlagBoundary(bool flagBoundary) { m_flagBoundary = flagBoundary; }
	inline bool GetFlagBoundary() { return m_flagBoundary; }

	inline void SetFlagFixed(bool flagFixed) { m_flagFixed = flagFixed; }
	inline bool GetFlagFixed() { return m_flagFixed; }

	inline int GetNumOfSamples() { return int(m_vecSamplePos.size()); }

	void RenderSoftBody(Vec3f trans = Vec3f(0.0f, 0.0f, 0.0f));

	void TranslateSoftBody(Vec3f trans);

	inline void ClearNeighboringSamples() { m_vecNeighboringSample.clear(); }
	inline void AddNeighboringSample(CNeighboringSample& neigh) { m_vecNeighboringSample.push_back(neigh); }
	inline vector<CNeighboringSample>& GetNeighboringSamples() { return m_vecNeighboringSample; }

	inline void ScaleParticleData(Flt scaling) { m_pos *= scaling; }

private:
	Vec3f m_pos;
	Vec3f m_vel;
	Vec4f m_quat;
	vector<Vec3f> m_vecSamplePos;
	vector<CNeighboringSample> m_vecNeighboringSample;
	bool m_flagBoundary;
	bool m_flagFixed;
};

#endif // PARTICLEDATA_H
