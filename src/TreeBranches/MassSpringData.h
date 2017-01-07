
#ifndef MASSSPRINGDATA_H
#define MASSSPRINGDATA_H
//#define HACK_WITH_RELATIVE_POSITION
#define HACK_WITH_NEIGHBOR_USER_COUNT

#include "MassSpringSynConfig.h"
#include "../DynamicElement/NURBSinterpolate.h"
#include "../DynamicElement/RenderFunc.h"

class CMass
{
public:
	CMass();

	inline int GetSrcIdx() { return m_srcIdx; }
	inline void SetSrcIdx(int srcIdx) { m_srcIdx = srcIdx; }

	inline Vec3f& GetPos() { return m_pos; }
	inline Vec3f GetPos() const { return m_pos; }
	inline void SetPos(Vec3f pos) { m_pos = pos; }

	inline Vec3f& GetVel() { return m_vel; }
	inline void SetVel(Vec3f vel) { m_vel = vel; }

	inline Vec3f& GetForce() { return m_force; }
	inline void SetForce(Vec3f force) { m_force = force; }
	inline void ResetForce() { SetZero(m_force); }
	inline void AddForce(Vec3f dForce) { m_force += dForce; }

	inline Vec4f& GetQuat() { return m_quat; }
	inline void SetQuat(Vec4f quat) { m_quat = quat; }

	inline bool GetFlagFixed() { return m_flagFixed; }
	inline void SetFlagFixed(bool flagFixed) { m_flagFixed = flagFixed; }

	void TimeIntegrationMass(const Flt dt);

private:
	int m_srcIdx;
	Vec3f m_pos;
	Vec3f m_vel;
	Vec3f m_force;
	Vec4f m_quat;
	bool m_flagFixed;
};

class CMassSpringData
{
public:
	CMassSpringData();

	void ComputeForce(const float springLength, const float stiffness, const float springFriction, const float gravity, const float airFriction);

	void ComputeSpringForce(const int idx, const float springLength, const float stiffness, const float springFriction);

	void TimeIntegrationMassSpring(const Flt dt);

	inline void ResizeMassSpringData(int sz) { m_vecMass.resize(sz); }

	inline int GetNumOfMasses() const { return int(m_vecMass.size()); }

	inline CMass& GetMass(int idx) { return m_vecMass[idx]; }
	inline const CMass& GetMass(int idx) const { return m_vecMass[idx]; }
	inline void AddMass(CMass& mass) { m_vecMass.push_back(mass); }

	inline int GetSrcIdx() { return m_srcIdx; }
	inline void SetSrcIdx(int srcIdx) { m_srcIdx = srcIdx; }

	inline int GetStartIdx() { return m_startIdx; }
	inline void SetStartIdx(int startIdx) { m_startIdx = startIdx; }

	void RenderMassSpringData(Vec3f trans = Vec3f(0.f, 0.f, 0.f));

	CMassSpringData NurbsInterpolateStrand(int interval = 5);

	CMassSpringData DownsampleStrand(int interval = 5);

	CMassSpringData TranslateStrand(Vec3f trans);

	CMassSpringData DeformStrand(Flt d);

	CMassSpringData FixLengthAsNonStretchable();

	static vector<CMassSpringData> LoadStrandsFromTXT(string fileName, int maxNumOfStrands = -1);

	static vector<CMassSpringData> NurbsInterpolateStrands(vector<CMassSpringData>& strands, int interval = 5);

	static vector<CMassSpringData> DownsampleStrands(vector<CMassSpringData>& strands, int interval = 5);

	static bool DumpStrandsToTXT(string fileName, vector<CMassSpringData>& strands, const Vec3f trans = Vec3f(0.0f, 0.0f, 0.0f));

	static bool DumpStrandsToOBJ(vector<CMassSpringData>& strands, string fileName = "lines.obj", Flt rad = 0.01f, int sampleNum = 8, int interval = 1, int interp = 0);

	static void DumpStrandToOBJ(CMassSpringData& strand, ofstream& fout, int& vAccum, Flt rad, int sampleNum, int interval);

	static bool DumpStrandsToOBJz(vector<CMassSpringData>& strands, string fileName = "lines.obj", Flt rad = 0.01f, int sampleNum = 8, int interval = 1, int interp = 0);

	static void DumpStrandToOBJz(CMassSpringData& strand, ofstream& fout, int& vAccum, Flt rad, int sampleNum, int interval);

private:
	int m_srcIdx;
	int m_startIdx;
	vector<CMass> m_vecMass;
};

typedef vector<CMassSpringData> MassSpringSequence;

typedef struct MassNeighborhoodExtended
{
	int m_massIdx;
#ifdef HACK_WITH_RELATIVE_POSITION
	Flt m_relativePos;
#endif
#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
	int m_neighUsedCount;
#endif
	Vec3f m_massPos; // Center position
	vector<Vec3f> m_vecPr0; // Temporal trajectory
	vector<Vec3f> m_vecNeighPr1;
	vector<Vec3f> m_vecNeighPr2;
	vector<Vec3f> m_vecNeighPr3;
	vector<vector<Vec3f> > m_vecNeighPr1Prev;
	vector<vector<Vec3f> > m_vecNeighPr2Prev;
	vector<vector<Vec3f> > m_vecNeighPr3Prev;
	vector<int> m_vecNeighIndex;
} MassNeighborhoodExtended;

typedef struct SegmentPatch
{
	vector<Vec3f> m_vecPos;
	vector<Vec3f> m_vecVel;
	bool m_flagStart;
} SegmentPatch;

// Spatial-temporal patch...
typedef struct SegmentPatchExtended
{
	vector<vector<Vec3f> > m_vecPosSequence;
	vector<vector<Vec3f> > m_vecVelSequence;
	bool m_flagStart;
} SegmentPatchExtended;

#endif // MASSSPRINGDATA_H
