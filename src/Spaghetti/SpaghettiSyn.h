
#ifndef SPAGHETTISYN_H
#define SPAGHETTISYN_H

//#define ONLY_COARSE_MOTION

#include "../DynamicElement/MathTypes.h"
#include "../DynamicElement/MathFunc.h"
using namespace machy_math;
#include "SpaghettiSynConfig.h"
#include "BulletSpaghetti.h"

class CSpaghettiSyn
{
public:
	CSpaghettiSyn();

	~CSpaghettiSyn();

	void LoadInputData();

	void UpdateOutput();

	void RenderOutput();

	void RestartSynthesis();

	inline int GetStepCount() { return m_stepCount; }

private:
	void ResetOutput(string configFileName);

	void GatherInputNeighborsExtended();

	Flt GetNearestInputNeighborExtended(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices);

	Flt NeighborhoodMetricExtended(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices);

	// Samples from different strands only
	Flt GetNearestInputNeighborExtendedNew(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices);

	Flt NeighborhoodMetricExtendedNew(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices);

	void InitializeOutput();

	void InitializeOutputViaPoissonDiskVertically();

	void InitializeOutputViaPoissonDiskHorizontally();

	void UpdateOutputViaOptimization();

	void UpdateFrameViaOptimization();

	void UpdateOneFrameViaOptimization(int frameIdx, int endFrameIdx);

	void UpdateFrameViaBulletSimulation();

	void ResetCoeffMatrix();

	void ResetFrameCoeffMatrix();

	void CollisionResponse(vector<Flt>& vecCx, vector<Flt>& vecCy, vector<Flt>& vecCz);

	MassNeighborhoodExtended GetMassNeighborhoodExtended(vector<CMassSpringData>& strands, int strandIdx, int massIdx, Flt neighDist);

	MassNeighborhoodExtended GetMassNeighborhoodExtended1(int frameIdx, int strandIdx, int massIdx, Flt neighDist, bool flagToroidal = false); // For input

	MassNeighborhoodExtended GetMassNeighborhoodExtended2(int frameIdx, int strandIdx, int massIdx, Flt neighDist); // For output

	MassNeighborhoodExtended GetMassNeighborhoodExtended2new(int strandIdx, int massIdx, Flt neighDist); // For output, frame-by-frame version

	vector<SegmentPatchExtended> GetSegmentPatchesExtended(vector<MassSpringSequence>& vecSequence);

	Flt SegmentPatchSpatialCompatibility(SegmentPatchExtended* ptrPatch1, SegmentPatchExtended* ptrPatch2);

	Flt SegmentPatchTemporalCompatibility(SegmentPatchExtended* ptrPatch1, SegmentPatchExtended* ptrPatch2);

	vector<vector<Vec3f>> GetConnectedInputPatchExtend(int segmentLength, int sequenceLength);

	void TimeIntegration();

	void CollisionResponse();

	void FinishCoeffsAccordingToBoundary();

	void ApplyGroundBoundary();

	void UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos, bool flagUpdateMat = true);

	void UpdateCoeffMatPairVals(int idxi, int idxj, Flt wt, Vec3f pr, bool flagUpdateMat = true);

	void ChangeInputSequence();

	void ChangeInputSequenceAccordingToCoarseMotion(Flt segmentLength);

	void ChangeCoarseSequence0(); // 12/29/2012: Fix the bug near the end

	void ChangeCoarseSequence1(int frameIdx = 63, int wd = 50, int segmentLength = 10, Vec3f pd = Vec3f(0.1f, 0.0f, 0.0f), int type = 0);

	void ChangeCoarseSequence();

	void ChangeCoarseSequenceAsNonStretchable();

	void DumpChopsticks(vector<MassSpringSequence>& coarseSequence);

	void DumpChopsticks(MassSpringSequence& coarseStrands, const char* fileName);

	vector<CMassSpringData> EditOneOutputFrame(vector<CMassSpringData>& strands);

	vector<CMassSpringData> EditOneOutputFrame1(vector<CMassSpringData>& strands, Flt wt = 1.0f);

public:
	vector<CMassSpringData> EditCoarseStrandsFromSynthesisOutput(vector<CMassSpringData>& vecCoarseStrand);

	vector<CMassSpringData> ResolveStrandsCollision(vector<CMassSpringData>& strands, Flt rad);

	vector<CMassSpringData> ResolveStrandsCollisionNew(vector<CMassSpringData>& strands, Flt rad);

#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
	void ResetInputNeighborUsedCount();
#endif

	inline Flt VecPrDifference(vector<Vec3f>& vecPr1, vector<Vec3f>& vecPr2)
	{
		if ( vecPr2.size() < vecPr1.size() ) return 1e10;
		Flt diff = 0.0f;
		//int sizeDiff = int(vecPr2.size() - vecPr1.size());
		for ( int i=0; i<int(vecPr1.size()); i++ )
		{
			diff += dist2(vecPr1[i], vecPr2[i]) * exp(CSpaghettiSynConfig::m_spatialSigma * i * i);
		}
		return diff;
	}

	inline Flt VecPrDifferenceTemporal(vector<Vec3f>& vecPr1, vector<Vec3f>& vecPr2)
	{
		if ( vecPr2.size() < vecPr1.size() ) return 1e10;
		Flt diff = 0.0f;
		//int sizeDiff = int(vecPr2.size() - vecPr1.size());
		for ( int i=0; i<int(vecPr1.size()); i++ )
		{
			Flt i2 = (i + 1) * (i + 1);
			diff += dist2(vecPr1[i], vecPr2[i]) * exp(CSpaghettiSynConfig::m_temporalSigma * i2);
		}
		return diff;
	}

	inline string GetLogFileName()
	{
		char fileName[MAX_PATH];
		sprintf_s(fileName, "%s\\log.txt", CSpaghettiSynConfig::m_outputPrefix.c_str());
		return string(fileName);
	}

	int m_stepCount;
	int m_serialCount;
	CSpaghettiSynConfig* m_ptrSynConfig;
	CCrossList* m_ptrCoeffMatrix;
	vector<Flt> m_vecCx, m_vecCy, m_vecCz;
	vector<MassNeighborhoodExtended> m_vecInputNeighExtended;
	vector<SegmentPatchExtended> m_vecInputPatchExtended;
	vector<MassSpringSequence> m_vecInputSequence;
	vector<CMassSpringData> m_vecOutputStrand;
	vector<MassSpringSequence> m_vecOutputSequence;
	vector<CMassSpringData> m_vecInitialStrand;
	vector<MassSpringSequence> m_vecInitialSequence;
	vector<MassSpringSequence> m_vecCoarseSequence;
};

#endif SPAGHETTISYN_H
