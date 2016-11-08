
#ifndef PARTICLESYSTEMSYN_H
#define PARTICLESYSTEMSYN_H

#include "ParticleSystem.h"
#include "../DynamicElement/HungarianAlgorithm.h"
#include "../DynamicElement/TAUCSsolver.h"
#include "../DynamicElement/BoundaryConstraint.h"

// For initialization by patch copy
typedef struct CubicPatch
{
	Vec3f m_patchCen;
	vector<int> m_patchIndices; // indices of elements in the arrangement
	int m_patchFrequency; // times appear in the destination
	int m_patchFrameIdx;
} CubicPatch;

class CSampleMatchingResult
{
public:
	CSampleMatchingResult()
	{
		m_frameIdx = -1;
		m_sampleIdx = -1;
		m_distSum = 1e10;
		m_weight = 1.0f;
	}

	int m_frameIdx;
	int m_sampleIdx;
	vector<int> m_matchingIndices;
	Flt m_distSum;
	Flt m_weight;
};

class CParticleSystemSyn
{
public:
	CParticleSystemSyn(const std::string& config_file_name);

	~CParticleSystemSyn();

	void LoadInputData();

	void LoadInputDataNew(); // 01/07/2013

	void UpdateOutput();

	void RenderOutput();

	void RestartSynthesis();

	inline int GetStepCount() { return m_stepCount; }

	static GLUquadricObj* m_ptrQuadricObj;

private:
	void SetInputNeighborhoods();

	void ResetOutput(string configFileName);

	void InitializeOutput();

	void InitializeOutputTrajectoriesViaPatchCopy(int numOfParticles = -1, int boundaryType = 0);

	void InitializeOutputTrajectoriesViaPatchCopy2(int numOfParticles = -1);

	void InitializeOutputTrajectoriesViaPatchCopy3(int numOfParticles = -1, CBoundaryConstraint* ptrBoundaryConstraint = NULL);

	void AdvectOutput();

	void AdvectOutputNew(); // 01/07/2013

	void UpdateOutputGroupViaOptimization(CParticleSystem& outputGroup, int boundaryType = 0);

	void UpdateOutputGroupViaOptimization3(CParticleSystem& outputGroup, CBoundaryConstraint* ptrBoundaryConstraint);

	void ResetCoeffMatrix(CParticleSystem& outputGroup);

	void FinishCoeffsAccordingToBoundary();

	void UpdateOutputFrameViaBlending(int frameIdx, CBoundaryConstraint* ptrBoundaryConstraint1, CBoundaryConstraint* ptrBoundaryConstraint2);

	bool ApplyBoundaryCondition(Vec3f& posOld, Vec3f& posNew);

	bool ApplyBoundaryCondition2(Vec3f& posOld, Vec3f& posNew);

	Flt GetNearestInputSoftBodies(CParticleData& softBodyOut, int outputIdx, CSampleMatchingResult& matchingRst1, CSampleMatchingResult& matchingRst2); // No temporal window

	Flt GetNearestInputSoftBodies(CParticleData& softBodyOut, int inputIdx, int outputIdx, CSampleMatchingResult& matchingRst1, CSampleMatchingResult& matchingRst2); // No temporal window

	Flt GetNearestInputSoftBodiesMorph(CParticleData& softBodyOut, vector<CSampleMatchingResult>& vecMatchingRst);

	Flt SoftBodyNeighDistance(CParticleData& softBody1, CParticleData& softBody2, vector<int>& matchIndices, Flt weight);

	void UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos);

	void UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos, int dim);

	void UpdateCoeffMatPairVals(int idxi, int idxj, Flt wt, Vec3f pr);

	CParticleSystem LoadInputExemplarFromTXT(string fileName);

	vector<CubicPatch> GetCubicPatches(vector<CParticleSystem>& groupSequence);

	vector<CubicPatch> GetCubicPatches(CParticleSystem& softBodyGroup);

	vector<int> GetCubicPatchIndices(CParticleSystem& softBodyGroup, Vec3f patchCen, Vec3f patchSize);

	vector<CParticleSystem> InterpolateParticleSystem(CParticleSystem& softBodyGroup1, CParticleSystem& softBodyGroup2, int sequenceLength);

	void LoadBoundaryConstraints();

	void SetSequencesVelocity(vector<CParticleSystem>& vecSequence, bool flagTemporallyToroidal = false);

	void SampleRepulsion();

	int m_stepCount;
	int m_serialCount;
	CParticleSystemConfig* m_ptrSynConfig;
	vector<CDenseMatrix*> m_vecPtrCoeffMatrix;
	vector<Flt> m_vecCx, m_vecCy, m_vecCz;
	vector<CParticleSystem> m_inputSequence;
	vector<CParticleSystem> m_vecInputExemplar;
	CParticleSystem m_outputGroup;
	CParticleSystem m_finalGroup; // target state
	vector<CParticleSystem> m_outputSequence;
	vector<Vec3f> m_vecGoalPos;
	vector<CBoundaryConstraint> m_vecBoundaryConstraint;
	vector<Vec3f> m_diffuseColorSequence;
};

#endif PARTICLESYSTEMSYN_H
