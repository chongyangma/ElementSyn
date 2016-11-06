
#ifndef SPAGHETTISYNCONFIG_H
#define SPAGHETTISYNCONFIG_H

#include "../DynamicElement/SynConfigBase.h"
#include "../DynamicElement/MassSpringData.h"
#ifdef _DEBUG
#pragma comment(lib, "../debug/DynamicElement.lib")
#else
#pragma comment(lib, "../release/DynamicElement.lib")
#endif

class CSpaghettiSynConfig : public CSynConfigBase
{
public:
	CSpaghettiSynConfig();

	bool ReloadConfigFromFile(string fileName);

	bool LoadFromSpaghettiSynConfig(string fileName = "SpaghettiSyn_config.txt");

	bool DumpToSpaghettiSynConfig();

	static bool m_flagTemporallyToroidal;
	static int m_numOfMasses;
	static int m_spatialNeighSize; // samples within the same strand
	static int m_temporalWindowSize;
	static int m_iterationNum;
	static int m_inputLoadInterval;
	static int m_inputLoadStart;
	static int m_inputLoadEnd;
	static int m_coarseLoadInterval;
	static int m_coarseLoadStart;
	static int m_coarseLoadEnd;
	static int m_interleavedSteps;
	static Flt m_springLength;
	static Flt m_springRadius;
	// Parameters for simulation...
	static Flt m_springStiffness;
	static Flt m_springFriction;
	static Flt m_springRepulsion;
	static Flt m_gravity;
	static Flt m_airFriction;
	static Flt m_groundHeight;
	static Flt m_groundRepulsion;
	static Flt m_groundFriction; // ground slide friction constant
	static Flt m_groundAbsorption;
	// Parameters for synthesis...
	static Flt m_spatialWt;
	static Flt m_spatialSigma;
	static Flt m_spatialNeighDist; // samples from different strands
	static Flt m_temporalWt;
	static Flt m_temporalSigma;
	static Flt m_compatibilityThresh;
	static Vec2f m_poisson2Dmin;
	static Vec2f m_poisson2Dmax;

private:
	void ResetConfig();

	void DumpParameters(FILE* file);
};

#endif SPAGHETTISYNCONFIG_H
