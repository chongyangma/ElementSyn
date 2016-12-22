
#ifndef MASSSPRINGSYNCONFIG_H
#define MASSSPRINGSYNCONFIG_H

#include "../DynamicElement/SynConfigBase.h"
#include "../DynamicElement/MassSpringData.h"
#ifdef _DEBUG
#pragma comment(lib, "../debug/DynamicElement.lib")
#else
#pragma comment(lib, "../release/DynamicElement.lib")
#endif

class CMassSpringSynConfig : public CSynConfigBase
{
public:
	CMassSpringSynConfig();

	bool ReloadConfigFromFile(string fileName);

	bool LoadFromMassSpringSynConfig(string fileName = "MassSpringSyn_config.txt");

	static bool m_flagTemporallyToroidal;
	static int m_numOfInputFrames;
	static int m_numOfMasses;
	static int m_numOfStrands;
	static int m_numOfInputStrands;
	static int m_neighSize;
	static int m_patchSize;
	static int m_inputLoadInterval;
	static int m_inputLoadStart;
	static int m_stepCountMax;
	static int m_synthesisWindowSize;
	static string m_inputPrefix;
	static string m_outputPrefix;
	static string m_exemplarName;
	static vector<string> m_vecInputPrefix;
	static vector<Flt> m_vecInputWt;
	static Flt m_springLength;
	static Flt m_wtIdxDiff;
	static Flt m_wtPrDiff;
	static Flt m_wtNeighPrDiff;
	static Flt m_shapeWt;
	static Flt m_smoothWt;
	static Flt m_alignWt;
	static Flt m_temporalWt;
	static Flt m_inputWtSumInv;
	static Flt m_strandRadius;
	static Flt m_strandRepulsionConstant;
	static Flt m_neighDist;
	static Flt m_gaussianSigma;
	static Flt m_temporalSigma;
	static Vec2f m_poisson2Dmin;
	static Vec2f m_poisson2Dmax;

private:
	void DumpParameters(FILE* file);
};

#endif MASSSPRINGSYNCONFIG_H
