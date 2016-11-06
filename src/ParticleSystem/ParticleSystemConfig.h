
#ifndef PARTICLESYSTEMCONFIG_H
#define PARTICLESYSTEMCONFIG_H

#include "../DynamicElement/SynConfigBase.h"
#ifdef _DEBUG
#pragma comment(lib, "../debug/DynamicElement.lib")
#else
#pragma comment(lib, "../release/DynamicElement.lib")
#endif

class CParticleSystemConfig : public CSynConfigBase
{
public:
	CParticleSystemConfig();

	bool ReloadConfigFromFile(string fileName);

	bool LoadFromParticleSystemConfig(string fileName = "ParticleSystemSyn_config.txt");

	static bool m_flagHybridSolver;
	static bool m_flagSmoothSynthesis;
	static int m_initMethod;
	static int m_minimumPatchSize;
	static int m_iterationNum;
	static int m_temporalWindowSize;
	static string m_trajectoriesFileName;
	static vector<string> m_vecInputFileName;
	static Flt m_shapeWt;
	static Flt m_neighWt;
	static Flt m_spatialWt; // weight of additional exemplar for spatial arrangement
	static Flt m_boundaryWt;
	static Flt m_temporalWt;
	static Flt m_neighDist;
	static Flt m_repulsionDist;
	static Flt m_repulsionWt;
	static Flt m_compatibilityThresh;
	static Flt m_boundaryThresh;
	static Flt m_gaussianSigma;
	static Vec3f m_inputTrans;
	static Vec3f m_outputTrans;
	static Vec3f m_cubicPosMin;
	static Vec3f m_cubicPosMax;
	static Vec3f m_patchSize;
	static Vec3f m_particleAmbient;
	static Vec3f m_particleDiffuse;
	static Vec3f m_particleSpecular;
	static vector<Vec3f> m_vecDiffuseColor;

private:
	void DumpParameters(FILE* file);
};

#endif PARTICLESYSTEMCONFIG_H
