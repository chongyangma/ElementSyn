
#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include "ParticleData.h"

class CParticleSystem
{
public:
	CParticleSystem();

	bool LoadParticleSystemFromTXT(string fileName);

	bool DumpParticleSystemToTXT(string fileName);

	inline int GetNumOfSoftBodies() { return int(m_vecParticleData.size()); }

	inline int GetNumOfSamples()
	{
		if ( m_vecParticleData.empty() == true )
		{
			return 0;
		}
		else
		{
			return int(m_vecParticleData.size()) * m_vecParticleData[0].GetNumOfSamples();
		}
	}

	inline CParticleData& GetParticleData(int idx) { return m_vecParticleData[idx]; }

	inline void ClearParticleSystem() { m_vecParticleData.clear(); }

	inline void ResizeParticleSystem(int sz) { m_vecParticleData.resize(sz); }

	inline void AddSoftBody(CParticleData& softBody) { m_vecParticleData.push_back(softBody); }

	void SetNeighboringSamples(Flt neighDist);

	void ScaleParticleSystem(Flt scaling);

private:
	vector<CParticleData> m_vecParticleData;
};

#endif // PARTICLESYSTEM_H
