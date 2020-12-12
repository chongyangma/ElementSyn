#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include "ParticleData.h"

class CParticleSystem
{
public:
    CParticleSystem();

    bool LoadParticleSystemFromTXT(const string& fileName);

    bool SaveParticleSystemAsTXT(const string& fileName);

    bool LoadParticleSystemFromCSV(const string& fileName);

    bool SaveParticleSystemAsCSV(const string& fileName);

    bool LoadParticleSystem(const string& fileName);

    bool SaveParticleSystem(const string& fileName);

    inline int GetNumOfSoftBodies() { return int(m_vecParticleData.size()); }

    inline int GetNumOfSamples()
    {
        if (m_vecParticleData.empty() == true)
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

    inline void AddParticle(CParticleData& particle) { m_vecParticleData.push_back(particle); }

    void SetNeighboringSamples(Flt neighDist);

    void ScaleParticleSystem(Flt scaling);

private:
    vector<CParticleData> m_vecParticleData;
};

#endif // PARTICLESYSTEM_H
