#include "ParticleSystem.h"

CParticleSystem::CParticleSystem()
{
}

bool CParticleSystem::LoadParticleSystemFromTXT(const string& fileName)
{
    ifstream fin(fileName.c_str());
    if (fin.fail() == true)
    {
        cout << "Failed to load particle system from " << fileName << "!\n";
        return false;
    }
    int numOfSoftBodies;
    fin >> numOfSoftBodies;
    m_vecParticleData.resize(numOfSoftBodies);
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        int numOfVertices;
        fin >> numOfVertices;
        vector<Vec3f> vertices(numOfVertices);
        for (int j = 0; j < numOfVertices; j++)
        {
            Vec3f pos;
            fin >> pos[0] >> pos[1] >> pos[2];
            vertices[j] = pos;
        }
        CParticleData& softBodyData = m_vecParticleData[i];
        softBodyData.SetVecSamplePos(vertices);
        softBodyData.SetPos(vertices[0]);
    }
    return true;
}

bool CParticleSystem::SaveParticleSystemAsTXT(const string& fileName)
{
    ofstream fout(fileName.c_str());
    if (fout.fail() == true)
    {
        cout << "Failed to save particle system into " << fileName << "!\n";
        return false;
    }
    int numOfSoftBodies = GetNumOfSoftBodies();
    fout << numOfSoftBodies << endl;
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        CParticleData& softBodyData = m_vecParticleData[i];
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        int numOfVertices = int(vecSamplePos.size());
        fout << numOfVertices << endl;
        for (int j = 0; j < numOfVertices; j++)
        {
            Vec3f& pos = vecSamplePos[j];
            fout << pos[0] << "\t" << pos[1] << "\t" << pos[2] << endl;
        }
    }
    return true;
}

bool CParticleSystem::LoadParticleSystemFromCSV(const string& fileName)
{
    std::ifstream fin(fileName.c_str());
    if (fin.fail())
    {
        std::cout << "Failed to load particle system from " << fileName << "!\n";
        return false;
    }
    std::string line_str;
    //std::getline(fin, line_str); // skip the first line
    while (std::getline(fin, line_str))
    {
        std::stringstream line_ss(line_str);
        std::string cell_str;

        std::getline(line_ss, cell_str, ',');
        std::stringstream ss1(cell_str);
        int numOfVertices;
        ss1 >> numOfVertices;
        vector<Vec3f> vertices(numOfVertices);
        for (int i = 0; i < numOfVertices; i++)
        {
            std::getline(line_ss, cell_str, ',');
            std::stringstream ss2(cell_str);
            float pos_x;
            ss2 >> pos_x;

            std::getline(line_ss, cell_str, ',');
            std::stringstream ss3(cell_str);
            float pos_y;
            ss3 >> pos_y;

            std::getline(line_ss, cell_str, ',');
            std::stringstream ss4(cell_str);
            float pos_z;
            ss4 >> pos_z;

            vertices[i] = Vec3f(pos_x, pos_y, pos_z);
        }

        CParticleData particle;
        particle.SetVecSamplePos(vertices);
        particle.SetPos(vertices[0]);
        AddParticle(particle);
    }
    std::cout << "Have loaded " << GetNumOfSoftBodies() << " particles from " << fileName << "!\n";
    return true;
}

bool CParticleSystem::SaveParticleSystemAsCSV(const string& fileName)
{
    ofstream fout(fileName.c_str());
    if (fout.fail() == true)
    {
        cout << "Failed to save particle system into " << fileName << "!\n";
        return false;
    }
    int numOfSoftBodies = GetNumOfSoftBodies();
    //fout << numOfSoftBodies << endl;
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        CParticleData& softBodyData = m_vecParticleData[i];
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        int numOfVertices = int(vecSamplePos.size());
        fout << numOfVertices << ",";
        for (int j = 0; j < numOfVertices; j++)
        {
            Vec3f& pos = vecSamplePos[j];
            fout << pos[0] << "," << pos[1] << "," << pos[2] << ",";
        }
        fout << std::endl;
    }
    return true;
}

bool CParticleSystem::LoadParticleSystem(const string& fileName)
{
    string suffix = fileName.substr(fileName.rfind(".") + 1, fileName.size());
    if (suffix == string("csv"))
    {
        return LoadParticleSystemFromCSV(fileName);
    }
    else
    {
        return LoadParticleSystemFromTXT(fileName);
    }
}

bool CParticleSystem::SaveParticleSystem(const string& fileName)
{
    string suffix = fileName.substr(fileName.rfind(".") + 1, fileName.size());
    if (suffix == string("csv"))
    {
        return SaveParticleSystemAsCSV(fileName);
    }
    else
    {
        return SaveParticleSystemAsTXT(fileName);
    }
}

void CParticleSystem::SetNeighboringSamples(Flt neighDist)
{
    const Flt neighDistSqr = neighDist * neighDist;
    const int numOfSoftBodies = GetNumOfSoftBodies();
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        m_vecParticleData[i].ClearNeighboringSamples();
    }
    int neighCount = 0;
    for (int i1 = 0; i1 < numOfSoftBodies; i1++)
    {
        CParticleData& softBodyData1 = m_vecParticleData[i1];
        vector<Vec3f>& vecSamplePos1 = softBodyData1.GetVecSamplePos();
        for (int i2 = i1 + 1; i2 < numOfSoftBodies; i2++)
        {
            CParticleData& softBodyData2 = m_vecParticleData[i2];
            vector<Vec3f>& vecSamplePos2 = softBodyData2.GetVecSamplePos();
            for (int j1 = 0; j1 < int(vecSamplePos1.size()); j1++)
            {
                for (int j2 = 0; j2 < int(vecSamplePos2.size()); j2++)
                {
                    Vec3f pr = vecSamplePos1[j1] - vecSamplePos2[j2];
                    if (mag2(pr) < neighDistSqr)
                    {
                        CNeighboringSample neigh1(pr, j1, i2, j2);
                        CNeighboringSample neigh2(-pr, j2, i1, j1);
                        softBodyData1.AddNeighboringSample(neigh1);
                        softBodyData2.AddNeighboringSample(neigh2);
                        neighCount += 2;
                    }
                } // End-For-j2
            }     // End-For-j1
        }         // End-For-i2
    }             // End-For-i1
    cout << "Average number of neighbors: " << neighCount / Flt(numOfSoftBodies) << endl;
}

void CParticleSystem::ScaleParticleSystem(Flt scaling)
{
    for (int i = 0; i < GetNumOfSoftBodies(); i++)
    {
        m_vecParticleData[i].ScaleParticleData(scaling);
    }
}
