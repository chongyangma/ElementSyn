
#include "ParticleSystem.h"

CParticleSystem::CParticleSystem()
{
}

bool CParticleSystem::LoadParticleSystemFromTXT(string fileName)
{
	ifstream fin(fileName.c_str());
	if ( fin.fail() == true )
	{
		cout << "Failed to load soft body group from file " << fileName << "!\n";
		return false;
	}
	int numOfSoftBodies;
	fin >> numOfSoftBodies;
	m_vecParticleData.resize(numOfSoftBodies);
	for ( int i=0; i<numOfSoftBodies; i++ )
	{
		int numOfVertices;
		fin >> numOfVertices;
		vector<Vec3f> vertices(numOfVertices);
		for ( int j=0; j<numOfVertices; j++ )
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

bool CParticleSystem::DumpParticleSystemToTXT(string fileName)
{
	ofstream fout(fileName.c_str());
	if ( fout.fail() == true )
	{
		cout << "Failed to dump soft body group to file " << fileName << "!\n";
		return false;
	}
	int numOfSoftBodies = GetNumOfSoftBodies();
	fout << numOfSoftBodies << endl;
	for ( int i=0; i<numOfSoftBodies; i++ )
	{
		CParticleData& softBodyData = m_vecParticleData[i];
		vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
		int numOfVertices = int(vecSamplePos.size());
		fout << numOfVertices << endl;
		for ( int j=0; j<numOfVertices; j++ )
		{
			Vec3f& pos = vecSamplePos[j];
			fout << pos[0] << "\t" << pos[1] << "\t" << pos[2] << endl;
		}
	}
	return true;
}

void CParticleSystem::SetNeighboringSamples(Flt neighDist)
{
	const Flt neighDistSqr = neighDist * neighDist;
	const int numOfSoftBodies = GetNumOfSoftBodies();
	for ( int i=0; i<numOfSoftBodies; i++ )
	{
		m_vecParticleData[i].ClearNeighboringSamples();
	}
	int neighCount = 0;
	for ( int i1=0; i1<numOfSoftBodies; i1++ )
	{
		CParticleData& softBodyData1 = m_vecParticleData[i1];
		vector<Vec3f>& vecSamplePos1 = softBodyData1.GetVecSamplePos();
		for ( int i2=i1+1; i2<numOfSoftBodies; i2++ )
		{
			CParticleData& softBodyData2 = m_vecParticleData[i2];
			vector<Vec3f>& vecSamplePos2 = softBodyData2.GetVecSamplePos();
			for ( int j1=0; j1<int(vecSamplePos1.size()); j1++ )
			{
				for ( int j2=0; j2<int(vecSamplePos2.size()); j2++ )
				{
					Vec3f pr = vecSamplePos1[j1] - vecSamplePos2[j2];
					if ( mag2(pr) < neighDistSqr )
					{
						CNeighboringSample neigh1( pr, j1, i2, j2);
						CNeighboringSample neigh2(-pr, j2, i1, j1);
						softBodyData1.AddNeighboringSample(neigh1);
						softBodyData2.AddNeighboringSample(neigh2);
						neighCount += 2;
					}
				} // End-For-j2
			} // End-For-j1
		} // End-For-i2
	} // End-For-i1
	cout << "Average number of neighbors: " << neighCount / Flt(numOfSoftBodies) << endl;
}

void CParticleSystem::ScaleParticleSystem(Flt scaling)
{
	for ( int i=0; i<GetNumOfSoftBodies(); i++ )
	{
		m_vecParticleData[i].ScaleParticleData(scaling);
	}
}
