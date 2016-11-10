
#include "MassSpringSyn.h"
#include "../DynamicElement/HungarianAlgorithm.h"
#include "../DynamicElement/PoissonDisk.h"

#define USE_TAUCS_SOLVER
#include "../DynamicElement/TAUCSsolver.h"

#include <omp.h>
#define NUM_THREADS 8

CMassSpringSyn::CMassSpringSyn()
{
	m_stepCount = 0;
	m_serialCount = 1;
	m_ptrSynConfig = new CMassSpringSynConfig;
	m_ptrCoeffMatrix = NULL;
	LoadInputData();
	InitializeOutput();
}

CMassSpringSyn::~CMassSpringSyn()
{
	DELETE_OBJECT(m_ptrSynConfig);
	DELETE_OBJECT(m_ptrCoeffMatrix);
}

void CMassSpringSyn::ResetOutput(string configFileName)
{
	m_stepCount = 0;
	DELETE_OBJECT(m_ptrCoeffMatrix);
	bool flag = m_ptrSynConfig->ReloadConfigFromFile(configFileName);
	if ( flag == false )
	{
		cout << "Have used all the configuration files!\n";
		exit(0);
	}
	LoadInputData();
	InitializeOutput();
}

void CMassSpringSyn::LoadInputData()
{
	m_vecInputSequence.clear();
	m_vecInputNeighExtended.clear();
	m_vecInputSegmentPatch.clear();
	int start = CMassSpringSynConfig::m_inputLoadStart;
	int interval = CMassSpringSynConfig::m_inputLoadInterval;
	for ( int i=start; i<=CMassSpringSynConfig::m_numOfInputFrames; i+=interval )
	{
		char fileName[MAX_PATH];
		sprintf_s(fileName, "%s%04d.txt", CMassSpringSynConfig::m_inputPrefix.c_str(), i);
		vector<CMassSpringData> strands = LoadStrandsFromTXT(fileName);
		vector<CMassSpringData> strandsNew;
#if 1 // Remove the first segment in each input branch...
		for ( int j=0; j<int(strands.size()); j++ )
		{
			CMassSpringData& strand = strands[j];
			int numOfMasses = strand.GetNumOfMasses();
			Vec3f pr = strand.GetMass(1).GetPos() - strand.GetMass(0).GetPos();
			CMassSpringData strandNew;
			strandNew.ResizeMassSpringData(numOfMasses-1);
			for ( int n=0; n<numOfMasses-1; n++ )
			{
				CMass& mass = strandNew.GetMass(n);
				mass.SetPos(strand.GetMass(n+1).GetPos()-pr);
			}
			strandsNew.push_back(strandNew);
		}
#endif
		if ( strands.empty() == true )
		{
			break;
		}
		else
		{
			//m_vecInputSequence.push_back(strands);
			m_vecInputSequence.push_back(strandsNew);
		}
	}
	// Set sequence velocity for detail exemplar...
	SetSequencesVelocity(m_vecInputSequence, CMassSpringSynConfig::m_flagTemporallyToroidal);
	// Get input segments...
	for ( int i=0; i<int(m_vecInputSequence.size()); i++ )
	{
		vector<CMassSpringData>& strands = m_vecInputSequence[i];
		for ( int j=0; j<int(strands.size()); j++ )
		{
			vector<SegmentPatch> segmentPatches = GetSegmentPatch(strands[j]);
			for ( int k=0; k<int(segmentPatches.size()); k++ )
			{
				m_vecInputSegmentPatch.push_back(segmentPatches[k]);
			}
		}
	}
	GatherInputNeighborsExtended();
	CMassSpringData& inputStrand0 = m_vecInputSequence[0][0];
	int numOfMasses = inputStrand0.GetNumOfMasses();
	CMassSpringSynConfig::m_springLength = abs(inputStrand0.GetMass(0).GetPos()[1] - inputStrand0.GetMass(numOfMasses-1).GetPos()[1]) / Flt(numOfMasses - 1);
}

void CMassSpringSyn::GatherInputNeighborsExtended()
{
	m_vecInputNeighExtended.clear();
	int neighCount = 0;
	for ( int n=0; n<int(m_vecInputSequence.size()); n++ )
	{
		MassSpringSequence& strands = m_vecInputSequence[n];
		for ( int i=0; i<int(strands.size()); i++ )
		{
			for ( int j=1; j<strands[i].GetNumOfMasses(); j++ )
			{
				MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended(n, i, j, CMassSpringSynConfig::m_neighDist * 1.5f);
				neighCount += int(neigh.m_vecNeighPr3.size());
				m_vecInputNeighExtended.push_back(neigh);
			}
		}
	}
	cout << "Average neighboring samples from other strands: " << neighCount / Flt(m_vecInputNeighExtended.size()) << endl;
}

void CMassSpringSyn::InitializeOutput()
{
	m_vecOutputSequence.clear();
	m_vecOutputStrand = InitializeStrandsViaPoissonDisk(CMassSpringSynConfig::m_exemplarName);
	for ( int n=0; n<int(m_vecOutputStrand.size()); n++ )
	{
		int numOfMasses = m_vecOutputStrand[n].GetNumOfMasses();
		vector<Vec3f> vecPr = GetConnectedInputPatch(numOfMasses - 1);
		for ( int i=1; i<numOfMasses; i++ )
		{
			CMass& mi = m_vecOutputStrand[n].GetMass(i);
			mi.SetPos(m_vecOutputStrand[n].GetMass(0).GetPos() + vecPr[i-1]);
		}
	}
	for ( int i=0; i<10; i++ )
	{
		UpdateStrandsExtended();
	}
}

void CMassSpringSyn::ResetCoeffMatrix()
{
	if ( CMassSpringSynConfig::m_neighDist <= 0.0f && CMassSpringSynConfig::m_strandRepulsionConstant <= 0.0f && m_ptrCoeffMatrix != NULL ) return;
	const int numOfStrands = int(m_vecOutputStrand.size());
	const int numOfUnknownsTotal = GetNumOfUnknownsTotal(m_vecOutputStrand);
	const int neighSize = CMassSpringSynConfig::m_neighSize;
	const Flt shapeWt = CMassSpringSynConfig::m_shapeWt;
	const Flt gaussianSigma = CMassSpringSynConfig::m_gaussianSigma;
	DELETE_OBJECT(m_ptrCoeffMatrix);
	m_ptrCoeffMatrix = new CCrossList(numOfUnknownsTotal, numOfUnknownsTotal);
	Flt temporalWtSum = 0.0f;
	for ( int i=1; i<=CMassSpringSynConfig::m_synthesisWindowSize; i++ )
	{
		temporalWtSum += CMassSpringSynConfig::m_temporalWt * exp(CMassSpringSynConfig::m_temporalSigma * i * i);
	}
	int n;
#pragma omp parallel default(shared) private(n)
	{
#pragma omp for
		for ( n=0; n<numOfStrands; n++ )
		{
			int startIdx = m_vecOutputStrand[n].GetStartIdx();
			int numOfUnknownsPerStrand = m_vecOutputStrand[n].GetNumOfMasses() - 1;
			for ( int i=1; i<=numOfUnknownsPerStrand; i++ )
			{
				m_ptrCoeffMatrix->UpdateDiagonalVal(startIdx+i-1, 1.0f);
				for ( int j=1; j<=neighSize; j++ )
				{
					Flt wt = shapeWt * exp(gaussianSigma * (j-1) * (j-1));
					if ( i-j > 0 )
					{
						m_ptrCoeffMatrix->UpdatePairVals(startIdx+i-1, startIdx+i-1-j, wt);
					}
					else if ( i-j == 0 )
					{
						m_ptrCoeffMatrix->UpdateDiagonalVal(startIdx+i-1, wt);
					}
					else
					{
						break;
					}
				}
				for ( int j=1; j<=neighSize; j++ )
				{
					if ( i+j <= numOfUnknownsPerStrand )
					{
						Flt wt = shapeWt * exp(gaussianSigma * (j-1) * (j-1));
						m_ptrCoeffMatrix->UpdatePairVals(startIdx+i-1, startIdx+i-1+j, wt);
					}
					else
					{
						break;
					}
				}
			}
		}
	}
}

void CMassSpringSyn::UpdateOutput()
{
	m_stepCount ++;
	//cout << "Step: " << m_stepCount << endl;
	UpdateStrandsExtended();
}

void CMassSpringSyn::UpdateStrandsExtended()
{
	ResetCoeffMatrix();
	const int numOfStrands = int(m_vecOutputStrand.size());
	const int numOfUnknownsTotal = GetNumOfUnknownsTotal(m_vecOutputStrand);
	const int neighSize = CMassSpringSynConfig::m_neighSize;
	const Flt shapeWt = CMassSpringSynConfig::m_shapeWt;
	const Flt neighWt = CMassSpringSynConfig::m_wtNeighPrDiff;
	const Flt smoothWt = CMassSpringSynConfig::m_smoothWt;
	const Flt temporalWt = CMassSpringSynConfig::m_temporalWt;
	const Flt gaussianSigma = CMassSpringSynConfig::m_gaussianSigma;
	const Flt temporalSigma = CMassSpringSynConfig::m_temporalSigma;
	const Flt dt = CMassSpringSynConfig::m_timeStep;
	vector<Flt> vecPx(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPy(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPz(numOfUnknownsTotal, 0.0f);
	m_vecCx.clear();
	m_vecCy.clear();
	m_vecCz.clear();
	m_vecCx.resize(numOfUnknownsTotal, 0.0f);
	m_vecCy.resize(numOfUnknownsTotal, 0.0f);
	m_vecCz.resize(numOfUnknownsTotal, 0.0f);
	int n;
	Flt distSum = 0.0f;
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel default(shared) private(n)
	{
#pragma omp for
		for ( n=0; n<numOfUnknownsTotal; n++ )
		{
			int strandIdx, sampleIdx;
			for ( int i=0; i<numOfStrands; i++ )
			{
				if ( m_vecOutputStrand[i].GetStartIdx() > n )
				{
					strandIdx = i - 1;
					sampleIdx = n - m_vecOutputStrand[strandIdx].GetStartIdx() + 1;
					break;
				}
				if ( i == numOfStrands - 1 )
				{
					strandIdx = i;
					sampleIdx = n - m_vecOutputStrand[strandIdx].GetStartIdx() + 1;
					break;
				}
			}
			CMassSpringData& outputStrand = m_vecOutputStrand[strandIdx];
			Vec3f p0 = outputStrand.GetMass(0).GetPos();
			MassNeighborhoodExtended neighOut = GetMassNeighborhoodExtended(strandIdx, sampleIdx, 0.0f);
			int nearestNeighIdx = -1;
			vector<int> matchIndices;
			Flt distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
			if ( nearestNeighIdx < 0 )
			{
				neighOut = GetMassNeighborhoodExtended(strandIdx, sampleIdx, 0.0f);
				distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
			}
			distSum += distMin;
			MassNeighborhoodExtended& neighIn = m_vecInputNeighExtended[nearestNeighIdx];
			CMass& mass = outputStrand.GetMass(sampleIdx);
			mass.SetSrcIdx(nearestNeighIdx);
			Vec3f pos = mass.GetPos();
			Vec4f quat = mass.GetQuat();
			mass.SetPos(pos);
			int idxi = n;
			vecPx[idxi] += pos[0];
			vecPy[idxi] += pos[1];
			vecPz[idxi] += pos[2];
			m_vecCx[idxi] += pos[0];
			m_vecCy[idxi] += pos[1];
			m_vecCz[idxi] += pos[2];
			vector<Vec3f>& vecNeighPr1 = neighIn.m_vecNeighPr1;
			vector<Vec3f>& vecNeighPr2 = neighIn.m_vecNeighPr2;
			vector<Vec3f>& vecNeighPrOut1 = neighOut.m_vecNeighPr1;
			for ( int j=1; j<=int(vecNeighPrOut1.size()); j++ )
			{
				int sizeDiff = int(vecNeighPr1.size() - vecNeighPrOut1.size());
				int idxj = n - j;
				Vec3f pr1 = vecNeighPr1[j-1+sizeDiff];
				Flt wt = shapeWt * exp(gaussianSigma * (j-1) * (j-1));
				if ( sampleIdx-j == 0 )
				{
					Vec3f pr1New = pr1 + p0;
					UpdateCoeffMatDiagonalVal(idxi, wt, pr1New, false);
				}
				else
				{
					UpdateCoeffMatPairVals(idxi, idxj, wt, pr1, false);
				}
			}
			for ( int j=1; j<=int(vecNeighPr2.size()); j++ )
			{
				int idxj = n + j;
				Vec3f pr2 = vecNeighPr2[j-1];
				Flt wt = shapeWt * exp(gaussianSigma * (j-1) * (j-1));
				if ( sampleIdx-1+j < m_vecOutputStrand[strandIdx].GetNumOfMasses()-1 )
				{
					UpdateCoeffMatPairVals(idxi, idxj, wt, pr2, false);
				}
			}
			// temporal part...
			vector<Vec3f>& vecPr0 = neighIn.m_vecPr0;
			for ( int i=1; i<=CMassSpringSynConfig::m_synthesisWindowSize; i++ )
			{
				Flt wt = temporalWt * exp(temporalSigma * i * i);
				if ( i <= int(neighOut.m_vecPr0.size()) )
				{
					Vec3f pr = vecPr0[i-1];
					Vec3f posPrev = m_vecOutputSequence[m_vecOutputSequence.size()-i][strandIdx].GetMass(sampleIdx).GetPos();
					Vec3f posPredict = posPrev + pr;
					UpdateCoeffMatDiagonalVal(idxi, wt, posPredict, false);
				}
				else
				{
					UpdateCoeffMatDiagonalVal(idxi, wt, pos, false);
				}
			}
			// spatial part from other strands...
			neighOut = GetMassNeighborhoodExtended(strandIdx, sampleIdx, CMassSpringSynConfig::m_neighDist);
			nearestNeighIdx = -1;
			distMin = GetNearestInputNeighborExtendedNew(neighOut, nearestNeighIdx, matchIndices);
			if ( nearestNeighIdx >= 0 )
			{
				distSum += distMin;
				MassNeighborhoodExtended& neighIn2 = m_vecInputNeighExtended[nearestNeighIdx];
				vector<Vec3f>& vecNeighPr3 = neighIn2.m_vecNeighPr3;
				for ( int j=0; j<int(matchIndices.size()); j++ )
				{
					Vec3f pr3 = vecNeighPr3[matchIndices[j]];
					int idxj = neighOut.m_vecNeighIndex[j];
					UpdateCoeffMatPairVals(idxi, idxj, neighWt, pr3);
				}
			}
		}
	} // End-OMP-Parallel
	cout << "distSum = " << distSum << endl;
	// Optimization with shape terms...
	CollisionResponse(m_vecCx, m_vecCy, m_vecCz);
	//CTAUCSsolver taucsSolver;
	vector<Flt> vecPxNew = machy_math::GetSolution(m_ptrCoeffMatrix, m_vecCx);
	vector<Flt> vecPyNew = machy_math::GetSolution(m_ptrCoeffMatrix, m_vecCy);
	vector<Flt> vecPzNew = machy_math::GetSolution(m_ptrCoeffMatrix, m_vecCz);
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[n];
		int startIdx = outputStrand.GetStartIdx();
		for ( int i=1; i<=outputStrand.GetNumOfMasses()-1; i++ )
		{
			int idxi = startIdx + i - 1;
			Vec3f pos = Vec3f(vecPxNew[idxi], vecPyNew[idxi], vecPzNew[idxi]);
			outputStrand.GetMass(i).SetPos(pos);
		}
	}
	m_vecOutputSequence.push_back(m_vecOutputStrand);
	while ( int(m_vecOutputSequence.size()) > CMassSpringSynConfig::m_synthesisWindowSize )
	{
		m_vecOutputSequence.erase(m_vecOutputSequence.begin());
	}
}

void CMassSpringSyn::CollisionResponse(vector<Flt>& vecCx, vector<Flt>& vecCy, vector<Flt>& vecCz)
{
	const Flt repulsionConstant = CMassSpringSynConfig::m_strandRepulsionConstant;
	if ( repulsionConstant <= 0.0f ) return;
	const Flt r = CMassSpringSynConfig::m_strandRadius;
	const Flt r2 = r * r;
	const int numOfUnknownsPerStrand = CMassSpringSynConfig::m_numOfMasses - 1;
	const int numOfStrands = int(m_vecOutputStrand.size());
	int repulsionCount = 0;
	for ( int i1=0; i1<numOfStrands; i1++ )
	{
		CMassSpringData& outputStrand1 = m_vecOutputStrand[i1];
		int numOfMass1 = outputStrand1.GetNumOfMasses();
		int startIdx1 = m_vecOutputStrand[i1].GetStartIdx();
		Vec3f p01 = outputStrand1.GetMass(0).GetPos();
		for ( int j1=0; j1<numOfMass1; j1++ )
		{
			Vec3f p1 = outputStrand1.GetMass(j1).GetPos();
			for ( int i2=i1+1; i2<numOfStrands; i2++ ) // No self-collision
			{
				CMassSpringData& outputStrand2 = m_vecOutputStrand[i2];
				int numOfMass2 = outputStrand2.GetNumOfMasses();
				int startIdx2 = m_vecOutputStrand[i2].GetStartIdx();
				Vec3f p02 = outputStrand2.GetMass(0).GetPos();
				for ( int j2=0; j2<numOfMass2; j2++ )
				{
					if ( j1 == 0 || j2 == 0 ) continue;
					Vec3f p2 = outputStrand2.GetMass(j2).GetPos();
					Flt distSqr = dist2(p1, p2);
					if ( distSqr < r2 )
					{
						repulsionCount ++;
						Vec3f repulsionPr = p1 - p2;
						if ( distSqr == 0.0f )
						{
							for ( int i=0; i<3; i++ )
							{
								repulsionPr[i] = rand() / Flt(RAND_MAX);
							}
						}
						normalize(repulsionPr);
						repulsionPr = repulsionPr * r;
						int idx1 = startIdx1 + j1 - 1;
						int idx2 = startIdx2 + j2 - 1;
						if ( j1 == 0 )
						{
							m_ptrCoeffMatrix->UpdateDiagonalVal(idx2, repulsionConstant);
							Vec3f prNew = p01 - repulsionPr;
							vecCx[idx2] += repulsionConstant * prNew[0];
							vecCy[idx2] += repulsionConstant * prNew[1];
							vecCz[idx2] += repulsionConstant * prNew[2];
						}
						else if ( j2 == 0 )
						{
							m_ptrCoeffMatrix->UpdateDiagonalVal(idx1, repulsionConstant);
							Vec3f prNew = p02 + repulsionConstant;
							vecCx[idx1] += repulsionConstant * prNew[0];
							vecCy[idx1] += repulsionConstant * prNew[1];
							vecCz[idx1] += repulsionConstant * prNew[2];
						}
						else
						{
							m_ptrCoeffMatrix->UpdatePairVals(idx1, idx2, repulsionConstant);
							vecCx[idx1] += repulsionConstant * repulsionPr[0];
							vecCy[idx1] += repulsionConstant * repulsionPr[1];
							vecCz[idx1] += repulsionConstant * repulsionPr[2];
							vecCx[idx2] -= repulsionConstant * repulsionPr[0];
							vecCy[idx2] -= repulsionConstant * repulsionPr[1];
							vecCz[idx2] -= repulsionConstant * repulsionPr[2];
						}
					}
				} // End-For-j2
			} // End-For-i2
		} // End-For-j1
	} // End-For-i1
	//cout << "Have detected " << repulsionCount << " pairs of collision.\n";
}

void CMassSpringSyn::RenderOutput()
{
	// Render output...
	glPushMatrix();
	glTranslatef(0.4f, -2.5f, 0.f);
	for ( int i=0; i<int(m_vecOutputStrand.size()); i++ )
	{
		m_vecOutputStrand[i].RenderMassSpringData();
	}
	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	glPopMatrix();
	// Render input exemplar...
	int winWd = glutGet(GLUT_WINDOW_WIDTH);
	int winHt = glutGet(GLUT_WINDOW_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, float(winWd)/float(winHt), 0.01, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	Vec3f posFrom = Vec3f(0.0f, 0.0f, 3.5f);
	Vec3f posAt = posFrom - Vec3f(0.0f, 0.0f, 1.0f);
	Vec3f up = Vec3f(0.0f, 1.0f, 0.0f);
	gluLookAt(posFrom[0], posFrom[1], posFrom[2], posAt[0], posAt[1], posAt[2], up[0], up[1], up[2]);
	glPushMatrix();
	glTranslatef(-1.6f, 0.f, 0.f);
	int idxCnt = min(m_stepCount, int(m_vecInputSequence.size())-1);
	idxCnt = m_stepCount % int(m_vecInputSequence.size());
	MassSpringSequence& strands = m_vecInputSequence[idxCnt];
	for ( int i=0; i<int(strands.size()); i++ )
	{
		strands[i].RenderMassSpringData(Vec3f(0.f, 0.f, 0.f));
	}
	glPopMatrix();
}

void CMassSpringSyn::RestartSynthesis()
{
	char fileName[MAX_PATH];
	sprintf_s(fileName, "MassSpringSyn_config%02d.txt", m_serialCount++);
	ResetOutput(fileName);
	m_stepCount = 0;
}

MassNeighborhoodExtended CMassSpringSyn::GetMassNeighborhoodExtended(vector<CMassSpringData>& strands, int strandIdx, int massIdx, Flt neighDist)
{
	CMassSpringData& msData = strands[strandIdx];
	const int numOfMasses = msData.GetNumOfMasses();
	const int neighSize = CMassSpringSynConfig::m_neighSize;
	CMass& mass = msData.GetMass(massIdx);
	MassNeighborhoodExtended neigh;
	neigh.m_massIdx = massIdx;
	Vec3f p0 = mass.GetPos();
	neigh.m_massPos = p0;
	for ( int i=1; i<=neighSize; i++ )
	{
		int idx1 = massIdx - i;
		if ( idx1 < 0 )
		{
			break;
		}
		Vec3f p1 = msData.GetMass(idx1).GetPos();
		Vec3f pr1 = p0 - p1;
		neigh.m_vecNeighPr1.push_back(pr1);
	}
	for ( int i=1; i<=neighSize; i++ )
	{
		int idx2 = massIdx + i;
		if ( idx2 > numOfMasses-1 )
		{
			break;
		}
		Vec3f p2 = msData.GetMass(idx2).GetPos();
		Vec3f pr2 = p0 - p2;
		neigh.m_vecNeighPr2.push_back(pr2);
	}
	if ( neighDist > 0.0f )
	{
		const Flt neighDistSqr = neighDist * neighDist;
		for ( int i=0; i<int(strands.size()); i++ )
		{
			if ( i == strandIdx )
			{
				continue;
			}
			CMassSpringData& strand = strands[i];
			for ( int j=1; j<int(strand.GetNumOfMasses()); j++ )
			{
				CMass& massNeigh = strand.GetMass(j);
				Vec3f p3 = massNeigh.GetPos();
				if ( dist2(p0, p3) < neighDistSqr )
				{
					Vec3f pr3 = p0 - p3;
					neigh.m_vecNeighPr3.push_back(pr3);
					int neighIndex = strand.GetStartIdx() + j - 1;
					neigh.m_vecNeighIndex.push_back(neighIndex);
				}
			}
		}
	}
	return neigh;
}

MassNeighborhoodExtended CMassSpringSyn::GetMassNeighborhoodExtended(int frameIdx, int strandIdx, int massIdx, Flt neighDist)
{
	vector<CMassSpringData>& strands = m_vecInputSequence[frameIdx];
	MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended(strands, strandIdx, massIdx, neighDist);
	const int sequenceLength = int(m_vecInputSequence.size());
	for ( int i=1; i<=CMassSpringSynConfig::m_synthesisWindowSize; i++ )
	{
		int frameIdxPrev = frameIdx - i;
		frameIdxPrev = (frameIdxPrev + sequenceLength) % sequenceLength;
		vector<CMassSpringData>& strandsPrev = m_vecInputSequence[frameIdxPrev];
		MassNeighborhoodExtended neighPrev = GetMassNeighborhoodExtended(strandsPrev, strandIdx, massIdx, 0.0f);
		neigh.m_vecNeighPr1Prev.push_back(neighPrev.m_vecNeighPr1);
		neigh.m_vecNeighPr2Prev.push_back(neighPrev.m_vecNeighPr2);
		neigh.m_vecNeighPr3Prev.push_back(neighPrev.m_vecNeighPr3);
		neigh.m_vecPr0.push_back(neigh.m_massPos - neighPrev.m_massPos);
	}
	return neigh;
}

MassNeighborhoodExtended CMassSpringSyn::GetMassNeighborhoodExtended(int strandIdx, int massIdx, Flt neighDist)
{
	MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended(m_vecOutputStrand, strandIdx, massIdx, neighDist);
	const int iMax = min(CMassSpringSynConfig::m_synthesisWindowSize, int(m_vecOutputSequence.size()));
	for ( int i=1; i<=iMax; i++ )
	{
		vector<CMassSpringData>& strandsPrev = m_vecOutputSequence[m_vecOutputSequence.size()-i];
		MassNeighborhoodExtended neighPrev = GetMassNeighborhoodExtended(strandsPrev, strandIdx, massIdx, 0.0f);
		neigh.m_vecNeighPr1Prev.push_back(neighPrev.m_vecNeighPr1);
		neigh.m_vecNeighPr2Prev.push_back(neighPrev.m_vecNeighPr2);
		neigh.m_vecNeighPr3Prev.push_back(neighPrev.m_vecNeighPr3);
		neigh.m_vecPr0.push_back(neigh.m_massPos - neighPrev.m_massPos);
	}
	return neigh;
}

Flt CMassSpringSyn::NeighborhoodMetricExtended(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices)
{
	Flt distSum = 0.0f;
	int idxDiff = neigh1.m_massIdx - neigh2.m_massIdx;
	distSum += idxDiff * idxDiff * CMassSpringSynConfig::m_wtIdxDiff;
	int size11 = int(neigh1.m_vecNeighPr1.size());
	int size21 = int(neigh2.m_vecNeighPr1.size());
	int size12 = int(neigh1.m_vecNeighPr2.size());
	int size22 = int(neigh2.m_vecNeighPr2.size());
	int size13 = int(neigh1.m_vecNeighPr3.size());
	int size23 = int(neigh2.m_vecNeighPr3.size());
	if ( (size11 != size21) || (size12 != size22) )
	{
		return 1e10;
	}
	if ( CMassSpringSynConfig::m_neighDist > 0.0f && CMassSpringSynConfig::m_wtNeighPrDiff > 0.0f && (size13 > size23) )
	{
		return 1e10;
	}
	if ( neigh1.m_vecNeighPr1Prev.size() > neigh2.m_vecNeighPr1Prev.size() )
	{
		return 1e10;
	}
	int neighPrNum1 = size11;
	int neighPrNum2 = size12;
	distSum += VecPrDifference(neigh1.m_vecNeighPr1, neigh2.m_vecNeighPr1) * CMassSpringSynConfig::m_wtPrDiff;
	distSum += VecPrDifference(neigh1.m_vecNeighPr2, neigh2.m_vecNeighPr2) * CMassSpringSynConfig::m_wtPrDiff;
	for ( int i=0; i<int(neigh1.m_vecNeighPr1Prev.size()); i++ )
	{
		Flt i2 = (i + 1) * (i + 1);
		Flt wt = CMassSpringSynConfig::m_wtPrDiff * exp(CMassSpringSynConfig::m_temporalSigma * i2);
		distSum += VecPrDifference(neigh1.m_vecNeighPr1Prev[i], neigh2.m_vecNeighPr1Prev[i]) * wt;
		distSum += VecPrDifference(neigh1.m_vecNeighPr2Prev[i], neigh2.m_vecNeighPr2Prev[i]) * wt;
	}
	distSum += VecPrDifferenceTemporal(neigh1.m_vecPr0, neigh2.m_vecPr0) * CMassSpringSynConfig::m_temporalWt;
	if ( CMassSpringSynConfig::m_neighDist <= 0.0f || CMassSpringSynConfig::m_wtNeighPrDiff <= 0.0f || size13 == 0 )
	{
		return distSum;
	}
	matchIndices.resize(size13);
	vector<double> distMatrix(size13 * size23);
	int cnt = 0;
	for ( int j=0; j<size23; j++ )
	{
		Vec3f pr2 = neigh2.m_vecNeighPr3[j];
		for ( int i=0; i<size13; i++ )
		{
			Vec3f pr1 = neigh1.m_vecNeighPr3[i];
			Vec3f prr = pr1 - pr2;
			Flt distTmp = mag2(prr);
			distMatrix[cnt++] = distTmp;
		}
	}
	vector<int> assignment(size13);
	double cost = 0.0f;
	CHungarianAlgorithm ha;
	ha.FindOptimalAssignment(assignment, &cost, distMatrix, size13, size23);
	distSum += cost * CMassSpringSynConfig::m_wtNeighPrDiff;
	for ( int i=0; i<size13; i++ )
	{
		matchIndices[i] = assignment[i] - 1;
	}
	return distSum;
}

Flt CMassSpringSyn::GetNearestInputNeighborExtended(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices)
{
	Flt distMin = 1e10;
	for ( int i=0; i<int(m_vecInputNeighExtended.size()); i++ )
	{
		vector<int> matchIndicesTmp;
		Flt distTmp = NeighborhoodMetricExtended(neighOut, m_vecInputNeighExtended[i], matchIndicesTmp);
		if ( distTmp < distMin )
		{
			distMin = distTmp;
			nearesNeighIdx = i;
			matchIndices = matchIndicesTmp;
		}
	}
	return distMin;
}

Flt CMassSpringSyn::NeighborhoodMetricExtendedNew(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices)
{
	Flt distSum = 0.0f;
	int size13 = int(neigh1.m_vecNeighPr3.size());
	int size23 = int(neigh2.m_vecNeighPr3.size());
	if ( CMassSpringSynConfig::m_neighDist > 0.0f && CMassSpringSynConfig::m_wtNeighPrDiff > 0.0f && (size13 > size23) )
	{
		return 1e10;
	}
	if ( CMassSpringSynConfig::m_neighDist <= 0.0f || CMassSpringSynConfig::m_wtNeighPrDiff <= 0.0f || size13 == 0 )
	{
		return distSum;
	}
	matchIndices.resize(size13);
	vector<double> distMatrix(size13 * size23);
	int cnt = 0;
	for ( int j=0; j<size23; j++ )
	{
		Vec3f pr2 = neigh2.m_vecNeighPr3[j];
		for ( int i=0; i<size13; i++ )
		{
			Vec3f pr1 = neigh1.m_vecNeighPr3[i];
			Vec3f prr = pr1 - pr2;
			Flt distTmp = mag2(prr);
			distMatrix[cnt++] = distTmp;
		}
	}
	vector<int> assignment(size13);
	double cost = 0.0f;
	CHungarianAlgorithm ha;
	ha.FindOptimalAssignment(assignment, &cost, distMatrix, size13, size23);
	distSum += cost * CMassSpringSynConfig::m_wtNeighPrDiff;
	for ( int i=0; i<size13; i++ )
	{
		matchIndices[i] = assignment[i] - 1;
	}
	return distSum;
}

Flt CMassSpringSyn::GetNearestInputNeighborExtendedNew(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices)
{
	Flt distMin = 1e10;
	for ( int i=0; i<int(m_vecInputNeighExtended.size()); i++ )
	{
		vector<int> matchIndicesTmp;
		Flt distTmp = NeighborhoodMetricExtendedNew(neighOut, m_vecInputNeighExtended[i], matchIndicesTmp);
		if ( distTmp < distMin )
		{
			distMin = distTmp;
			nearesNeighIdx = i;
			matchIndices = matchIndicesTmp;
		}
	}
	return distMin;
}

vector<SegmentPatch> CMassSpringSyn::GetSegmentPatch(CMassSpringData& msData)
{
	vector<SegmentPatch> segmentPatches;
	const int patchSize = CMassSpringSynConfig::m_patchSize;
	int numOfMasses = msData.GetNumOfMasses();
	for ( int i=0; i<numOfMasses-patchSize+1; i++ )
	{
		SegmentPatch patch;
		patch.m_flagStart = (i == 0) ? true : false;
		for ( int j=0; j<patchSize; j++ )
		{
			CMass& mass = msData.GetMass(i+j);
			patch.m_vecPos.push_back(mass.GetPos());
			patch.m_vecVel.push_back(mass.GetVel());
		}
		segmentPatches.push_back(patch);
	}
	return segmentPatches;
}

Flt CMassSpringSyn::SegmentPatchCompatibility(SegmentPatch* ptrPatch1, SegmentPatch* ptrPatch2)
{
	Flt compatibility = 0.0f;
	int patchSize1 = int(ptrPatch1->m_vecPos.size());
	Vec3f pos11 = ptrPatch1->m_vecPos[patchSize1-2];
	Vec3f pos12 = ptrPatch1->m_vecPos[patchSize1-1];
	Vec3f pos21 = ptrPatch2->m_vecPos[0];
	Vec3f pos22 = ptrPatch2->m_vecPos[1];
	Vec3f pr1 = pos12 - pos11;
	Vec3f pr2 = pos22 - pos21;
	compatibility += dist2(pr1, pr2);
	Vec3f vel11 = ptrPatch1->m_vecVel[patchSize1-2];
	Vec3f vel12 = ptrPatch1->m_vecVel[patchSize1-1];
	Vec3f vel21 = ptrPatch2->m_vecVel[0];
	Vec3f vel22 = ptrPatch2->m_vecVel[1];
	const Flt velDiffWt = 1.f;
	compatibility += dist2(vel11, vel21) * velDiffWt;
	compatibility += dist2(vel12, vel22) * velDiffWt;
	return compatibility;
}

vector<Vec3f> CMassSpringSyn::GetConnectedInputPatch(int patchLength)
{
	vector<Vec3f> vecPatchPr;
	const int patchSize = CMassSpringSynConfig::m_patchSize;
	int numOfInputPatches = int(m_vecInputSegmentPatch.size());
	bool flagStart = false;
	int patchIdx;
	SegmentPatch patch0;
	while ( flagStart != true )
	{
		patchIdx = int(rand() / Flt(RAND_MAX) * numOfInputPatches);
		patchIdx = patchIdx % numOfInputPatches;
		patch0 = m_vecInputSegmentPatch[patchIdx];
		flagStart = patch0.m_flagStart;
	}
	Vec3f p0 = patch0.m_vecPos[0];
	for ( int i=1; i<patchSize; i++ )
	{
		vecPatchPr.push_back(patch0.m_vecPos[i] - p0);
	}
	SegmentPatch* ptrPatch1 = &patch0;
	while ( int(vecPatchPr.size()) < patchLength )
	{
		Flt compatibility = 1e10;
		Flt compatibilityThresh = 0.01f;
		int trialCount = 0;
		SegmentPatch patchNew;
		while ( flagStart != false || compatibility > compatibilityThresh )
		{
			trialCount ++;
			patchIdx = int(rand() / Flt(RAND_MAX) * numOfInputPatches);
			patchIdx = patchIdx % numOfInputPatches;
			patchNew = m_vecInputSegmentPatch[patchIdx];
			flagStart = patchNew.m_flagStart;
			compatibility = SegmentPatchCompatibility(ptrPatch1, &patchNew);
			if ( trialCount > 100 )
			{
				compatibilityThresh += 0.01f;
				trialCount = 0;
			}
		}
		Vec3f pend = vecPatchPr.back();
		Vec3f p1 = patchNew.m_vecPos[1];
		for ( int i=2; i<patchSize; i++ )
		{
			vecPatchPr.push_back(patchNew.m_vecPos[i] + pend - p1);
		}
		ptrPatch1 = &patchNew;
	}
	return vecPatchPr;
}

void CMassSpringSyn::UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos, bool flagUpdateMat /* = true */)
{
	if ( flagUpdateMat == true )
	{
		m_ptrCoeffMatrix->UpdateDiagonalVal(idx, wt);
	}
	m_vecCx[idx] += wt * pos[0];
	m_vecCy[idx] += wt * pos[1];
	m_vecCz[idx] += wt * pos[2];
}

void CMassSpringSyn::UpdateCoeffMatPairVals(int idxi, int idxj, Flt wt, Vec3f pr, bool flagUpdateMat /* = true */)
{
	if ( flagUpdateMat == true )
	{
		m_ptrCoeffMatrix->UpdatePairVals(idxi, idxj, wt);
	}
	m_vecCx[idxi] += wt * pr[0];
	m_vecCy[idxi] += wt * pr[1];
	m_vecCz[idxi] += wt * pr[2];
	m_vecCx[idxj] -= wt * pr[0];
	m_vecCy[idxj] -= wt * pr[1];
	m_vecCz[idxj] -= wt * pr[2];
}

vector<CMassSpringData> CMassSpringSyn::LoadStrandsFromTXT(string fileName)
{
	vector<CMassSpringData> strands;
	ifstream fin(fileName.c_str());
	if ( fin.fail() == true )
	{
		cout << "Failed to load strands from TXT file " << fileName << "!\n";
		return strands;
	}
	int numOfStrands;
	int numOfMasses;
	fin >> numOfStrands;
	//numOfStrands = 2; // reduce input branches for debug purpose
	numOfStrands = min(numOfStrands, CMassSpringSynConfig::m_numOfInputStrands);
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData strand;
		fin >> numOfMasses;
		strand.ResizeMassSpringData(numOfMasses);
		for ( int i=0; i<numOfMasses; i++ )
		{
			CMass& mass = strand.GetMass(i);
			Flt px, py, pz;
			fin >> px >> py >> pz;
			mass.SetPos(Vec3f(px, py, pz));
			Flt vx, vy, vz;
			fin >> vx >> vy >> vz;
			mass.SetVel(Vec3f(vx, vy, vz));
			Flt fx, fy, fz;
			fin >> fx >> fy >> fz;
		}
		strands.push_back(strand);
	}
	return strands;
}

vector<CMassSpringData> CMassSpringSyn::InitializeStrandsViaPoissonDisk(string fileName)
{
	PoissonDisk2f pd2f;
	Flt distSqrMinThresh = 0.04f;
	//distSqrMinThresh *= 9.0f; // reduce branches for debug purpose
	const Vec2f posMin = CMassSpringSynConfig::m_poisson2Dmin;
	const Vec2f posMax = CMassSpringSynConfig::m_poisson2Dmax;
	vector<Vec2f> vecOrigin = pd2f.GeneratePoints(posMin, posMax, distSqrMinThresh);
	const int numOfStrands = int(vecOrigin.size());
	cout << "Have initialized " << numOfStrands << " branches!\n";
	int numOfMasses = CMassSpringSynConfig::m_numOfMasses;
	Flt springLength = CMassSpringSynConfig::m_springLength;
	int startIdx = 0;
	vector<CMassSpringData> strands;
	strands.resize(numOfStrands);
	for ( int n=0; n<numOfStrands; n++ )
	{
		Vec3f pn = Vec3f(vecOrigin[n][0], 2.0f, vecOrigin[n][1]);
		Flt totalLength = pn[1] - 0.5f;
		totalLength += rand() / Flt(RAND_MAX) * 0.25f;
		numOfMasses = totalLength / springLength;
		strands[n].ResizeMassSpringData(numOfMasses);
		strands[n].SetStartIdx(startIdx);
		for ( int i=0; i<numOfMasses; i++ )
		{
			CMass& mass = strands[n].GetMass(i);
			Flt py = -i * springLength;
			mass.SetPos(pn + Vec3f(0.0f, py, 0.0f));
		}
		startIdx += numOfMasses - 1;
	}
	return strands;
}

int CMassSpringSyn::GetNumOfUnknownsTotal(vector<CMassSpringData>& strands)
{
	int numOfUnknownsTotal = 0;
	for ( int i=0; i<int(strands.size()); i++ )
	{
		numOfUnknownsTotal += strands[i].GetNumOfMasses() - 1;
	}
	return numOfUnknownsTotal;
}

void CMassSpringSyn::SetSequencesVelocity(vector<MassSpringSequence>& vecSequence, bool flagTemporallyToroidal /* = false */)
{
	const int numOfFrames = int(vecSequence.size());
	for ( int t=1; t<numOfFrames; t++ )
	{
		MassSpringSequence& sequence1 = vecSequence[t];
		MassSpringSequence& sequence0 = vecSequence[t-1];
		for ( int i=0; i<int(sequence1.size()); i++ )
		{
			CMassSpringData& msData1 = sequence1[i];
			CMassSpringData& msData0 = sequence0[i];
			for ( int j=0; j<msData1.GetNumOfMasses(); j++ )
			{
				CMass& sample1 = msData1.GetMass(j);
				CMass& sample0 = msData0.GetMass(j);
				Vec3f pos1 = sample1.GetPos();
				Vec3f pos0 = sample0.GetPos();
				Vec3f vel1 = pos1 - pos0;
				sample1.SetVel(vel1);
			}
		}
	}
	// For the first frame...
	if ( vecSequence.size() <= 1 ) return;
	if ( flagTemporallyToroidal == true )
	{
		MassSpringSequence& sequence1 = vecSequence[0];
		MassSpringSequence& sequence0 = vecSequence[vecSequence.size()-1];
		for ( int i=0; i<int(sequence1.size()); i++ )
		{
			CMassSpringData& msData1 = sequence1[i];
			CMassSpringData& msData0 = sequence0[i];
			for ( int j=0; j<msData1.GetNumOfMasses(); j++ )
			{
				CMass& sample1 = msData1.GetMass(j);
				CMass& sample0 = msData0.GetMass(j);
				Vec3f pos1 = sample1.GetPos();
				Vec3f pos0 = sample0.GetPos();
				Vec3f vel1 = pos1 - pos0;
				sample1.SetVel(vel1);
			}
		}
		return;
	}
	MassSpringSequence& sequence1 = vecSequence[1];
	MassSpringSequence& sequence0 = vecSequence[0];
	for ( int i=0; i<int(sequence1.size()); i++ )
	{
		CMassSpringData& msData1 = sequence1[i];
		CMassSpringData& msData0 = sequence0[i];
		for ( int j=0; j<msData1.GetNumOfMasses(); j++ )
		{
			CMass& sample1 = msData1.GetMass(j);
			CMass& sample0 = msData0.GetMass(j);
			sample0.SetVel(sample1.GetVel());
		}
	}
}
