
#include "SpaghettiSyn.h"
#include <omp.h>
#define USE_TAUCS_SOLVER
#ifdef USE_TAUCS_SOLVER
#include "../DynamicElement/TAUCSsolver.h"
#else
#include "../DynamicElement/MKLsolver.h"
#endif
#include "../DynamicElement/PoissonDisk.h"
#include "../DynamicElement/HungarianAlgorithm.h"
//#define RUN_SIMULATION
#define SYNTHESIZE_FRAME_BY_FRAME

CSpaghettiSyn::CSpaghettiSyn()
{
	m_stepCount = 0;
	m_serialCount = 1;
	m_ptrSynConfig = new CSpaghettiSynConfig;
	m_ptrCoeffMatrix = NULL;
	LoadInputData();
#ifndef RUN_SIMULATION
#ifndef ONLY_COARSE_MOTION
	InitializeOutput();
#endif
#endif
}

CSpaghettiSyn::~CSpaghettiSyn()
{
	DELETE_OBJECT(m_ptrSynConfig);
	DELETE_OBJECT(m_ptrCoeffMatrix);
}

void CSpaghettiSyn::LoadInputData()
{
#ifdef RUN_SIMULATION
	m_vecOutputStrand = CMassSpringData::LoadStrandsFromTXT("init.txt");
	return;
#endif
	m_vecInputSequence.clear();
	int frameCnt = 0;
	for ( int i=CSpaghettiSynConfig::m_inputLoadStart; i<=CSpaghettiSynConfig::m_inputLoadEnd; i+=CSpaghettiSynConfig::m_inputLoadInterval )
	{
		char fileName[MAX_PATH];
		sprintf_s(fileName, "%s%04d.txt", CSpaghettiSynConfig::m_inputPrefix.c_str(), i);
		vector<CMassSpringData> strands = CMassSpringData::LoadStrandsFromTXT(fileName); //, 5);
		if ( strands.empty() == true )
		{
			break;
		}
		else
		{
			m_vecInputSequence.push_back(strands);
		}
		frameCnt ++;
	}
	cout << "Have loaded " << frameCnt << " input frames!\n";
	m_vecCoarseSequence.clear();
	frameCnt = 0;
	for ( int i=CSpaghettiSynConfig::m_coarseLoadStart; i<=CSpaghettiSynConfig::m_coarseLoadEnd; i+=CSpaghettiSynConfig::m_coarseLoadInterval )
	{
		char fileName[MAX_PATH];
		sprintf_s(fileName, "%s%04d.txt", CSpaghettiSynConfig::m_coarsePrefix.c_str(), i);
		vector<CMassSpringData> strands = CMassSpringData::LoadStrandsFromTXT(fileName);
		if ( strands.empty() == true )
		{
			break;
		}
		else
		{
			m_vecCoarseSequence.push_back(strands);
		}
		frameCnt ++;
	}
	cout << "Have loaded " << frameCnt << " coarse frames!\n";
	ChangeCoarseSequenceAsNonStretchable();
	ChangeCoarseSequence0();
	ChangeCoarseSequence1();
	ChangeCoarseSequenceAsNonStretchable();
	ChangeCoarseSequence();
	CMassSpringData& strandCoarse = m_vecCoarseSequence[0][0];
	Vec3f p0 = strandCoarse.GetMass(strandCoarse.GetNumOfMasses()/2).GetPos();
	Vec3f p1 = strandCoarse.GetMass(strandCoarse.GetNumOfMasses()/2-1).GetPos();
	Flt segmentLength = dist(p0, p1);
	cout << "Segment length = " << segmentLength << endl;
	ChangeInputSequence();
	ChangeInputSequenceAccordingToCoarseMotion(segmentLength);
	//DumpChopsticks(m_vecCoarseSequence);
	if ( CSpaghettiSynConfig::m_flagTemporallyToroidal == true )
	{
		// Make input temporally toroidal...
		int originalLength = int(m_vecInputSequence.size());
		for ( int i=1; i<originalLength; i++ )
		{
			m_vecInputSequence.push_back(m_vecInputSequence[originalLength-i]);
		}
	}
	// Get input segments...
	m_vecInputPatchExtended = GetSegmentPatchesExtended(m_vecInputSequence);
	//InitializeOutputViaPoissonDiskVertically(); m_vecInputSequence = m_vecOutputSequence;
	GatherInputNeighborsExtended();
}

void CSpaghettiSyn::UpdateOutput()
{
	m_stepCount ++;
#ifdef RUN_SIMULATION
	TimeIntegration();
#else // Run synthesis...
#ifdef SYNTHESIZE_FRAME_BY_FRAME
#ifndef ONLY_COARSE_MOTION
	if ( m_stepCount >= 120 )
	{
		if ( m_stepCount == 120 )
		{
			CSpaghettiSynConfig::m_gravity = -98.0f;
			MassSpringSequence vecStrandLast = m_vecOutputSequence[m_vecOutputSequence.size()-2];
			for ( int i=0; i<int(m_vecOutputStrand.size()); i++ )
			{
				m_vecOutputStrand[i].GetMass(0).SetFlagFixed(true);
				for ( int j=0; j<m_vecOutputStrand[i].GetNumOfMasses(); j++ )
				{
					Vec3f pj = m_vecOutputStrand[i].GetMass(j).GetPos();
					Vec3f pjLast = vecStrandLast[i].GetMass(j).GetPos();
					Vec3f vj = (pj - pjLast) * (1.0f / CSpaghettiSynConfig::m_timeStep);
					m_vecOutputStrand[i].GetMass(j).SetVel(vj);
					//m_vecOutputStrand[i].GetMass(j).SetVel(Vec3f(0.f, 0.f, 0.f));
					m_vecOutputStrand[i].GetMass(j).SetForce(Vec3f(0.f, 0.f, 0.f));
				}
			}
		}
		TimeIntegration();
	}
	else
	{
		UpdateFrameViaOptimization();
	}
	//UpdateFrameViaBulletSimulation();
#endif
#else
	int frameIdx = min(m_stepCount, int(m_vecOutputSequence.size())-1);
	m_vecOutputStrand = m_vecOutputSequence[frameIdx];
	//m_vecInitialStrand = m_vecInitialSequence[frameIdx];
#endif
#endif
	char fileName[MAX_PATH];
	sprintf_s(fileName, "%sDumped\\dumped_%04d.txt", CSpaghettiSynConfig::m_outputPrefix.c_str(), GetStepCount());
	CMassSpringData::DumpStrandsToTXT(fileName, m_vecOutputStrand);
	//sprintf_s(fileName, "%sDumped\\coarse_%04d.obj", CSpaghettiSynConfig::m_outputPrefix.c_str(), GetStepCount());
	//CMassSpringData::DumpStrandsToOBJ(m_vecCoarseSequence[GetStepCount()], fileName, CSpaghettiSynConfig::m_springRadius * 0.5f, 16, 1, 5);
	//sprintf_s(fileName, "%sDumped\\exemplar_%04d.obj", CSpaghettiSynConfig::m_outputPrefix.c_str(), GetStepCount());
	//CMassSpringData::DumpStrandsToOBJ(m_vecInputSequence[GetStepCount()], fileName, CSpaghettiSynConfig::m_springRadius * 0.5f, 16, 1, 5);
	int coarseIdx = GetStepCount();
#ifndef SYNTHESIZE_FRAME_BY_FRAME
	//coarseIdx = coarseIdx % CSpaghettiSynConfig::m_iterationNum;
	coarseIdx = (coarseIdx - 1) % CSpaghettiSynConfig::m_iterationNum + 1;
#endif
#ifndef ONLY_COARSE_MOTION
	vector<CMassSpringData> finalStrands = EditCoarseStrandsFromSynthesisOutput(m_vecCoarseSequence[coarseIdx]);
	sprintf_s(fileName, "%soutput_%04d.obj", CSpaghettiSynConfig::m_outputPrefix.c_str(), GetStepCount());
	CMassSpringData::DumpStrandsToOBJ(finalStrands, fileName, CSpaghettiSynConfig::m_springRadius * 0.5f, 16, 1, 5);
	sprintf_s(fileName, "%sDumped\\syn_%04d.obj", CSpaghettiSynConfig::m_outputPrefix.c_str(), GetStepCount());
	//CMassSpringData::DumpStrandsToOBJ(m_vecOutputStrand, fileName, CSpaghettiSynConfig::m_springRadius * 0.5f, 16, 1, 5);
	sprintf_s(fileName, "%scoarse_%04d.obj", CSpaghettiSynConfig::m_outputPrefix.c_str(), coarseIdx);
	CMassSpringData::DumpStrandsToOBJ(m_vecCoarseSequence[coarseIdx], fileName, CSpaghettiSynConfig::m_springRadius * 0.5f, 16, 1, 5);
#endif
}

void CSpaghettiSyn::RenderOutput()
{
	// Render output...
	glPushMatrix();
	glTranslatef(-1.0f, 0.0f, 0.0f);
	vector<CMassSpringData> outputStrands = CMassSpringData::NurbsInterpolateStrands(m_vecOutputStrand, 5);
	for ( int i=0; i<int(m_vecOutputStrand.size()); i++ )
	{
		outputStrands[i].RenderMassSpringData();
	}
	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	glPopMatrix();
#ifndef RUN_SIMULATION
	// Render coarse motion...
	glPushMatrix();
	//glTranslatef(1.0f, 0.5f, 0.0f);
	int coarseIdx = m_stepCount % int(m_vecCoarseSequence.size());
#ifndef SYNTHESIZE_FRAME_BY_FRAME
	//coarseIdx = coarseIdx % CSpaghettiSynConfig::m_iterationNum;
	coarseIdx = (coarseIdx - 1) % CSpaghettiSynConfig::m_iterationNum + 1;
#endif
#ifdef ONLY_COARSE_MOTION
	vector<CMassSpringData>& coarseStrandsOld = m_vecCoarseSequence[coarseIdx];
#else
	vector<CMassSpringData> coarseStrandsOld = EditCoarseStrandsFromSynthesisOutput(m_vecCoarseSequence[coarseIdx]);
#endif
	vector<CMassSpringData> coarseStrands = CMassSpringData::NurbsInterpolateStrands(coarseStrandsOld, 5);
	for ( int i=0; i<int(coarseStrands.size()); i++ )
	{
		coarseStrands[i].RenderMassSpringData();
	}
	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	glPopMatrix();
	// Render input
	glPushMatrix();
	glTranslatef(-0.5f, -0.f, 0.0f);
	//int inputIdx = min(m_stepCount, int(m_vecInputSequence.size())-1);
	int inputIdx = m_stepCount % int(m_vecInputSequence.size());
	vector<CMassSpringData>& inputStrandsOld = m_vecInputSequence[inputIdx];
	vector<CMassSpringData> inputStrands = CMassSpringData::NurbsInterpolateStrands(inputStrandsOld, 5);
	for ( int i=0; i<int(inputStrands.size()); i++ )
	{
		inputStrands[i].RenderMassSpringData();
	}
	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	glPopMatrix();
#ifndef SYNTHESIZE_FRAME_BY_FRAME
	// Render initialization result...
	glPushMatrix();
	glTranslatef(0.5f, 0.0f, 0.0f);
	for ( int i=0; i<int(m_vecInitialStrand.size()); i++ )
	{
		//m_vecInitialStrand[i].RenderMassSpringData();
	}
	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	glPopMatrix();
#endif
#endif
}

void CSpaghettiSyn::RestartSynthesis()
{
	char fileName[MAX_PATH];
	sprintf_s(fileName, "SpaghettiSyn_Config%02d.txt", m_serialCount++);
	ResetOutput(fileName);
	m_stepCount = 0;
}

void CSpaghettiSyn::ResetOutput(string configFileName)
{
	m_stepCount = 0;
	bool flag = m_ptrSynConfig->ReloadConfigFromFile(configFileName);
	if ( flag == false )
	{
		cout << "Have used all the configuration files!\n";
		exit(0);
	}
	LoadInputData();
#ifndef RUN_SIMULATION
	InitializeOutput();
#endif
}

void CSpaghettiSyn::GatherInputNeighborsExtended()
{
	m_vecInputNeighExtended.clear();
	int neighCount = 0;
	for ( int n=0; n<int(m_vecInputSequence.size()); n++ )
	{
		MassSpringSequence& strands = m_vecInputSequence[n];
		for ( int i=0; i<int(strands.size()); i++ )
		{
			for ( int j=0; j<strands[i].GetNumOfMasses(); j++ )
			{
				MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended1(n, i, j, CSpaghettiSynConfig::m_spatialNeighDist * 1.5f, CSpaghettiSynConfig::m_flagTemporallyToroidal);
				neighCount += int(neigh.m_vecNeighPr3.size());
				m_vecInputNeighExtended.push_back(neigh);
			}
		}
	}
	cout << "Average neighboring samples from other strands: " << neighCount / Flt(m_vecInputNeighExtended.size()) << endl;
}

Flt CSpaghettiSyn::GetNearestInputNeighborExtended(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices)
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

Flt CSpaghettiSyn::NeighborhoodMetricExtended(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices)
{
	Flt distSum = 0.0f;
#ifdef HACK_WITH_RELATIVE_POSITION
	Flt relativePos1 = neigh1.m_relativePos;
	Flt relativePos2 = neigh2.m_relativePos;
	relativePos2 = abs(relativePos2 - 0.5f) * 2.0f;
	Flt relativePosDiff = relativePos1 - relativePos2;
	distSum += relativePosDiff * relativePosDiff * 0.01f;
#endif
#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
	distSum += neigh2.m_neighUsedCount * 0.01f;
#endif
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
	if ( CSpaghettiSynConfig::m_spatialNeighDist > 0.0f && CSpaghettiSynConfig::m_spatialWt > 0.0f && (size13 > size23) )
	{
		return 1e10;
	}
	if ( neigh1.m_vecNeighPr1Prev.size() > neigh2.m_vecNeighPr1Prev.size() )
	{
		return 1e10;
	}
	int neighPrNum1 = size11;
	int neighPrNum2 = size12;
	distSum += VecPrDifference(neigh1.m_vecNeighPr1, neigh2.m_vecNeighPr1) * CSpaghettiSynConfig::m_spatialWt;
	distSum += VecPrDifference(neigh1.m_vecNeighPr2, neigh2.m_vecNeighPr2) * CSpaghettiSynConfig::m_spatialWt;
	for ( int i=0; i<int(neigh1.m_vecNeighPr1Prev.size()); i++ )
	{
		Flt i2 = (i + 1) * (i + 1);
		Flt wt = CSpaghettiSynConfig::m_temporalWt * exp(CSpaghettiSynConfig::m_temporalSigma * i2);
		distSum += VecPrDifference(neigh1.m_vecNeighPr1Prev[i], neigh2.m_vecNeighPr1Prev[i]) * wt;
		distSum += VecPrDifference(neigh1.m_vecNeighPr2Prev[i], neigh2.m_vecNeighPr2Prev[i]) * wt;
	}
	distSum += VecPrDifferenceTemporal(neigh1.m_vecPr0, neigh2.m_vecPr0) * CSpaghettiSynConfig::m_temporalWt;
	if ( CSpaghettiSynConfig::m_spatialNeighDist <= 0.0f || CSpaghettiSynConfig::m_spatialWt <= 0.0f || size13 == 0 )
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
	distSum += cost * CSpaghettiSynConfig::m_spatialWt;
	for ( int i=0; i<size13; i++ )
	{
		matchIndices[i] = assignment[i] - 1;
	}
	return distSum;
}

Flt CSpaghettiSyn::GetNearestInputNeighborExtendedNew(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices)
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

Flt CSpaghettiSyn::NeighborhoodMetricExtendedNew(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices)
{
	Flt distSum = 0.0f;
	int size13 = int(neigh1.m_vecNeighPr3.size());
	int size23 = int(neigh2.m_vecNeighPr3.size());
	if ( CSpaghettiSynConfig::m_spatialNeighDist > 0.0f && CSpaghettiSynConfig::m_spatialWt > 0.0f && (size13 > size23) )
	{
		return 1e10;
	}
	if ( CSpaghettiSynConfig::m_spatialNeighDist <= 0.0f || CSpaghettiSynConfig::m_spatialWt <= 0.0f || size13 == 0 )
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
	distSum += cost * CSpaghettiSynConfig::m_spatialWt;
	for ( int i=0; i<size13; i++ )
	{
		matchIndices[i] = assignment[i] - 1;
	}
	return distSum;
}

void CSpaghettiSyn::InitializeOutput()
{
	m_vecOutputSequence.clear();
#if 0
	//InitializeOutputViaPoissonDiskVertically();
	Vec3f pos0 = Vec3f(0.0f, -0.5f, 0.0f);
	vector<vector<Vec3f>> vecPrs = GetConnectedInputPatchExtend(10, 20);
	for ( int i=0; i<int(vecPrs.size()); i++ )
	{
		vector<CMassSpringData> vecOutputStrand;
		CMassSpringData strand;
		strand.ResizeMassSpringData(vecPrs[i].size()+1);
		strand.GetMass(0).SetPos(pos0);
		for ( int j=1; j<strand.GetNumOfMasses(); j++ )
		{
			strand.GetMass(j).SetPos(pos0+vecPrs[i][j-1]);
		}
		vecOutputStrand.push_back(strand);
		m_vecOutputSequence.push_back(vecOutputStrand);
	}
#else
#if 1 // Initialize based on coarse motion
	MassSpringSequence& coarseStrands = m_vecCoarseSequence.back();
	vector<CMassSpringData> vecOutputStrand;
	int startIdx = 0;
	for ( int i=0; i<int(coarseStrands.size()); i++ )
	{
		CMassSpringData& strand = coarseStrands[i];
		int numOfMasses = strand.GetNumOfMasses();
		int midIdx = numOfMasses / 2;
		// Strand 1...
		int idx1 = midIdx - 1;
		int segmentLength1 = idx1;
		vector<vector<Vec3f>> vecPrs1 = GetConnectedInputPatchExtend(segmentLength1, 1);
		Vec3f p01 = strand.GetMass(idx1).GetPos();
		CMassSpringData strand1;
		strand1.ResizeMassSpringData(segmentLength1+1);
		strand1.GetMass(0).SetPos(p01);
		strand1.GetMass(0).SetFlagFixed(true);
		strand1.SetStartIdx(startIdx);
		for ( int j=1; j<strand1.GetNumOfMasses(); j++ )
		{
			strand1.GetMass(j).SetPos(p01+vecPrs1[0][j-1]);
		}
		startIdx += strand1.GetNumOfMasses();
		vecOutputStrand.push_back(strand1);
		// Strand 2...
		int idx2 = midIdx + 1;
		int segmentLength2 = numOfMasses-1-idx2;
		vector<vector<Vec3f>> vecPrs2 = GetConnectedInputPatchExtend(segmentLength2, 1);
		Vec3f p02 = strand.GetMass(idx2).GetPos();
		CMassSpringData strand2;
		strand2.ResizeMassSpringData(segmentLength2+1);
		strand2.GetMass(0).SetPos(p02);
		strand2.GetMass(0).SetFlagFixed(true);
		strand2.SetStartIdx(startIdx);
		for ( int j=1; j<strand2.GetNumOfMasses(); j++ )
		{
			strand2.GetMass(j).SetPos(p02+vecPrs2[0][j-1]);
		}
		startIdx += strand2.GetNumOfMasses();
		vecOutputStrand.push_back(strand2);
	}
	m_vecOutputStrand = vecOutputStrand;
#else
	PoissonDisk3f pd3f;
	Flt distSqrMinThresh = CSpaghettiSynConfig::m_springRadius * CSpaghettiSynConfig::m_springRadius;
	const Vec3f posMin = Vec3f(-0.1f, 0.5f,-0.1f);
	const Vec3f posMax = Vec3f( 0.1f, 0.6f, 0.1f);
	vector<Vec3f> vecOrigin = pd3f.GeneratePoints(posMin, posMax, distSqrMinThresh);
	vector<MassSpringSequence> vecStrandSequence;
	int segmentLength = 20;
	int sequenceLength = 1;
	int startIdx = 0;
	for ( int n=0; n<int(vecOrigin.size()); n++ )
	{
		Vec3f pos0 = vecOrigin[n];
		vector<vector<Vec3f>> vecPrs = GetConnectedInputPatchExtend(segmentLength, sequenceLength);
		for ( int i=0; i<sequenceLength; i++ )
		{
			vector<CMassSpringData> vecOutputStrand;
			CMassSpringData strand;
			strand.ResizeMassSpringData(segmentLength+1);
			strand.GetMass(0).SetPos(pos0);
			strand.GetMass(0).SetFlagFixed(true);
			strand.SetStartIdx(startIdx);
			for ( int j=1; j<strand.GetNumOfMasses(); j++ )
			{
				strand.GetMass(j).SetPos(pos0+vecPrs[i][j-1]);
			}
			vecOutputStrand.push_back(strand);
			vecStrandSequence.push_back(vecOutputStrand);
		}
		startIdx += int(segmentLength + 1);
	}
	for ( int i=0; i<sequenceLength; i++ )
	{
		MassSpringSequence sequence;
		for ( int n=0; n<int(vecOrigin.size()); n++ )
		{
			sequence.push_back(vecStrandSequence[n*sequenceLength+i][0]);
		}
		m_vecOutputSequence.push_back(sequence);
	}
	m_vecOutputStrand = m_vecOutputSequence[0];
	m_vecOutputSequence.clear();
#endif
#endif
#if 0 //ndef SYNTHESIZE_FRAME_BY_FRAME
	m_vecInitialSequence = m_vecOutputSequence;
	for ( int i=0; i<CSpaghettiSynConfig::m_iterationNum; i++ )
	{
		ofstream foutLog;
		foutLog.open(GetLogFileName().c_str(), std::ios_base::app);
		foutLog << "Iteration " << i+1 << "...\n";
		foutLog.close();
		UpdateOutputViaOptimization();
	}
#else
#ifndef ONLY_COARSE_MOTION
	//m_vecOutputStrand = CMassSpringData::LoadStrandsFromTXT("dumped_0001.txt");
	int interleavedSteps = CSpaghettiSynConfig::m_interleavedSteps;
	//CSpaghettiSynConfig::m_interleavedSteps = 0;
	for ( int i=0; i<100; i++ )
	{
		cout << "Init " << i << endl;
		UpdateFrameViaOptimization();
	}
#ifdef SYNTHESIZE_FRAME_BY_FRAME
	return;
#endif
	m_vecOutputStrand  = m_vecOutputSequence.back();
	m_vecOutputSequence.clear();
	for ( int i=0; i<CSpaghettiSynConfig::m_iterationNum; i++ )
	{
		cout << "Init again: " << i << endl;
		UpdateFrameViaOptimization();
	}
	vector<MassSpringSequence> vecOutputSequence = m_vecOutputSequence;
	int idxC2 = 25; //CSpaghettiSynConfig::m_iterationNum;
	int idxC1 = idxC2 - 3;
	int idxC3 = idxC2 + 10;
	int idxH1 = idxC1 - 1;
	int idxH2 = idxC1 - 10;
	Flt wtC = 1.0f;
	for ( int i=idxC1; i<idxC2; i++ )
	{
		m_vecOutputSequence[i] = EditOneOutputFrame1(m_vecOutputSequence[i], wtC);
	}
	for ( int i=idxC2; i<idxC3; i++ )
	{
		Flt wt = 1.f - (i - idxC2) / Flt(idxC3 - idxC2) * wtC;
		m_vecOutputSequence[i] = EditOneOutputFrame1(m_vecOutputSequence[i], wt);
	}
	for ( int i=idxH1; i>idxH2; i-- )
	{
		Flt wt = (i - idxH2) / Flt(idxH1 - idxH2) * wtC;
		m_vecOutputSequence[i] = EditOneOutputFrame1(m_vecOutputSequence[i], wt);
	}
	vector<MassSpringSequence> vecOutputSequence1 = m_vecOutputSequence;
	for ( int n=0; n<5; n++ )
	{
		cout << "Iter " << n << endl;
		for ( int i=idxH1; i>idxH2; i-- )
		{
			cout << "Frame " << i << endl;
			int endFrameIdx = i + 3; //min(80, i + CSpaghettiSynConfig::m_temporalWindowSize);
			UpdateOneFrameViaOptimization(i, endFrameIdx);
		}
	}
	for ( int i=0; i<int(vecOutputSequence1.size()); i++ )
	{
		m_vecOutputSequence.push_back(vecOutputSequence1[i]);
	}
	for ( int i=0; i<int(vecOutputSequence.size()); i++ )
	{
		m_vecOutputSequence.push_back(vecOutputSequence[i]);
	}
	CSpaghettiSynConfig::m_interleavedSteps = interleavedSteps;
#endif
#endif
}

void CSpaghettiSyn::InitializeOutputViaPoissonDiskVertically()
{
	PoissonDisk2f pd2f;
	Flt distSqrMinThresh = CSpaghettiSynConfig::m_springRadius * 0.25f;
	const Vec2f posMin = CSpaghettiSynConfig::m_poisson2Dmin;
	const Vec2f posMax = CSpaghettiSynConfig::m_poisson2Dmax;
	vector<Vec2f> vecOrigin = pd2f.GeneratePoints(posMin, posMax, distSqrMinThresh);
	const int numOfStrands = int(vecOrigin.size());
	cout << "Have initialized " << numOfStrands << " strands!\n";
	int numOfMasses = CSpaghettiSynConfig::m_numOfMasses;
	Flt springLength = CSpaghettiSynConfig::m_springLength;
	int startIdx = 0;
	vector<CMassSpringData> strands;
	strands.resize(numOfStrands);
	for ( int n=0; n<numOfStrands; n++ )
	{
		Vec3f pn = Vec3f(vecOrigin[n][0], 0.5f, vecOrigin[n][1]);
		strands[n].ResizeMassSpringData(numOfMasses);
		strands[n].SetStartIdx(startIdx);
		// mass 0...
		int idx0 = numOfMasses/2;
		CMass& mass0 = strands[n].GetMass(idx0);
		Vec3f pos0 = pn;
		mass0.SetPos(pos0);
		mass0.SetFlagFixed(true);
		// mass 1...
		Vec3f pr = Vec3f(1.0f,-1.0f, 0.0);
		normalize(pr);
		pr *= CSpaghettiSynConfig::m_springLength;
		int idx1 = idx0 - 1;
		CMass& mass1 = strands[n].GetMass(idx1);
		Vec3f pos1 = pos0 + pr;
		mass1.SetPos(pos1);
		mass1.SetFlagFixed(true);
		for ( int i=0; i<idx1; i++ )
		{
			CMass& mass = strands[n].GetMass(i);
			Flt dy = -i * springLength;
			mass.SetPos(pos1 + Vec3f(0.0f, dy, 0.0f));
		}
		// mass 2...
		pr[0] *= -1.0f;
		int idx2 = idx0 + 1;
		CMass& mass2 = strands[n].GetMass(idx2);
		Vec3f pos2 = pos0 + pr;
		mass2.SetPos(pos2);
		mass2.SetFlagFixed(true);
		for ( int i=idx2+1; i<numOfMasses; i++ )
		{
			CMass& mass = strands[n].GetMass(i);
			Flt dy = -(i-idx2) * springLength;
			mass.SetPos(pos2 + Vec3f(0.0f, dy, 0.0f));
		}
		startIdx += numOfMasses;
	}
	m_vecOutputStrand = strands;
	for ( int n=0; n<5; n++ )
	{
		m_vecOutputSequence.push_back(m_vecOutputStrand);
	}
}

void CSpaghettiSyn::InitializeOutputViaPoissonDiskHorizontally()
{
	PoissonDisk2f pd2f;
	Flt distSqrMinThresh = CSpaghettiSynConfig::m_springRadius * 0.25f;
	const Vec2f posMin = CSpaghettiSynConfig::m_poisson2Dmin;
	const Vec2f posMax = CSpaghettiSynConfig::m_poisson2Dmax;
	vector<Vec2f> vecOrigin = pd2f.GeneratePoints(posMin, posMax, distSqrMinThresh);
	const int numOfStrands = int(vecOrigin.size());
	cout << "Have initialized " << numOfStrands << " strands!\n";
	int numOfMasses = CSpaghettiSynConfig::m_numOfMasses;
	Flt springLength = CSpaghettiSynConfig::m_springLength;
	int startIdx = 0;
	vector<CMassSpringData> strands;
	strands.resize(numOfStrands);
	for ( int n=0; n<numOfStrands; n++ )
	{
		Vec3f pn = Vec3f(0.5f, vecOrigin[n][0], vecOrigin[n][1]);
		strands[n].ResizeMassSpringData(numOfMasses);
		strands[n].SetStartIdx(startIdx);
		for ( int i=0; i<numOfMasses; i++ )
		{
			CMass& mass = strands[n].GetMass(i);
			Flt px = -i * springLength;
			mass.SetPos(pn + Vec3f(px, 0.0f, 0.0f));
		}
		strands[n].GetMass(numOfMasses/2).SetFlagFixed(true);
		startIdx += numOfMasses;
	}
	m_vecOutputStrand = strands;
}

void CSpaghettiSyn::UpdateOutputViaOptimization()
{
	ResetCoeffMatrix();
	const int numOfFrames = int(m_vecOutputSequence.size());
	const int numOfStrands = int(m_vecOutputSequence[0].size());
	int numOfUnknownsPerFrame = 0;
	for ( int i=0; i<numOfStrands; i++ )
	{
		numOfUnknownsPerFrame += m_vecOutputSequence[0][i].GetNumOfMasses();
	}
	int numOfUnknownsTotal = numOfUnknownsPerFrame * numOfFrames;
	const int spatialSize = CSpaghettiSynConfig::m_spatialNeighSize;
	const Flt spatialWt = CSpaghettiSynConfig::m_spatialWt;
	const Flt spatialSigma = CSpaghettiSynConfig::m_spatialSigma;
	const Flt spatialNeighDist = CSpaghettiSynConfig::m_spatialNeighDist;
	const Flt temporalWt = CSpaghettiSynConfig::m_temporalWt;
	const Flt temporalSigma = CSpaghettiSynConfig::m_temporalSigma;
	vector<Flt> vecPx(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPy(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPz(numOfUnknownsTotal, 0.0f);
	m_vecCx.clear();
	m_vecCy.clear();
	m_vecCz.clear();
	m_vecCx.resize(numOfUnknownsTotal, 0.0f);
	m_vecCy.resize(numOfUnknownsTotal, 0.0f);
	m_vecCz.resize(numOfUnknownsTotal, 0.0f);
	Flt distSum = 0.0f;
#pragma omp parallel for
	for ( int d=0; d<numOfFrames; d++ )
	{
		MassSpringSequence& outputStrands = m_vecOutputSequence[d];
		for ( int n=0; n<numOfStrands; n++ )
		{
			CMassSpringData& outputStrand = outputStrands[n];
			int startIdx = outputStrand.GetStartIdx() + d * numOfUnknownsPerFrame;
			int numOfUnknownsPerStrand = outputStrands[n].GetNumOfMasses();
			for ( int i=0; i<numOfUnknownsPerStrand; i++ )
			{
				// Find the nearest input neighbor...
				MassNeighborhoodExtended neighOut = GetMassNeighborhoodExtended2(d, n, i, spatialNeighDist);
				int nearestNeighIdx = -1;
				vector<int> matchIndices;
				Flt distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
				if ( nearestNeighIdx < 0 )
				{
					neighOut = GetMassNeighborhoodExtended2(d, n, i, 0.0f);
					distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
				}
				distSum += distMin;
				MassNeighborhoodExtended& neighIn = m_vecInputNeighExtended[nearestNeighIdx];
				CMass& mass = outputStrand.GetMass(i);
				mass.SetSrcIdx(nearestNeighIdx);
				Vec3f pos = mass.GetPos();
				Vec4f quat = mass.GetQuat();
				int idxi = startIdx + i;
				vecPx[idxi] += pos[0];
				vecPy[idxi] += pos[1];
				vecPz[idxi] += pos[2];
				m_vecCx[idxi] += pos[0];
				m_vecCy[idxi] += pos[1];
				m_vecCz[idxi] += pos[2];
				// Spatial part...
				vector<Vec3f>& vecNeighPr1 = neighIn.m_vecNeighPr1;
				vector<Vec3f>& vecNeighPr2 = neighIn.m_vecNeighPr2;
				vector<Vec3f>& vecNeighPrOut1 = neighOut.m_vecNeighPr1;
				for ( int j=1; j<=int(vecNeighPrOut1.size()); j++ )
				{
					//int sizeDiff = int(vecNeighPr1.size() - vecNeighPrOut1.size());
					int idxj = idxi - j;
					Vec3f pr1 = vecNeighPr1[j-1];
					Flt wt = spatialWt * exp(spatialSigma * (j-1) * (j-1));
					UpdateCoeffMatPairVals(idxi, idxj, wt, pr1);
				}
				for ( int j=1; j<=int(vecNeighPr2.size()); j++ )
				{
					int idxj = idxi + j;
					Vec3f pr2 = vecNeighPr2[j-1];
					Flt wt = spatialWt * exp(spatialSigma * (j-1) * (j-1));
					if ( i+j < numOfUnknownsPerStrand )
					{
						UpdateCoeffMatPairVals(idxi, idxj, wt, pr2);
					}
				}
				// Temporal part...
				vector<Vec3f>& vecPr0 = neighIn.m_vecPr0;
				for ( int j=1; j<=CSpaghettiSynConfig::m_temporalWindowSize; j++ )
				{
					Flt wt = temporalWt * exp(temporalSigma * j * j);
					if ( j <= int(neighOut.m_vecPr0.size()) )
					{
						Vec3f pr = vecPr0[j-1];
						//Vec3f posPrev = m_vecOutputSequence[d-j][n].GetMass(i).GetPos();
						//Vec3f posPredict = posPrev + pr;
						//UpdateCoeffMatDiagonalVal(idxi, wt, posPredict, false);
						int idxj = idxi - j * numOfUnknownsPerFrame;
						UpdateCoeffMatPairVals(idxi, idxj, wt, pr);
					}
					else
					{
						//UpdateCoeffMatDiagonalVal(idxi, wt, pos, false);
					}
				}
				// TO-DO: Spatial part from other strands...
			}
		}
	}
	cout << "Step " << GetStepCount() << ": distSum = " << distSum << endl;
	ofstream foutLog;
	foutLog.open(GetLogFileName().c_str(), std::ios_base::app);
	foutLog << "Step " << GetStepCount() << ": distSum = " << distSum << endl;
	foutLog.close();
#ifdef USE_TAUCS_SOLVER
	CTAUCSsolver taucsSolver;
	vector<Flt> vecPxNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCx);
	vector<Flt> vecPyNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCy);
	vector<Flt> vecPzNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCz);
#else
	CMKLsolver mklSolver;
	mklSolver.SetCoeffsFromCoeffMatrix(m_ptrCoeffMatrix);
	vector<Flt> vecPxNew = mklSolver.GetSolutionViaDCG(m_vecCx, vecPx);
	vector<Flt> vecPyNew = mklSolver.GetSolutionViaDCG(m_vecCy, vecPy);
	vector<Flt> vecPzNew = mklSolver.GetSolutionViaDCG(m_vecCz, vecPz);
#endif
	for ( int d=0; d<numOfFrames; d++ )
	{
		MassSpringSequence& outputStrands = m_vecOutputSequence[d];
		for ( int n=0; n<numOfStrands; n++ )
		{
			CMassSpringData& outputStrand = outputStrands[n];
			int startIdx = outputStrand.GetStartIdx() + d * numOfUnknownsPerFrame;
			for ( int i=0; i<outputStrand.GetNumOfMasses(); i++ )
			{
				int idxi = startIdx + i;
				Vec3f pos = Vec3f(vecPxNew[idxi], vecPyNew[idxi], vecPzNew[idxi]);
				if ( outputStrand.GetMass(i).GetFlagFixed() == true )
				{
					continue;
				}
				outputStrand.GetMass(i).SetPos(pos);
			}
		}
	}
}

void CSpaghettiSyn::UpdateFrameViaOptimization()
{
#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
	ResetInputNeighborUsedCount();
#endif
	ResetFrameCoeffMatrix();
	const int numOfStrands = int(m_vecOutputStrand.size());
	int numOfUnknownsTotal = 0;
	for ( int i=0; i<numOfStrands; i++ )
	{
		numOfUnknownsTotal += m_vecOutputStrand[i].GetNumOfMasses();
	}
	const int spatialSize = CSpaghettiSynConfig::m_spatialNeighSize;
	const Flt spatialWt = CSpaghettiSynConfig::m_spatialWt;
	const Flt spatialSigma = CSpaghettiSynConfig::m_spatialSigma;
	const Flt spatialNeighDist = CSpaghettiSynConfig::m_spatialNeighDist;
	const Flt temporalWt = CSpaghettiSynConfig::m_temporalWt;
	const Flt temporalSigma = CSpaghettiSynConfig::m_temporalSigma;
	vector<Flt> vecPx(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPy(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPz(numOfUnknownsTotal, 0.0f);
	m_vecCx.clear();
	m_vecCy.clear();
	m_vecCz.clear();
	m_vecCx.resize(numOfUnknownsTotal, 0.0f);
	m_vecCy.resize(numOfUnknownsTotal, 0.0f);
	m_vecCz.resize(numOfUnknownsTotal, 0.0f);
	Flt distSum = 0.0f;
#pragma omp parallel for
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[n];
		int startIdx = outputStrand.GetStartIdx();
		Vec3f p0 = outputStrand.GetMass(0).GetPos();
		for ( int i=0; i<outputStrand.GetNumOfMasses(); i++ )
		{
			MassNeighborhoodExtended neighOut = GetMassNeighborhoodExtended2new(n, i, 0.0f); //spatialNeighDist);
			int nearestNeighIdx = -1;
			vector<int> matchIndices;
			Flt distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
			if ( nearestNeighIdx < 0 )
			{
				neighOut = GetMassNeighborhoodExtended2new(n, i, 0.0f);
				distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
			}
#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
			m_vecInputNeighExtended[nearestNeighIdx].m_neighUsedCount ++;
#endif
			distSum += distMin;
			MassNeighborhoodExtended& neighIn = m_vecInputNeighExtended[nearestNeighIdx];
			CMass& mass = outputStrand.GetMass(i);
			mass.SetSrcIdx(nearestNeighIdx);
			Vec3f pos = mass.GetPos();
			Vec4f quat = mass.GetQuat();
			mass.SetPos(pos);
			int idxi = startIdx + i;
			vecPx[idxi] += pos[0];
			vecPy[idxi] += pos[1];
			vecPz[idxi] += pos[2];
			m_vecCx[idxi] += pos[0];
			m_vecCy[idxi] += pos[1];
			m_vecCz[idxi] += pos[2];
			if ( mass.GetFlagFixed() == true )
			{
				UpdateCoeffMatDiagonalVal(idxi, 10.0f, pos);
			}
			vector<Vec3f>& vecNeighPr1 = neighIn.m_vecNeighPr1;
			vector<Vec3f>& vecNeighPr2 = neighIn.m_vecNeighPr2;
			vector<Vec3f>& vecNeighPrOut1 = neighOut.m_vecNeighPr1;
			vector<Vec3f>& vecNeighPrOut2 = neighOut.m_vecNeighPr2;
			for ( int j=1; j<=int(vecNeighPrOut1.size()); j++ )
			{
				int idxj = idxi - j;
				Vec3f pr1 = vecNeighPr1[j-1];
				Flt wt = spatialWt * exp(spatialSigma * (j-1) * (j-1));
				UpdateCoeffMatPairVals(idxi, idxj, wt, pr1);
			}
			for ( int j=1; j<=int(vecNeighPrOut2.size()); j++ )
			{
				int idxj = idxi + j;
				Vec3f pr2 = vecNeighPr2[j-1];
				Flt wt = spatialWt * exp(spatialSigma * (j-1) * (j-1));
				UpdateCoeffMatPairVals(idxi, idxj, wt, pr2);
			}
			// temporal part...
			vector<Vec3f>& vecPr0 = neighIn.m_vecPr0;
			for ( int j=1; j<=CSpaghettiSynConfig::m_temporalWindowSize; j++ )
			{
				Flt wt = temporalWt * exp(temporalSigma * j * j);
				if ( j <= int(neighOut.m_vecPr0.size()) )
				{
					Vec3f pr = vecPr0[j-1];
					Vec3f posPrev = m_vecOutputSequence[m_vecOutputSequence.size()-j][n].GetMass(i).GetPos();
					Vec3f posPredict = posPrev + pr;
					UpdateCoeffMatDiagonalVal(idxi, wt, posPredict);
				}
				else
				{
					//UpdateCoeffMatDiagonalVal(idxi, wt, pos);
				}
			}
			// Spatial part from other strands...
			neighOut = GetMassNeighborhoodExtended2new(n, i, spatialNeighDist);
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
					UpdateCoeffMatPairVals(idxi, idxj, spatialWt, pr3);
				}
			}
		}
	} // End-OMP-Parallel
	cout << "Step " << GetStepCount() << ": distSum = " << distSum << endl;
	CollisionResponse(m_vecCx, m_vecCy, m_vecCz);
	ofstream foutLog;
	foutLog.open(GetLogFileName().c_str(), std::ios_base::app);
	foutLog << "Step " << GetStepCount() << ": distSum = " << distSum << endl;;
	foutLog.close();
#ifdef USE_TAUCS_SOLVER
	CTAUCSsolver taucsSolver;
	vector<Flt> vecPxNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCx);
	vector<Flt> vecPyNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCy);
	vector<Flt> vecPzNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCz);
#else
	CMKLsolver mklSolver;
	mklSolver.SetCoeffsFromCoeffMatrix(m_ptrCoeffMatrix);
	vector<Flt> vecPxNew = mklSolver.GetSolutionViaDCG(m_vecCx, vecPx);
	vector<Flt> vecPyNew = mklSolver.GetSolutionViaDCG(m_vecCy, vecPy);
	vector<Flt> vecPzNew = mklSolver.GetSolutionViaDCG(m_vecCz, vecPz);
#endif
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[n];
		int startIdx = outputStrand.GetStartIdx();
		for ( int i=0; i<outputStrand.GetNumOfMasses(); i++ )
		{
			if ( outputStrand.GetMass(i).GetFlagFixed() == true )
			{
				continue;
			}
			int idxi = startIdx + i;
			Vec3f pos = Vec3f(vecPxNew[idxi], vecPyNew[idxi], vecPzNew[idxi]);
			outputStrand.GetMass(i).SetPos(pos);
		}
	}
#if 0 // My own interleaved solver...
	int interval = 5;
	for ( int i=0; i<5; i++ )
	{
		MassSpringSequence vecOutputStrandNew = CMassSpringData::NurbsInterpolateStrands(m_vecOutputStrand, interval);
		vecOutputStrandNew = ResolveStrandsCollision(vecOutputStrandNew, CSpaghettiSynConfig::m_springRadius * 1.1f);
		vecOutputStrandNew = CMassSpringData::DownsampleStrands(vecOutputStrandNew, interval);
		for ( int j=0; j<int(m_vecOutputStrand.size()); j++ )
		{
			for ( int i=0; i<m_vecOutputStrand[j].GetNumOfMasses(); i++ )
			{
				Vec3f pos = vecOutputStrandNew[j].GetMass(i).GetPos();
				m_vecOutputStrand[j].GetMass(i).SetPos(pos);
			}
		}
	}
#endif
	if ( CSpaghettiSynConfig::m_interleavedSteps > 0 )
	{
		CBulletSpaghetti bulletSolver;
		bulletSolver.InitBulletSpaghetti(m_vecOutputStrand);
		for ( int i=0; i<CSpaghettiSynConfig::m_interleavedSteps; i++ )
		{
			bulletSolver.clientMoveAndDisplay();
		}
		bulletSolver.EditSpaghettiFromBullet(m_vecOutputStrand);
	}
	m_vecOutputSequence.push_back(m_vecOutputStrand);
	//while ( int(m_vecOutputSequence.size()) > CSpaghettiSynConfig::m_temporalWindowSize )
	//{
	//	m_vecOutputSequence.erase(m_vecOutputSequence.begin());
	//}
}

void CSpaghettiSyn::UpdateOneFrameViaOptimization(int frameIdx, int endFrameIdx)
{
#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
	ResetInputNeighborUsedCount();
#endif
	ResetFrameCoeffMatrix();
	m_vecOutputStrand = m_vecOutputSequence[frameIdx];
	const int numOfStrands = int(m_vecOutputStrand.size());
	int numOfUnknownsTotal = 0;
	for ( int i=0; i<numOfStrands; i++ )
	{
		numOfUnknownsTotal += m_vecOutputStrand[i].GetNumOfMasses();
	}
	const int spatialSize = CSpaghettiSynConfig::m_spatialNeighSize;
	const Flt spatialWt = CSpaghettiSynConfig::m_spatialWt;
	const Flt spatialSigma = CSpaghettiSynConfig::m_spatialSigma;
	const Flt spatialNeighDist = CSpaghettiSynConfig::m_spatialNeighDist;
	const Flt temporalWt = CSpaghettiSynConfig::m_temporalWt;
	const Flt temporalSigma = CSpaghettiSynConfig::m_temporalSigma;
	vector<Flt> vecPx(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPy(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPz(numOfUnknownsTotal, 0.0f);
	m_vecCx.clear();
	m_vecCy.clear();
	m_vecCz.clear();
	m_vecCx.resize(numOfUnknownsTotal, 0.0f);
	m_vecCy.resize(numOfUnknownsTotal, 0.0f);
	m_vecCz.resize(numOfUnknownsTotal, 0.0f);
	Flt distSum = 0.0f;
	int frameShift = endFrameIdx - frameIdx;
#pragma omp parallel for
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[n];
		int startIdx = outputStrand.GetStartIdx();
		Vec3f p0 = outputStrand.GetMass(0).GetPos();
		for ( int i=0; i<outputStrand.GetNumOfMasses(); i++ )
		{
			MassNeighborhoodExtended neighOut = GetMassNeighborhoodExtended2(endFrameIdx, n, i, 0.0f); //spatialNeighDist);
			int nearestNeighIdx = -1;
			vector<int> matchIndices;
			Flt distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
			if ( nearestNeighIdx < 0 )
			{
				cout << "No matching!\n";
				neighOut = GetMassNeighborhoodExtended2new(n, i, 0.0f);
				distMin = GetNearestInputNeighborExtended(neighOut, nearestNeighIdx, matchIndices);
			}
#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
			m_vecInputNeighExtended[nearestNeighIdx].m_neighUsedCount ++;
#endif
			distSum += distMin;
			MassNeighborhoodExtended& neighIn = m_vecInputNeighExtended[nearestNeighIdx];
			CMass& mass = outputStrand.GetMass(i);
			mass.SetSrcIdx(nearestNeighIdx);
			Vec3f pos = mass.GetPos();
			Vec4f quat = mass.GetQuat();
			mass.SetPos(pos);
			int idxi = startIdx + i;
			vecPx[idxi] += pos[0];
			vecPy[idxi] += pos[1];
			vecPz[idxi] += pos[2];
			m_vecCx[idxi] += pos[0];
			m_vecCy[idxi] += pos[1];
			m_vecCz[idxi] += pos[2];
			if ( mass.GetFlagFixed() == true )
			{
				UpdateCoeffMatDiagonalVal(idxi, 10.0f, pos);
			}
			//vector<Vec3f>& vecNeighPr1 = neighIn.m_vecNeighPr1;
			//vector<Vec3f>& vecNeighPr2 = neighIn.m_vecNeighPr2;
			vector<Vec3f>& vecNeighPr1 = neighIn.m_vecNeighPr1Prev[frameShift-1];
			vector<Vec3f>& vecNeighPr2 = neighIn.m_vecNeighPr2Prev[frameShift-1];
			vector<Vec3f>& vecNeighPrOut1 = neighOut.m_vecNeighPr1;
			vector<Vec3f>& vecNeighPrOut2 = neighOut.m_vecNeighPr2;
			for ( int j=1; j<=int(vecNeighPrOut1.size()); j++ )
			{
				int idxj = idxi - j;
				Vec3f pr1 = vecNeighPr1[j-1];
				Flt wt = spatialWt * exp(spatialSigma * (j-1) * (j-1));
				//UpdateCoeffMatPairVals(idxi, idxj, wt, pr1);
			}
			for ( int j=1; j<=int(vecNeighPrOut2.size()); j++ )
			{
				int idxj = idxi + j;
				Vec3f pr2 = vecNeighPr2[j-1];
				Flt wt = spatialWt * exp(spatialSigma * (j-1) * (j-1));
				//UpdateCoeffMatPairVals(idxi, idxj, wt, pr2);
			}
			// temporal part...
			vector<Vec3f>& vecPr0 = neighIn.m_vecPr0;
			for ( int j=0; j<=CSpaghettiSynConfig::m_temporalWindowSize; j++ )
			{
				//Flt wt = temporalWt * exp(temporalSigma * j * j);
				//if ( j <= int(neighOut.m_vecPr0.size()) )
				//{
				//	Vec3f pr = vecPr0[j-1];
				//	Vec3f posPrev = m_vecOutputSequence[m_vecOutputSequence.size()-j][n].GetMass(i).GetPos();
				//	Vec3f posPredict = posPrev + pr;
				//	UpdateCoeffMatDiagonalVal(idxi, wt, posPredict);
				//}
				int dj = j - frameShift;
				Flt wt = temporalWt * exp(temporalSigma * dj * dj);
				if ( dj == 0 )
				{
					continue;
				}
				Vec3f pr;
				if ( j == 0 )
				{
					pr = -vecPr0[frameShift-1];
				}
				else
				{
					pr = vecPr0[j-1] - vecPr0[frameShift-1];
				}
				Vec3f posPrev = m_vecOutputSequence[endFrameIdx-j][n].GetMass(i).GetPos();
				Vec3f posPredict = posPrev + pr;
				UpdateCoeffMatDiagonalVal(idxi, wt, posPredict);
			}
		}
	} // End-OMP-Parallel
	cout << "Step " << GetStepCount() << ": distSum = " << distSum << endl;
	CollisionResponse(m_vecCx, m_vecCy, m_vecCz);
	ofstream foutLog;
	foutLog.open(GetLogFileName().c_str(), std::ios_base::app);
	foutLog << "Step " << GetStepCount() << ": distSum = " << distSum << endl;;
	foutLog.close();
#ifdef USE_TAUCS_SOLVER
	CTAUCSsolver taucsSolver;
	vector<Flt> vecPxNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCx);
	vector<Flt> vecPyNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCy);
	vector<Flt> vecPzNew = taucsSolver.GetSolution(m_ptrCoeffMatrix, m_vecCz);
#else
	CMKLsolver mklSolver;
	mklSolver.SetCoeffsFromCoeffMatrix(m_ptrCoeffMatrix);
	vector<Flt> vecPxNew = mklSolver.GetSolutionViaDCG(m_vecCx, vecPx);
	vector<Flt> vecPyNew = mklSolver.GetSolutionViaDCG(m_vecCy, vecPy);
	vector<Flt> vecPzNew = mklSolver.GetSolutionViaDCG(m_vecCz, vecPz);
#endif
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[n];
		int startIdx = outputStrand.GetStartIdx();
		for ( int i=0; i<outputStrand.GetNumOfMasses(); i++ )
		{
			if ( outputStrand.GetMass(i).GetFlagFixed() == true )
			{
				continue;
			}
			int idxi = startIdx + i;
			Vec3f pos = Vec3f(vecPxNew[idxi], vecPyNew[idxi], vecPzNew[idxi]);
			outputStrand.GetMass(i).SetPos(pos);
		}
	}
	if ( CSpaghettiSynConfig::m_interleavedSteps > 0 )
	{
		CBulletSpaghetti bulletSolver;
		bulletSolver.InitBulletSpaghetti(m_vecOutputStrand);
		for ( int i=0; i<CSpaghettiSynConfig::m_interleavedSteps; i++ )
		{
			bulletSolver.clientMoveAndDisplay();
		}
		bulletSolver.EditSpaghettiFromBullet(m_vecOutputStrand);
	}
	m_vecOutputSequence[frameIdx] = m_vecOutputStrand;
}

void CSpaghettiSyn::UpdateFrameViaBulletSimulation()
{
	CBulletSpaghetti bulletSolver;
	bulletSolver.InitBulletSpaghetti(m_vecOutputStrand);
	//for ( int i=0; i<CSpaghettiSynConfig::m_interleavedSteps; i++ )
	{
		bulletSolver.clientMoveAndDisplay();
	}
	bulletSolver.EditSpaghettiFromBullet(m_vecOutputStrand);
}

void CSpaghettiSyn::ResetCoeffMatrix()
{
	const int numOfFrames = int(m_vecOutputSequence.size());
	const int numOfStrands = int(m_vecOutputSequence[0].size());
	int numOfUnknownsPerFrame = 0;
	for ( int i=0; i<numOfStrands; i++ )
	{
		numOfUnknownsPerFrame += m_vecOutputSequence[0][i].GetNumOfMasses();
	}
	int numOfUnknownsTotal = numOfUnknownsPerFrame * numOfFrames;
	DELETE_OBJECT(m_ptrCoeffMatrix);
	m_ptrCoeffMatrix = new CCrossList(numOfUnknownsTotal, numOfUnknownsTotal);
	for ( int i=0; i<numOfUnknownsTotal; i++ )
	{
		m_ptrCoeffMatrix->UpdateDiagonalVal(i, 1.0f);
	}
	return;
}

void CSpaghettiSyn::ResetFrameCoeffMatrix()
{
	const int numOfStrands = int(m_vecOutputStrand.size());
	int numOfUnknownsTotal = 0;
	for ( int i=0; i<numOfStrands; i++ )
	{
		numOfUnknownsTotal += m_vecOutputStrand[i].GetNumOfMasses();
	}
	DELETE_OBJECT(m_ptrCoeffMatrix);
	m_ptrCoeffMatrix = new CCrossList(numOfUnknownsTotal, numOfUnknownsTotal);
	for ( int i=0; i<numOfUnknownsTotal; i++ )
	{
		m_ptrCoeffMatrix->UpdateDiagonalVal(i, 1.0f);
	}
}

void CSpaghettiSyn::CollisionResponse(vector<Flt>& vecCx, vector<Flt>& vecCy, vector<Flt>& vecCz)
{
	const Flt repulsionConstant = CSpaghettiSynConfig::m_springRepulsion;
	if ( repulsionConstant <= 0.0f ) return;
	const Flt r = CSpaghettiSynConfig::m_springRadius;
	const Flt r2 = r * r;
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
					if ( j1 == 0 && j2 == 0 ) continue;
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
						int idx1 = startIdx1 + j1;
						int idx2 = startIdx2 + j2;
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
	cout << "Have detected " << repulsionCount << " pairs of collision.\n";
}

MassNeighborhoodExtended CSpaghettiSyn::GetMassNeighborhoodExtended(vector<CMassSpringData>& strands, int strandIdx, int massIdx, Flt neighDist)
{
	CMassSpringData& msData = strands[strandIdx];
	const int numOfMasses = msData.GetNumOfMasses();
	const int neighSize = CSpaghettiSynConfig::m_spatialNeighSize;
	CMass& mass = msData.GetMass(massIdx);
	MassNeighborhoodExtended neigh;
	neigh.m_massIdx = massIdx;
#ifdef HACK_WITH_RELATIVE_POSITION
	neigh.m_relativePos = massIdx / Flt(numOfMasses);
#endif
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

MassNeighborhoodExtended CSpaghettiSyn::GetMassNeighborhoodExtended1(int frameIdx, int strandIdx, int massIdx, Flt neighDist, bool flagToroidal /* = false */)
{
	vector<CMassSpringData>& strands = m_vecInputSequence[frameIdx];
	MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended(strands, strandIdx, massIdx, neighDist);
	const int iMax = (flagToroidal == true) ? CSpaghettiSynConfig::m_temporalWindowSize : min(CSpaghettiSynConfig::m_temporalWindowSize, frameIdx);
	const int sequenceLength = int(m_vecInputSequence.size());
	for ( int i=1; i<=iMax; i++ )
	{
		int frameIdxPrev = frameIdx - i;
		if ( flagToroidal == true )
		{
			frameIdxPrev = (frameIdxPrev + sequenceLength) % sequenceLength;
		}
		vector<CMassSpringData>& strandsPrev = m_vecInputSequence[frameIdxPrev];
		MassNeighborhoodExtended neighPrev = GetMassNeighborhoodExtended(strandsPrev, strandIdx, massIdx, 0.0f);
		neigh.m_vecNeighPr1Prev.push_back(neighPrev.m_vecNeighPr1);
		neigh.m_vecNeighPr2Prev.push_back(neighPrev.m_vecNeighPr2);
		neigh.m_vecNeighPr3Prev.push_back(neighPrev.m_vecNeighPr3);
		neigh.m_vecPr0.push_back(neigh.m_massPos - neighPrev.m_massPos);
	}
	return neigh;
}

MassNeighborhoodExtended CSpaghettiSyn::GetMassNeighborhoodExtended2(int frameIdx, int strandIdx, int massIdx, Flt neighDist)
{
	MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended(m_vecOutputSequence[frameIdx], strandIdx, massIdx, neighDist);
	const int iMax = min(CSpaghettiSynConfig::m_temporalWindowSize, frameIdx);
	for ( int i=1; i<=iMax; i++ )
	{
		vector<CMassSpringData>& strandsPrev = m_vecOutputSequence[frameIdx-i];
		MassNeighborhoodExtended neighPrev = GetMassNeighborhoodExtended(strandsPrev, strandIdx, massIdx, 0.0f);
		neigh.m_vecNeighPr1Prev.push_back(neighPrev.m_vecNeighPr1);
		neigh.m_vecNeighPr2Prev.push_back(neighPrev.m_vecNeighPr2);
		neigh.m_vecNeighPr3Prev.push_back(neighPrev.m_vecNeighPr3);
		neigh.m_vecPr0.push_back(neigh.m_massPos - neighPrev.m_massPos);
	}
	return neigh;
}

MassNeighborhoodExtended CSpaghettiSyn::GetMassNeighborhoodExtended2new(int strandIdx, int massIdx, Flt neighDist)
{
	MassNeighborhoodExtended neigh = GetMassNeighborhoodExtended(m_vecOutputStrand, strandIdx, massIdx, neighDist);
	const int iMax = min(CSpaghettiSynConfig::m_temporalWindowSize, int(m_vecOutputSequence.size()));
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

vector<SegmentPatchExtended> CSpaghettiSyn::GetSegmentPatchesExtended(vector<MassSpringSequence>& vecSequence)
{
	vector<SegmentPatchExtended> segmentPatchesExtended;
	const int patchSpatialSize = CSpaghettiSynConfig::m_spatialNeighSize;
	const int patchTemporalSize = 10; //CSpaghettiSynConfig::m_temporalWindowSize;
	for ( int i=0; i<int(vecSequence.size())-patchTemporalSize; i++ )
	{
		for ( int k=0; k<int(vecSequence[i].size()); k++ )
		{
			for ( int j=vecSequence[i][k].GetNumOfMasses()/2+1; j<vecSequence[i][k].GetNumOfMasses()-patchSpatialSize+1; j++ )
			{
				SegmentPatchExtended patchExtended;
				for ( int di=0; di<patchTemporalSize; di++ )
				{
					vector<Vec3f> vecPos;
					for ( int dj=0; dj<patchSpatialSize; dj++ )
					{
						Vec3f pos = vecSequence[i+di][k].GetMass(j+dj).GetPos();
						vecPos.push_back(pos);
					}
					patchExtended.m_vecPosSequence.push_back(vecPos);
				}
				segmentPatchesExtended.push_back(patchExtended);
			}
		}
	}
	return segmentPatchesExtended;
}

Flt CSpaghettiSyn::SegmentPatchSpatialCompatibility(SegmentPatchExtended* ptrPatch1, SegmentPatchExtended* ptrPatch2)
{
	Flt compatibility = 0.0f;
	for ( int i=0; i<int(ptrPatch1->m_vecPosSequence.size()); i++ )
	{
		vector<Vec3f>& vecPos1 = ptrPatch1->m_vecPosSequence[i];
		vector<Vec3f>& vecPos2 = ptrPatch2->m_vecPosSequence[i];
		int patchSize1 = int(vecPos1.size());
		Vec3f pos11 = vecPos1[patchSize1-2];
		Vec3f pos12 = vecPos1[patchSize1-1];
		Vec3f pos21 = vecPos2[0];
		Vec3f pos22 = vecPos2[1];
		Vec3f pr1 = pos12 - pos11;
		Vec3f pr2 = pos22 - pos21;
		compatibility += dist2(pr1, pr2);
	}
	return compatibility;
}

Flt CSpaghettiSyn::SegmentPatchTemporalCompatibility(SegmentPatchExtended* ptrPatch1, SegmentPatchExtended* ptrPatch2)
{
	Flt compatibility = 0.0f;
	int patchSize1 = int(ptrPatch1->m_vecPosSequence.size());
	vector<Vec3f>& vecPos11 = ptrPatch1->m_vecPosSequence[patchSize1-2];
	vector<Vec3f>& vecPos12 = ptrPatch1->m_vecPosSequence[patchSize1-1];
	vector<Vec3f>& vecPos21 = ptrPatch2->m_vecPosSequence[0];
	vector<Vec3f>& vecPos22 = ptrPatch2->m_vecPosSequence[1];
	for ( int i=0; i<int(vecPos11.size()); i++ )
	{
		Vec3f pr1 = vecPos12[i] - vecPos11[i];
		Vec3f pr2 = vecPos22[i] - vecPos21[i];
		compatibility += dist2(pr1, pr2);
	}
	return compatibility;
}

vector<vector<Vec3f>> CSpaghettiSyn::GetConnectedInputPatchExtend(int segmentLength, int sequenceLength)
{
	vector<vector<Vec3f>> vecPatchPrOutput;
	int numOfInputPatches = int(m_vecInputPatchExtended.size());
	int patchIdx = int(rand() / Flt(RAND_MAX) * numOfInputPatches);
	patchIdx = patchIdx % numOfInputPatches;
	SegmentPatchExtended patch0 = m_vecInputPatchExtended[patchIdx];
	vecPatchPrOutput.resize(patch0.m_vecPosSequence.size()); // Resize output sequence
	for ( int j=0; j<int(vecPatchPrOutput.size()); j++ )
	{
		Vec3f p0 = patch0.m_vecPosSequence[j][0];
		for ( int i=1; i<int(patch0.m_vecPosSequence[j].size()); i++ )
		{
			vecPatchPrOutput[j].push_back(patch0.m_vecPosSequence[j][i] - p0);
		}
	}
	SegmentPatchExtended* ptrPatch1 = &patch0;
	vector<SegmentPatchExtended*> ptrPatches;
	ptrPatches.push_back(ptrPatch1);
	while ( int(vecPatchPrOutput[0].size()) < segmentLength )
	{
		Flt compatibility = 1e10;
		Flt compatibilityThresh = CSpaghettiSynConfig::m_compatibilityThresh;
		int trialCount = 0;
		SegmentPatchExtended patch0New;
		while ( compatibility > compatibilityThresh )
		{
			trialCount ++;
			patchIdx = int(rand() / Flt(RAND_MAX) * numOfInputPatches);
			patchIdx = patchIdx % numOfInputPatches;
			patch0New = m_vecInputPatchExtended[patchIdx];
			compatibility = SegmentPatchSpatialCompatibility(ptrPatch1, &patch0New);
			if ( trialCount > 100 )
			{
				compatibilityThresh += CSpaghettiSynConfig::m_compatibilityThresh * 0.1f;
				trialCount = 0;
			}
		}
		for ( int j=0; j<int(vecPatchPrOutput.size()); j++ )
		{
			Vec3f pend = vecPatchPrOutput[j].back();
			Vec3f p1 = patch0New.m_vecPosSequence[j][1];
			for ( int i=2; i<int(patch0New.m_vecPosSequence[j].size()); i++ )
			{
				vecPatchPrOutput[j].push_back(patch0New.m_vecPosSequence[j][i] + pend - p1);
			}

		}
		ptrPatch1 = &patch0New;
		ptrPatches.push_back(ptrPatch1);
	}
	// Temporal extension...
	while ( int(vecPatchPrOutput.size()) < sequenceLength )
	{
		Flt compatibility = 1e10;
		Flt compatibilityThresh = CSpaghettiSynConfig::m_compatibilityThresh;
		int trialCount = 0;
		SegmentPatchExtended patchNew;
		while ( compatibility > compatibilityThresh )
		{
			trialCount ++;
			patchIdx = int(rand() / Flt(RAND_MAX) * numOfInputPatches);
			patchIdx = patchIdx % numOfInputPatches;
			patchNew = m_vecInputPatchExtended[patchIdx];
			compatibility = SegmentPatchTemporalCompatibility(ptrPatches[0], &patchNew);
			if ( trialCount > 100 )
			{
				compatibilityThresh += CSpaghettiSynConfig::m_compatibilityThresh * 0.1f;
				trialCount = 0;
			}
		}
		//cout << "Number of trials: " << trialCount << endl;
		vector<vector<Vec3f>> vecPatchPrNew;
		vecPatchPrNew.resize(patchNew.m_vecPosSequence.size());
		for ( int j=0; j<int(vecPatchPrNew.size()); j++ )
		{
			Vec3f p0 = patchNew.m_vecPosSequence[j][0];
			for ( int i=1; i<int(patchNew.m_vecPosSequence[j].size()); i++ )
			{
				vecPatchPrNew[j].push_back(patchNew.m_vecPosSequence[j][i] - p0);
			}
		}
		SegmentPatchExtended* ptrPatch1 = &patchNew;
		vector<SegmentPatchExtended*> ptrPatchesNew;
		ptrPatchesNew.push_back(ptrPatch1);
		int idx = 1;
		while ( int(vecPatchPrNew[0].size()) < segmentLength )
		{
			Flt compatibility = 1e10;
			Flt compatibilityThresh = CSpaghettiSynConfig::m_compatibilityThresh;
			int trialCount = 0;
			SegmentPatchExtended patchNew;
			while ( compatibility > compatibilityThresh )
			{
				trialCount ++;
				patchIdx = int(rand() / Flt(RAND_MAX) * numOfInputPatches);
				patchIdx = patchIdx % numOfInputPatches;
				patchNew = m_vecInputPatchExtended[patchIdx];
				compatibility = SegmentPatchSpatialCompatibility(ptrPatch1, &patchNew) + SegmentPatchTemporalCompatibility(ptrPatches[idx], &patchNew);
				if ( trialCount > 100 )
				{
					compatibilityThresh += CSpaghettiSynConfig::m_compatibilityThresh * 0.1f;
					trialCount = 0;
				}
			}
			//cout << "Number of trials: " << trialCount << endl;
			for ( int j=0; j<int(vecPatchPrNew.size()); j++ )
			{
				Vec3f pend = vecPatchPrNew[j].back();
				Vec3f p1 = patchNew.m_vecPosSequence[j][1];
				for ( int i=2; i<int(patchNew.m_vecPosSequence[j].size()); i++ )
				{
					vecPatchPrNew[j].push_back(patchNew.m_vecPosSequence[j][i] + pend - p1);
				}
			}
			ptrPatch1 = &patchNew;
			ptrPatchesNew.push_back(ptrPatch1);
		}
		for ( int i=0; i<int(vecPatchPrNew.size()); i++ )
		{
			vecPatchPrOutput.push_back(vecPatchPrNew[i]); // Add to output sequence
		}
		ptrPatches = ptrPatchesNew;
		idx ++;
	}
	return vecPatchPrOutput;
}

void CSpaghettiSyn::TimeIntegration()
{
	//FinishCoeffsAccordingToBoundary();
	const float dt = CSpaghettiSynConfig::m_timeStep;
	const float springLength = CSpaghettiSynConfig::m_springLength;
	const float stiffness = CSpaghettiSynConfig::m_springStiffness;
	const float springFriction = CSpaghettiSynConfig::m_springFriction;
	const float gravity = CSpaghettiSynConfig::m_gravity;
	const float airFriction = CSpaghettiSynConfig::m_airFriction;
	const int numOfMassSprings = int(m_vecOutputStrand.size());
#pragma omp parallel for
	for ( int i=0; i<numOfMassSprings; i++ )
	{
		m_vecOutputStrand[i].ComputeForce(springLength, stiffness, springFriction, gravity, airFriction);
	}
	ApplyGroundBoundary();
	CollisionResponse();
#pragma omp parallel for
	for ( int i=0; i<numOfMassSprings; i++ )
	{
		m_vecOutputStrand[i].TimeIntegrationMassSpring(dt);
	}
}

void CSpaghettiSyn::CollisionResponse()
{
	if ( CSpaghettiSynConfig::m_springRepulsion <= 0.0f ) return;
	const int numOfMassSprings = int(m_vecOutputStrand.size());
	const float r = CSpaghettiSynConfig::m_springRadius;
	const float r2 = r * r;
	int repulsionCount = 0;
#pragma omp parallel for
	for ( int i1=0; i1<numOfMassSprings; i1++ )
	{
		CMassSpringData& strand1 = m_vecOutputStrand[i1];
		int numOfMasses1 = strand1.GetNumOfMasses();
		for ( int j1=0; j1<numOfMasses1; j1++ )
		{
			Vec3f p1 = strand1.GetMass(j1).GetPos();
			for ( int i2=i1; i2<numOfMassSprings; i2++ )
			{
				CMassSpringData& strand2 = m_vecOutputStrand[i2];
				int numOfMasses2 = strand2.GetNumOfMasses();
				for ( int j2=0; j2<numOfMasses2; j2++ )
				{
					if ( (i2 == i1) ) // No self-collision
					{
						continue;
					}
					Vec3f p2 = strand2.GetMass(j2).GetPos();
					Flt distSqr = dist2(p1, p2);
					if ( distSqr < r2 )
					{
						Vec3f replusionForce = p2 - p1;
						if ( distSqr == 0.0f )
						{
							for ( int i=0; i<3; i++ )
							{
								replusionForce[i] = rand() / Flt(RAND_MAX);
							}
						}
						Flt d = normalize(replusionForce);
						replusionForce = replusionForce * CSpaghettiSynConfig::m_springRepulsion * (r - d);
						strand1.GetMass(j1).AddForce(-replusionForce);
						strand2.GetMass(j2).AddForce(replusionForce);
						repulsionCount ++;
					}
				} // End-For-j2
			} // End-For-i2
		} // End-For-j1
	} // End-For-i1
	//cout << "Number of repulsion pairs: " << repulsionCount << endl;
}

void CSpaghettiSyn::FinishCoeffsAccordingToBoundary()
{
	const float dt = CSpaghettiSynConfig::m_timeStep;
	const int numOfMassSprings = int(m_vecOutputStrand.size());
	int n = 100;
	float nInv = 1.0 / float(n);
	float vy = (m_stepCount < n) ? 3.0 * sin(m_stepCount * nInv * M_PI) : 0.0f;
	Vec3f vel = Vec3f(0.0f, vy, 0.0f);
#pragma omp parallel for
	for ( int i=0; i<numOfMassSprings; i++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[i];
		const int numOfMasses = outputStrand.GetNumOfMasses();
		// mass 0...
		int idx0 = numOfMasses/2;
		CMass& mass0 = outputStrand.GetMass(idx0);
		Vec3f pos0 = mass0.GetPos();
		pos0 = pos0 + dt * vel;
		mass0.SetPos(pos0);
		mass0.SetVel(vel);
		mass0.SetFlagFixed(true);
		// mass 1...
		Vec3f pr = Vec3f(1.0f,-1.0f, 0.0);
		normalize(pr);
		pr *= CSpaghettiSynConfig::m_springLength;
		int idx1 = idx0 - 1;
		CMass& mass1 = outputStrand.GetMass(idx1);
		Vec3f pos1 = pos0 + pr;
		mass1.SetPos(pos1);
		mass1.SetVel(vel);
		mass1.SetFlagFixed(true);
		// mass 2...
		pr[0] *= -1.0f;
		int idx2 = idx0 + 1;
		CMass& mass2 = outputStrand.GetMass(idx2);
		Vec3f pos2 = pos0 + pr;
		mass2.SetPos(pos2);
		mass2.SetVel(vel);
		mass2.SetFlagFixed(true);
	}
}

void CSpaghettiSyn::ApplyGroundBoundary()
{
	const Flt groundHeight = CSpaghettiSynConfig::m_groundHeight;
	const int numOfMassSprings = int(m_vecOutputStrand.size());
#pragma omp parallel for
	for ( int n=0; n<numOfMassSprings; n++ )
	{
		CMassSpringData& outputStrand = m_vecOutputStrand[n];
		const int numOfMasses = outputStrand.GetNumOfMasses();
		for ( int i=0; i<numOfMasses; i++ )
		{
			CMass& mass = outputStrand.GetMass(i);
			Vec3f pi = mass.GetPos();
			if ( pi[1] < groundHeight ) // A mass collides with the ground...
			{
				// Apply the friction force to create a sliding effect...
				Vec3f vel = mass.GetVel();
				vel[1] = 0.0f;
				mass.AddForce(-vel * CSpaghettiSynConfig::m_groundFriction);
				// Apply the absorption force vertical to the ground...
				vel = mass.GetVel();
				vel[0] = 0.0f;
				vel[2] = 0.0f;
				if ( vel[1] < 0.0f )
				{
					mass.AddForce(-vel * CSpaghettiSynConfig::m_groundAbsorption);
				}
				// Apply the ground repulsion force...
				Vec3f force = Vec3f(0.0f, CSpaghettiSynConfig::m_groundRepulsion, 0.0f) * (groundHeight - pi[1]);
				mass.AddForce(force);
			}
		}
	}
}

void CSpaghettiSyn::UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos, bool flagUpdateMat /* = true */)
{
	if ( flagUpdateMat == true )
	{
		m_ptrCoeffMatrix->UpdateDiagonalVal(idx, wt);
	}
	m_vecCx[idx] += wt * pos[0];
	m_vecCy[idx] += wt * pos[1];
	m_vecCz[idx] += wt * pos[2];
}

void CSpaghettiSyn::UpdateCoeffMatPairVals(int idxi, int idxj, Flt wt, Vec3f pr, bool flagUpdateMat /* = true */)
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

void CSpaghettiSyn::ChangeInputSequence()
{
	vector<MassSpringSequence> vecInputSequenceNew;
	for ( int i=0; i<int(m_vecInputSequence.size()); i++ )
	{
		MassSpringSequence vecStrandNew;
		MassSpringSequence vecStrandOld = m_vecInputSequence[i];
		for ( int j=0; j<int(vecStrandOld.size()); j++ )
		{
			CMassSpringData& strand = vecStrandOld[j];
			int numOfMass = strand.GetNumOfMasses();
			int midIdx = numOfMass / 2;
			CMassSpringData strandNew1;
			CMassSpringData strandNew2;
			for ( int k=midIdx-1; k>=0; k-- )
			{
				strandNew1.AddMass(strand.GetMass(k));
			}
			for ( int k=midIdx+1; k<numOfMass; k++ )
			{
				strandNew2.AddMass(strand.GetMass(k));
			}
			vecStrandNew.push_back(strandNew1);
			vecStrandNew.push_back(strandNew2);
		}
		vecInputSequenceNew.push_back(vecStrandNew);
	}
	m_vecInputSequence = vecInputSequenceNew;
}

void CSpaghettiSyn::ChangeInputSequenceAccordingToCoarseMotion(Flt segmentLength)
{
	vector<MassSpringSequence> vecInputSequenceNew = m_vecInputSequence;
	int segmentCount = 0;
	Flt segmentLengthSum = 0.f;
	for ( int i=0; i<int(m_vecInputSequence.size()); i++ )
	{
		MassSpringSequence& strands = m_vecInputSequence[i];
		for ( int j=0; j<strands.size(); j++ )
		{
			CMassSpringData& strand = strands[j];
			for ( int k=1; k<strand.GetNumOfMasses(); k++ )
			{
				Vec3f p0 = strand.GetMass(k-1).GetPos();
				Vec3f p1 = strand.GetMass(k).GetPos();
				Flt segmentLengthTmp = dist(p0, p1);
				segmentLengthSum += segmentLengthTmp;
				segmentCount ++;
			}
		}
	}
	Flt segmentLengthAve = segmentLengthSum / Flt(segmentCount);
	Flt scaling = segmentLength / segmentLengthAve;
	scaling = (scaling + 1.0f) / 2.0f;
	cout << "Scale the input sequence with a factor " << scaling << "!\n";
	for ( int i=0; i<int(m_vecInputSequence.size()); i++ )
	{
		MassSpringSequence& strands = m_vecInputSequence[i];
		MassSpringSequence& strandsNew = vecInputSequenceNew[i];
		for ( int j=0; j<strands.size(); j++ )
		{
			Vec3f pLast = strands[j].GetMass(0).GetPos();
			CMassSpringData& strand = strands[j];
			for ( int k=1; k<strand.GetNumOfMasses(); k++ )
			{
				Vec3f p0 = strand.GetMass(k-1).GetPos();
				Vec3f p1 = strand.GetMass(k).GetPos();
				Vec3f pr = p1 - p0;
				pr *= scaling;
				Vec3f pk = pLast + pr;
				strandsNew[j].GetMass(k).SetPos(pk);
				pLast = pk;
			}
		}
	}
	m_vecInputSequence = vecInputSequenceNew;
}

void CSpaghettiSyn::ChangeCoarseSequence0()
{
	for ( int i=0; i<int(m_vecCoarseSequence.size()); i++ )
	{
		CMassSpringData& strand = m_vecCoarseSequence[i][0];
		int numOfMasses = strand.GetNumOfMasses();
		int idx1 = 0; //numOfMasses - 1;
		int idx2 = 1; //numOfMasses - 2;
		int idx3 = 2; //numOfMasses - 3;
		Vec3f p2 = strand.GetMass(idx2).GetPos();
		Vec3f p3 = strand.GetMass(idx3).GetPos();
		Vec3f p1 = p2 * 2.0f - p3;
		strand.GetMass(idx1).SetPos(p1);
	}
}

void CSpaghettiSyn::ChangeCoarseSequence1(int frameIdx /* = 63 */, int wd /* = 50 */, int segmentLength /* = 10 */, Vec3f pd /* = Vec3f */, int type /* = 0 */)
{
	//int frameIdx = 63;
	//int wd = 50;
	//int segmentLength = 10;
	//Vec3f pd = Vec3f(0.1f, 0.0f, 0.0f);
	for ( int i=frameIdx-wd; i<=frameIdx+wd; i++ )
	{
		Flt frameWt = 1.0f - abs(i - frameIdx) / Flt(wd);
		frameWt = pow(frameWt, 5.0f);
		frameWt = frameWt * HALF_PI;
		frameWt = sin(frameWt);
		CMassSpringData& strand = m_vecCoarseSequence[i][0];
		for ( int j=0; j<segmentLength; j++ )
		{
			Flt massWt = 1.0f - j / Flt(segmentLength);
			massWt = pow(massWt, 5.0f);
			//massWt = massWt * HALF_PI;
			//massWt = sin(massWt);
			int nj = (type == 0) ? j : strand.GetNumOfMasses() - 1 - j;
			Vec3f pdTmp = pd * frameWt * massWt;
			CMass& mass = strand.GetMass(nj);
			Vec3f pj = mass.GetPos();
			mass.SetPos(pj + pdTmp);
		}
	}
}

void CSpaghettiSyn::ChangeCoarseSequence()
{
	PoissonDisk3f pd3f;
	Flt distSqrMinThresh = CSpaghettiSynConfig::m_springRadius * CSpaghettiSynConfig::m_springRadius;
	//const Vec3f posMin = Vec3f(-0.05f,-0.05f,-0.2f);
	//const Vec3f posMax = Vec3f( 0.05f, 0.05f, 0.2f);
	//vector<Vec3f> vecOrigin = pd3f.GeneratePoints(posMin, posMax, distSqrMinThresh);
	vector<MassSpringSequence> vecCoarseSequenceNew(m_vecCoarseSequence.size());
	vector<Vec3f> vecOrigin;
	int kMax = 1; // Thickness
	Flt r = CSpaghettiSynConfig::m_springRadius;
	for ( int i=0; i<=1; i++ )
	{
		for ( int k=0; k<=kMax; k++ )
		{
			int jMin = -4;
			int jMax = (k == 0) ? 4 : 3;
			Flt jd = (k == 0) ? 0.f : r*0.5f;
			for ( int j=jMin; j<=jMax; j++ )
			{
				Vec3f pos = Vec3f(i*r*3.0f + i*0.05f + i*kMax*r*2.0f, 0, j*r*1.05f + jd);
				for ( int n=0; n<int(vecCoarseSequenceNew.size()); n++ )
				{
					CMassSpringData strand = m_vecCoarseSequence[n][0];
					strand = strand.TranslateStrand(pos);
					strand = strand.DeformStrand(k * r);
					vecCoarseSequenceNew[n].push_back(strand);
				}
				vecOrigin.push_back(pos);
			}
		}
	}
	m_vecCoarseSequence = vecCoarseSequenceNew;
}

void CSpaghettiSyn::ChangeCoarseSequenceAsNonStretchable()
{
	vector<MassSpringSequence> vecCoarseSequenceNew;
	for ( int i=0; i<int(m_vecCoarseSequence.size()); i++ )
	{
		MassSpringSequence vecStrandNew;
		CMassSpringData strand = m_vecCoarseSequence[i][0];
		CMassSpringData strandNew = strand.FixLengthAsNonStretchable();
		vecStrandNew.push_back(strandNew);
		vecCoarseSequenceNew.push_back(vecStrandNew);
	}
	m_vecCoarseSequence = vecCoarseSequenceNew;
}

void CSpaghettiSyn::DumpChopsticks(vector<MassSpringSequence>& coarseSequence)
{
	char fileName[MAX_PATH];
	for ( int n=0; n<int(coarseSequence.size()); n++ )
	{
		sprintf_s(fileName, "%schopsticks_%04d.obj", CSpaghettiSynConfig::m_outputPrefix.c_str(), n);
		DumpChopsticks(coarseSequence[n], fileName);
	}
	exit(0);
}

void CSpaghettiSyn::DumpChopsticks(MassSpringSequence& coarseStrands, const char* fileName)
{
	vector<Vec3f> vecPosCen;
	CMassSpringData& strand0 = coarseStrands[0];
	int numOfMass = strand0.GetNumOfMasses();
	Vec3f p1 = strand0.GetMass(numOfMass/2-1).GetPos();
	Vec3f p2 = strand0.GetMass(numOfMass/2+1).GetPos();
	Vec3f p0 = (p1 + p2) * 0.5f;
	Flt radius = dist(p2, p1) * 0.5f - CSpaghettiSynConfig::m_springRadius * 0.5f;
	vecPosCen.push_back(p0);
	if ( coarseStrands.size() > 1 )
	{
		CMassSpringData& strand1 = coarseStrands[coarseStrands.size()/2];
		int numOfMass = strand1.GetNumOfMasses();
		Vec3f p1 = strand1.GetMass(numOfMass/2-1).GetPos();
		Vec3f p2 = strand1.GetMass(numOfMass/2+1).GetPos();
		Vec3f p0 = (p1 + p2) * 0.5f;
		vecPosCen.push_back(p0);
	}
	vector<CMassSpringData> vecChopsticks;
	for ( int i=0; i<int(vecPosCen.size()); i++ )
	{
		CMassSpringData chopstick;
		int chopstickLength = 20;
		chopstick.ResizeMassSpringData(chopstickLength);
		Vec3f p0 = vecPosCen[i];
		for ( int j=0; j<chopstickLength; j++ )
		{
			Vec3f pj = p0 - Flt(j - 4) * Vec3f(0.f, 0.f, 0.05f);
			chopstick.GetMass(j).SetPos(pj);
		}
		vecChopsticks.push_back(chopstick);
	}
	CMassSpringData::DumpStrandsToOBJz(vecChopsticks, fileName, radius, 16, 1, 0);
}

vector<CMassSpringData> CSpaghettiSyn::EditOneOutputFrame(vector<CMassSpringData>& strands)
{
	vector<CMassSpringData> strandsNew;
	for ( int n=0; n<int(strands.size()); n++ )
	{
		CMassSpringData strand = strands[n];
		for ( int i=1; i<strand.GetNumOfMasses(); i++ )
		{
			Vec3f pLast = strands[n].GetMass(i-1).GetPos();
			Vec3f pi = strands[n].GetMass(i).GetPos();
			Vec3f pr = pi - pLast;
			Flt pr1 = pr[1];
			pr1 *= 0.95f;
			//pr1 *= 1.f - sqrt(i / Flt(strand.GetNumOfMasses())) * 0.1f;
			Vec3f prNew = Vec3f(pr[0], 0.f, pr[2]);
			if ( mag2(prNew) > 0.f )
			{
				normalize(prNew);
			}
			Flt d = sqrt(mag2(pr) - pr1 * pr1);
			prNew *= d;
			prNew[1] = pr1;
			pi = strand.GetMass(i-1).GetPos() + prNew;
			strand.GetMass(i).SetPos(pi);
		}
		strandsNew.push_back(strand);
	}
	return strandsNew;
}

vector<CMassSpringData> CSpaghettiSyn::EditOneOutputFrame1(vector<CMassSpringData>& strands, Flt wt /* = 1.0f */)
{
	int segmentLength = 10;
	Vec3f pd = Vec3f(0.15f, 0.0f, 0.0f) * wt;
	vector<CMassSpringData> strandsNew = strands;
	for ( int n=0; n<int(strandsNew.size()); n++ )
	{
		CMassSpringData& strand = strandsNew[n];
		for ( int j=0; j<segmentLength; j++ )
		{
			Flt massWt = 1.0f - j / Flt(segmentLength);
			massWt = pow(massWt, 3.0f);
			//massWt = massWt * HALF_PI;
			//massWt = sin(massWt);
			Vec3f pdTmp = pd * massWt;
			CMass& mass = strand.GetMass(strand.GetNumOfMasses()-1-j);
			Vec3f pj = mass.GetPos();
			if ( n % 2 == 0 && n >= strandsNew.size()/2 )
			{
				mass.SetPos(pj + pdTmp);
			}
			else if ( n % 2 == 1 && n < strandsNew.size()/2 )
			{
				mass.SetPos(pj - pdTmp);
			}
		}
	}
	return strandsNew;
}

vector<CMassSpringData> CSpaghettiSyn::EditCoarseStrandsFromSynthesisOutput(vector<CMassSpringData>& vecCoarseStrand)
{
	vector<CMassSpringData> vecCoarseStrandNew = vecCoarseStrand;
	vector<CMassSpringData>& vecSynthesisOutput = m_vecOutputStrand;
	Flt mixWeight = 1.0f - m_stepCount / Flt(CSpaghettiSynConfig::m_stepCountMax);
	mixWeight = max(0.0f, mixWeight);
	mixWeight = sqrt(mixWeight);
	mixWeight = 1.0f;
	for ( int i=0; i<int(vecCoarseStrandNew.size()); i++ )
	{
		CMassSpringData& coarseStrand = vecCoarseStrand[i];
		CMassSpringData& coarseStrandNew = vecCoarseStrandNew[i];
		CMassSpringData& outputStrand1 = vecSynthesisOutput[2*i];
		CMassSpringData& outputStrand2 = vecSynthesisOutput[2*i+1];
		int numOfMasses = coarseStrand.GetNumOfMasses();
		int midIdx = numOfMasses / 2;
#if 0
		Vec3f posLast = coarseStrand.GetMass(midIdx-1).GetPos();
		for ( int j=1; j<outputStrand1.GetNumOfMasses(); j++ )
		{
			Vec3f pr = outputStrand1.GetMass(j).GetPos() - outputStrand1.GetMass(j-1).GetPos();
			Flt len = dist(outputStrand1.GetMass(j).GetPos(), outputStrand1.GetMass(j-1).GetPos());
			Vec3f prOld = coarseStrand.GetMass(midIdx-j-1).GetPos() - coarseStrand.GetMass(midIdx-j).GetPos();
			Flt lenOld = mag(prOld);
			normalize(pr);
			normalize(prOld);
			Flt d = dot(prOld, pr);
			prOld *= lenOld * (1.0f - d) + len * d;
			Vec3f pj = posLast + prOld;
			coarseStrandNew.GetMass(midIdx-j-1).SetPos(pj);
			posLast = pj;
		}
		posLast = coarseStrand.GetMass(midIdx+1).GetPos();
		for ( int j=1; j<outputStrand2.GetNumOfMasses(); j++ )
		{
			Vec3f pr = outputStrand2.GetMass(j).GetPos() - outputStrand2.GetMass(j-1).GetPos();
			Flt len = dist(outputStrand2.GetMass(j).GetPos(), outputStrand2.GetMass(j-1).GetPos());
			Vec3f prOld = coarseStrand.GetMass(midIdx+j+1).GetPos() - coarseStrand.GetMass(midIdx+j).GetPos();
			Flt lenOld = mag(prOld);
			normalize(pr);
			normalize(prOld);
			Flt d = dot(prOld, pr);
			prOld *= lenOld * (1.0f - d) + len * d;
			Vec3f pj = posLast + prOld;
			coarseStrandNew.GetMass(midIdx+j+1).SetPos(pj);
			posLast = pj;
		}
#else
		Vec3f posLast = coarseStrand.GetMass(midIdx-1).GetPos();
		for ( int j=1; j<outputStrand1.GetNumOfMasses(); j++ )
		{
			Vec3f pr = outputStrand1.GetMass(j).GetPos() - outputStrand1.GetMass(j-1).GetPos();
			Vec3f prOld = coarseStrand.GetMass(midIdx-j-1).GetPos() - coarseStrand.GetMass(midIdx-j).GetPos();
			Vec3f y = prOld;
			normalize(y);
			Vec3f z = Vec3f(0.0f, 0.0f, 1.0);
			z = z - y * dot(y, z);
			normalize(z);
			Vec3f x = cross(y, z);
			Flt mat[16];
			for ( int n=0; n<16; n++ ) mat[n] = 0.0f;
			mat[15] = 1.0f;
			mat[0] = x[0]; mat[1] = x[1]; mat[2] = x[2];
			mat[4] = y[0]; mat[5] = y[1]; mat[6] = y[2];
			mat[8] = z[0]; mat[9] = z[1]; mat[10]= z[2];
			Vec4f quat;
			SetQuatFromMatrix(mat, quat);
			prOld = QuatMultiplyVec(quat, -pr);
			Vec3f pj = posLast + prOld;
			Vec3f pjOld = coarseStrandNew.GetMass(midIdx-j-1).GetPos();
			pj = pj * mixWeight + pjOld * (1.0f - mixWeight);
			coarseStrandNew.GetMass(midIdx-j-1).SetPos(pj);
			posLast = pj;
		}
		posLast = coarseStrand.GetMass(midIdx+1).GetPos();
		for ( int j=1; j<outputStrand2.GetNumOfMasses(); j++ )
		{
			Vec3f pr = outputStrand2.GetMass(j).GetPos() - outputStrand2.GetMass(j-1).GetPos();
			Vec3f prOld = coarseStrand.GetMass(midIdx+j+1).GetPos() - coarseStrand.GetMass(midIdx+j).GetPos();
			Vec3f y = prOld;
			normalize(y);
			Vec3f z = Vec3f(0.0f, 0.0f, 1.0);
			z = z - y * dot(y, z);
			normalize(z);
			Vec3f x = cross(y, z);
			Flt mat[16];
			for ( int n=0; n<16; n++ ) mat[n] = 0.0f;
			mat[15] = 1.0f;
			mat[0] = x[0]; mat[1] = x[1]; mat[2] = x[2];
			mat[4] = y[0]; mat[5] = y[1]; mat[6] = y[2];
			mat[8] = z[0]; mat[9] = z[1]; mat[10]= z[2];
			Vec4f quat;
			SetQuatFromMatrix(mat, quat);
			prOld = QuatMultiplyVec(quat, -pr);
			Vec3f pj = posLast + prOld;
			Vec3f pjOld = coarseStrandNew.GetMass(midIdx+j+1).GetPos();
			pj = pj * mixWeight + pjOld * (1.0f - mixWeight);
			coarseStrandNew.GetMass(midIdx+j+1).SetPos(pj);
			posLast = pj;
		}
#endif
	}
	return vecCoarseStrandNew;
}

vector<CMassSpringData> CSpaghettiSyn::ResolveStrandsCollision(vector<CMassSpringData>& strands, Flt rad)
{
	int numOfUnknownsTotal = 0;
	int numOfStrands = int(strands.size());
	for ( int i=0; i<numOfStrands; i++ )
	{
		strands[i].SetStartIdx(numOfUnknownsTotal);
		numOfUnknownsTotal += strands[i].GetNumOfMasses();
	}
	vector<CMassSpringData> strandsNew = strands;
	DELETE_OBJECT(m_ptrCoeffMatrix);
	m_ptrCoeffMatrix = new CCrossList(numOfUnknownsTotal, numOfUnknownsTotal);
	for ( int i=0; i<numOfUnknownsTotal; i++ )
	{
		m_ptrCoeffMatrix->UpdateDiagonalVal(i, 1.0f);
	}
	vector<Flt> vecPx(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPy(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPz(numOfUnknownsTotal, 0.0f);
	m_vecCx.clear();
	m_vecCy.clear();
	m_vecCz.clear();
	m_vecCx.resize(numOfUnknownsTotal, 0.0f);
	m_vecCy.resize(numOfUnknownsTotal, 0.0f);
	m_vecCz.resize(numOfUnknownsTotal, 0.0f);
#pragma omp parallel for
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& strand = strands[n];
		int startIdx = strand.GetStartIdx();
		for ( int i=0; i<strand.GetNumOfMasses(); i++ )
		{
			Vec3f pos = strand.GetMass(i).GetPos();
			int idxi = startIdx + i;
			vecPx[idxi] += pos[0];
			vecPy[idxi] += pos[1];
			vecPz[idxi] += pos[2];
			m_vecCx[idxi] += pos[0];
			m_vecCy[idxi] += pos[1];
			m_vecCz[idxi] += pos[2];
		}
	}
	const Flt r2 = rad * rad;
	const Flt repulsionConstant = CSpaghettiSynConfig::m_springRepulsion;
	int repulsionCount = 0;
#pragma omp parallel for
	for ( int i1=0; i1<numOfStrands; i1++ )
	{
		CMassSpringData& strand1 = strands[i1];
		int startIdx1 = strand1.GetStartIdx();
		for ( int i2=i1+1; i2<numOfStrands; i2++ )
		{
			CMassSpringData& strand2 = strands[i2];
			int startIdx2 = strand2.GetStartIdx();
			for ( int j1=0; j1<strand1.GetNumOfMasses(); j1++ )
			{
				Vec3f pos1 = strand1.GetMass(j1).GetPos();
				for ( int j2=0; j2<strand2.GetNumOfMasses(); j2++ )
				{
					if ( j1 == 0 && j2 == 0 ) continue;
					Vec3f pos2 = strand2.GetMass(j2).GetPos();
					Flt distSqr = dist2(pos1, pos2);
					if ( distSqr < r2 )
					{
						repulsionCount ++;
						Vec3f repulsionPr = pos1 - pos2;
						while ( distSqr == 0.0f )
						{
							for ( int i=0; i<3; i++ )
							{
								repulsionPr[i] = rand() / Flt(RAND_MAX);
							}
							distSqr = mag2(repulsionPr);
						}
						normalize(repulsionPr);
						repulsionPr = repulsionPr * rad;
						int idx1 = startIdx1 + j1;
						int idx2 = startIdx2 + j2;
						if ( j1 == 0 )
						{
							m_ptrCoeffMatrix->UpdateDiagonalVal(idx2, repulsionConstant);
							Vec3f prNew = pos1 - repulsionPr;
							m_vecCx[idx2] += repulsionConstant * prNew[0];
							m_vecCy[idx2] += repulsionConstant * prNew[1];
							m_vecCz[idx2] += repulsionConstant * prNew[2];
						}
						else if ( j2 == 0 )
						{
							m_ptrCoeffMatrix->UpdateDiagonalVal(idx1, repulsionConstant);
							Vec3f prNew = pos2 + repulsionPr;
							m_vecCx[idx1] += repulsionConstant * prNew[0];
							m_vecCy[idx1] += repulsionConstant * prNew[1];
							m_vecCz[idx1] += repulsionConstant * prNew[2];
						}
						else
						{
							m_ptrCoeffMatrix->UpdatePairVals(idx1, idx2, repulsionConstant);
							m_vecCx[idx1] += repulsionConstant * repulsionPr[0];
							m_vecCy[idx1] += repulsionConstant * repulsionPr[1];
							m_vecCz[idx1] += repulsionConstant * repulsionPr[2];
							m_vecCx[idx2] -= repulsionConstant * repulsionPr[0];
							m_vecCy[idx2] -= repulsionConstant * repulsionPr[1];
							m_vecCz[idx2] -= repulsionConstant * repulsionPr[2];
						}
					}
				} // End-For-j2
			} // End-For-i2
		} // End-For-j1
	} // End-For-i1
	cout << "Have detected " << repulsionCount << " pairs of collision.\n";
	CMKLsolver mklSolver;
	mklSolver.SetCoeffsFromCoeffMatrix(m_ptrCoeffMatrix);
	vector<Flt> vecPxNew = mklSolver.GetSolutionViaDCG(m_vecCx, vecPx);
	vector<Flt> vecPyNew = mklSolver.GetSolutionViaDCG(m_vecCy, vecPy);
	vector<Flt> vecPzNew = mklSolver.GetSolutionViaDCG(m_vecCz, vecPz);
#pragma omp parallel for
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& strand = strandsNew[n];
		int startIdx = strand.GetStartIdx();
		for ( int i=0; i<strand.GetNumOfMasses(); i++ )
		{
			int idxi = startIdx + i;
			Vec3f pos = Vec3f(vecPxNew[idxi], vecPyNew[idxi], vecPzNew[idxi]);
			strand.GetMass(i).SetPos(pos);
		}
	}
	return strandsNew;
}

vector<CMassSpringData> CSpaghettiSyn::ResolveStrandsCollisionNew(vector<CMassSpringData>& strands, Flt rad)
{
	int numOfUnknownsTotal = 0;
	int numOfStrands = int(strands.size());
	for ( int i=0; i<numOfStrands; i++ )
	{
		strands[i].SetStartIdx(numOfUnknownsTotal);
		numOfUnknownsTotal += strands[i].GetNumOfMasses();
	}
	vector<CMassSpringData> strandsNew = strands;
	DELETE_OBJECT(m_ptrCoeffMatrix);
	m_ptrCoeffMatrix = new CCrossList(numOfUnknownsTotal, numOfUnknownsTotal);
	vector<Flt> vecPx(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPy(numOfUnknownsTotal, 0.0f);
	vector<Flt> vecPz(numOfUnknownsTotal, 0.0f);
	m_vecCx.clear();
	m_vecCy.clear();
	m_vecCz.clear();
	m_vecCx.resize(numOfUnknownsTotal, 0.0f);
	m_vecCy.resize(numOfUnknownsTotal, 0.0f);
	m_vecCz.resize(numOfUnknownsTotal, 0.0f);
#pragma omp parallel for
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& strand = strands[n];
		int startIdx = strand.GetStartIdx();
		m_ptrCoeffMatrix->UpdateDiagonalVal(startIdx, 1.0f);
		Vec3f pos = strand.GetMass(0).GetPos();
		vecPx[startIdx] += pos[0];
		vecPy[startIdx] += pos[1];
		vecPz[startIdx] += pos[2];
		m_vecCx[startIdx] += pos[0];
		m_vecCy[startIdx] += pos[1];
		m_vecCz[startIdx] += pos[2];
		for ( int i=1; i<strand.GetNumOfMasses(); i++ )
		{
			Vec3f pr = strand.GetMass(i-1).GetPos() - strand.GetMass(i).GetPos();
			int idx1 = startIdx + i - 1;
			int idx2 = startIdx + i;
			Flt wt = 0.1f;
			m_ptrCoeffMatrix->UpdatePairVals(idx1, idx2, wt);
			m_vecCx[idx1] += wt * pr[0];
			m_vecCy[idx1] += wt * pr[1];
			m_vecCz[idx1] += wt * pr[2];
			m_vecCx[idx2] -= wt * pr[0];
			m_vecCy[idx2] -= wt * pr[1];
			m_vecCz[idx2] -= wt * pr[2];
		}
	}
	const Flt r2 = rad * rad;
	const Flt repulsionConstant = CSpaghettiSynConfig::m_springRepulsion;
	int repulsionCount = 0;
#pragma omp parallel for
	for ( int i1=0; i1<numOfStrands; i1++ )
	{
		CMassSpringData& strand1 = strands[i1];
		int startIdx1 = strand1.GetStartIdx();
		for ( int i2=i1+1; i2<numOfStrands; i2++ )
		{
			CMassSpringData& strand2 = strands[i2];
			int startIdx2 = strand2.GetStartIdx();
			for ( int j1=0; j1<strand1.GetNumOfMasses(); j1++ )
			{
				Vec3f pos1 = strand1.GetMass(j1).GetPos();
				for ( int j2=0; j2<strand2.GetNumOfMasses(); j2++ )
				{
					if ( j1 == 0 && j2 == 0 ) continue;
					Vec3f pos2 = strand2.GetMass(j2).GetPos();
					Flt distSqr = dist2(pos1, pos2);
					if ( distSqr < r2 )
					{
						repulsionCount ++;
						Vec3f repulsionPr = pos1 - pos2;
						if ( distSqr == 0.0f )
						{
							for ( int i=0; i<3; i++ )
							{
								repulsionPr[i] = rand() / Flt(RAND_MAX);
							}
						}
						normalize(repulsionPr);
						repulsionPr = repulsionPr * rad;
						int idx1 = startIdx1 + j1;
						int idx2 = startIdx2 + j2;
						if ( j1 == 0 )
						{
							m_ptrCoeffMatrix->UpdateDiagonalVal(idx2, repulsionConstant);
							Vec3f prNew = pos1 - repulsionPr;
							m_vecCx[idx2] += repulsionConstant * prNew[0];
							m_vecCy[idx2] += repulsionConstant * prNew[1];
							m_vecCz[idx2] += repulsionConstant * prNew[2];
						}
						else if ( j2 == 0 )
						{
							m_ptrCoeffMatrix->UpdateDiagonalVal(idx1, repulsionConstant);
							Vec3f prNew = pos2 + repulsionPr;
							m_vecCx[idx1] += repulsionConstant * prNew[0];
							m_vecCy[idx1] += repulsionConstant * prNew[1];
							m_vecCz[idx1] += repulsionConstant * prNew[2];
						}
						else
						{
							m_ptrCoeffMatrix->UpdatePairVals(idx1, idx2, repulsionConstant);
							m_vecCx[idx1] += repulsionConstant * repulsionPr[0];
							m_vecCy[idx1] += repulsionConstant * repulsionPr[1];
							m_vecCz[idx1] += repulsionConstant * repulsionPr[2];
							m_vecCx[idx2] -= repulsionConstant * repulsionPr[0];
							m_vecCy[idx2] -= repulsionConstant * repulsionPr[1];
							m_vecCz[idx2] -= repulsionConstant * repulsionPr[2];
						}
					}
				} // End-For-j2
			} // End-For-i2
		} // End-For-j1
	} // End-For-i1
	cout << "Have detected " << repulsionCount << " pairs of collision.\n";
	CMKLsolver mklSolver;
	mklSolver.SetCoeffsFromCoeffMatrix(m_ptrCoeffMatrix);
	vector<Flt> vecPxNew = mklSolver.GetSolutionViaDCG(m_vecCx, vecPx);
	vector<Flt> vecPyNew = mklSolver.GetSolutionViaDCG(m_vecCy, vecPy);
	vector<Flt> vecPzNew = mklSolver.GetSolutionViaDCG(m_vecCz, vecPz);
#pragma omp parallel for
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& strand = strandsNew[n];
		int startIdx = strand.GetStartIdx();
		for ( int i=0; i<strand.GetNumOfMasses(); i++ )
		{
			int idxi = startIdx + i;
			Vec3f pos = Vec3f(vecPxNew[idxi], vecPyNew[idxi], vecPzNew[idxi]);
			strand.GetMass(i).SetPos(pos);
		}
	}
	return strandsNew;
}

#ifdef HACK_WITH_NEIGHBOR_USER_COUNT
void CSpaghettiSyn::ResetInputNeighborUsedCount()
{
	for ( int i=0; i<int(m_vecInputNeighExtended.size()); i++ )
	{
		m_vecInputNeighExtended[i].m_neighUsedCount = 0;
	}
}
#endif
