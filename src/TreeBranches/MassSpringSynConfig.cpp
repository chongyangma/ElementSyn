
#include "MassSpringSynConfig.h"

bool CMassSpringSynConfig::m_flagTemporallyToroidal = true;
int CMassSpringSynConfig::m_numOfInputFrames = 100;
int CMassSpringSynConfig::m_numOfMasses = 20;
int CMassSpringSynConfig::m_numOfStrands = 1;
int CMassSpringSynConfig::m_numOfInputStrands = 1e4;
int CMassSpringSynConfig::m_neighSize = 1;
int CMassSpringSynConfig::m_patchSize = 10;
int CMassSpringSynConfig::m_inputLoadInterval = 1;
int CMassSpringSynConfig::m_inputLoadStart = 1;
int CMassSpringSynConfig::m_stepCountMax = 2000;
int CMassSpringSynConfig::m_synthesisWindowSize = 0;
string CMassSpringSynConfig::m_inputPrefix;
string CMassSpringSynConfig::m_outputPrefix;
string CMassSpringSynConfig::m_exemplarName;
vector<string> CMassSpringSynConfig::m_vecInputPrefix;
vector<Flt> CMassSpringSynConfig::m_vecInputWt;
Flt CMassSpringSynConfig::m_springLength = 0.05f;
Flt CMassSpringSynConfig::m_wtIdxDiff = 0.001f;
Flt CMassSpringSynConfig::m_wtPrDiff = 1.0f;
Flt CMassSpringSynConfig::m_wtNeighPrDiff = 1.0f;
Flt CMassSpringSynConfig::m_shapeWt = 0.0f;
Flt CMassSpringSynConfig::m_smoothWt = 1.0f;
Flt CMassSpringSynConfig::m_alignWt = 0.0f;
Flt CMassSpringSynConfig::m_temporalWt = 0.0f;
Flt CMassSpringSynConfig::m_inputWtSumInv = 1.0f;
Flt CMassSpringSynConfig::m_strandRadius = 0.05f;
Flt CMassSpringSynConfig::m_strandRepulsionConstant = 0.0f;
Flt CMassSpringSynConfig::m_neighDist = 0.0f;
Flt CMassSpringSynConfig::m_gaussianSigma = 0.0f;
Flt CMassSpringSynConfig::m_temporalSigma = 0.0f;
Vec2f CMassSpringSynConfig::m_poisson2Dmin = Vec2f(-0.1f,-0.1f);
Vec2f CMassSpringSynConfig::m_poisson2Dmax = Vec2f( 0.1f, 0.1f);

CMassSpringSynConfig::CMassSpringSynConfig()
{
	LoadFromMassSpringSynConfig();
	ResetConfig();
}

bool CMassSpringSynConfig::ReloadConfigFromFile(string fileName)
{
	bool flag = LoadFromMassSpringSynConfig(fileName);
	if ( flag == false ) return false;
	ResetConfig();
	return true;
}

bool CMassSpringSynConfig::LoadFromMassSpringSynConfig(string fileName /* = */ )
{
	CSynConfigBase::LoadFromSynConfig(fileName);
	ifstream fin(fileName.c_str());
	if ( fin.fail() == true )
	{
		cout << "Failed to load parameters from config file " << fileName << "!\n";
		return false;
	}
	string param;
	while ( fin >> param )
	{
		if ( param == string("FLAG_TEMPORALLY_TOROIDAL") )
		{
			fin >> m_flagTemporallyToroidal;
		}
		else if ( param == string("NUM_OF_INPUT_FRAMES") )
		{
			fin >> m_numOfInputFrames;
		}
		else if ( param == string("NUMBER_OF_MASS") )
		{
			fin >> m_numOfMasses;
		}
		else if ( param == string("NUMBER_OF_STRAND") )
		{
			fin >> m_numOfStrands;
		}
		else if ( param == string("NEIGHBOR_SIZE") )
		{
			fin >> m_neighSize;
		}
		else if ( param == string("INPUT_LOAD_INTERVAL") )
		{
			fin >> m_inputLoadInterval;
		}
		else if ( param == string("INPUT_LOAD_START") )
		{
			fin >> m_inputLoadStart;
		}
		else if ( param == string("NUMBER_OF_INPUT_STRANDS") )
		{
			fin >> m_numOfInputStrands;
		}
		else if ( param == string("STEP_COUNT_MAX") )
		{
			fin >> m_stepCountMax;
		}
		else if ( param == string("SYNTHESIS_WINDOW_SIZE") )
		{
			fin >> m_synthesisWindowSize;
		}
		else if ( param == string("EXEMPLAR_NAME") )
		{
			fin >> m_exemplarName;
		}
		else if ( param == string("INPUT_WEIGHT") )
		{
			Flt inputWt;
			fin >> inputWt;
			m_vecInputWt.push_back(inputWt);
		}
		else if ( param == string("SPRING_LENGTH") )
		{
			fin >> m_springLength;
		}
		else if ( param == string("WEIGHT_INDEX_DIFFERENCE") )
		{
			fin >> m_wtIdxDiff;
		}
		else if ( param == string("WEIGHT_PR_DIFFERENCE") )
		{
			fin >> m_wtPrDiff;
		}
		else if ( param == string("WEIGHT_NEIGHBOR_PR_DIFFERENCE") )
		{
			fin >> m_wtNeighPrDiff;
		}
		else if ( param == string("SHAPE_WEIGHT") )
		{
			fin >> m_shapeWt;
		}
		else if ( param == string("SMOOTH_WEIGHT") )
		{
			fin >> m_smoothWt;
		}
		else if ( param == string("ALIGN_WEIGHT") )
		{
			fin >> m_alignWt;
		}
		else if ( param == string("TEMPORAL_WEIGHT") )
		{
			fin >> m_temporalWt;
		}
		else if ( param == string("STRAND_RADIUS") )
		{
			fin >> m_strandRadius;
		}
		else if ( param == string("STRAND_REPULSION_CONST") )
		{
			fin >> m_strandRepulsionConstant;
		}
		else if ( param == string("NEIGHBOR_DISTANCE") )
		{
			fin >> m_neighDist;
		}
		else if ( param == string("GAUSSIAN_SIGMA") )
		{
			fin >> m_gaussianSigma;
		}
		else if ( param == string("TEMPORAL_SIGMA") )
		{
			fin >> m_temporalSigma;
		}
		else if ( param == string("POISSON_2D_MIN") )
		{
			fin >> m_poisson2Dmin[0] >> m_poisson2Dmin[1];
		}
		else if ( param == string("POISSON_2D_MAX") )
		{
			fin >> m_poisson2Dmax[0] >> m_poisson2Dmax[1];
		}
	}
	if ( (m_vecInputPrefix.size() > 0) )
	{
		if ( m_vecInputWt.size() != m_vecInputPrefix.size() )
		{
			m_vecInputWt.resize(m_vecInputPrefix.size(), 1.0f);
		}
		Flt inputWtSum = 0.0f;
		for ( int i=0; i<int(m_vecInputWt.size()); i++ )
		{
			inputWtSum += m_vecInputWt[i];
		}
		m_inputWtSumInv = 1.0f / inputWtSum;
	}
	return true;
}

bool CMassSpringSynConfig::DumpToMassSpringSynConfig()
{
	ostringstream oss;
	oss << m_outputPrefix << "MassSpringSynConfig.txt";
	string fileName = oss.str();
	FILE* file;
	file = fopen(fileName.c_str(), "w");
	if ( !file )
	{
		cout << "Failed to dump parameters into config file " << fileName << "!\n";
		return false;
	}
	DumpTimeAndDate(file);
	DumpParameters(file);
	fclose(file);

	return true;
}

string CMassSpringSynConfig::AddOutputPrefix(string str)
{
	if ( m_outputPrefix.empty() == true )
	{
		return str;
	}
	ostringstream oss;
	oss << m_outputPrefix << str;
	string strNew = oss.str();
	return strNew;
}

void CMassSpringSynConfig::ResetConfig()
{
	if ( m_flagRandomness == true )
	{
		srand((unsigned int)time(0));
	}
	else
	{
		srand(0); // reset rand seed!
	}
	UpdateOutputPrefix();
	DumpToMassSpringSynConfig();
}

void CMassSpringSynConfig::DumpTimeAndDate(FILE* file)
{
	time_t myTime = time(NULL);
	tm* ptrTime = localtime(&myTime);
	fprintf(file, "%02d:%02d:%02d ", ptrTime->tm_hour, ptrTime->tm_min, ptrTime->tm_sec);
	fprintf(file, "%02d/%02d/%04d\n\n", ptrTime->tm_mon+1, ptrTime->tm_mday, ptrTime->tm_year+1900);
}

void CMassSpringSynConfig::DumpParameters(FILE* file)
{
	CSynConfigBase::DumpParameters(file);
	fprintf(file, "%s\t%d\n", "FLAG_TEMPORALLY_TOROIDAL", m_flagTemporallyToroidal);
	fprintf(file, "%s\t%d\n", "NUM_OF_INPUT_FRAMES", m_numOfInputFrames);
	fprintf(file, "%s\t%d\n", "NUMBER_OF_MASS", m_numOfMasses);
	fprintf(file, "%s\t%d\n", "NUMBER_OF_STRAND", m_numOfStrands);
	fprintf(file, "%s\t%d\n", "NUMBER_OF_INPUT_STRANDS", m_numOfInputStrands);
	fprintf(file, "%s\t%d\n", "NEIGHBOR_SIZE", m_neighSize);
	fprintf(file, "%s\t%d\n", "INPUT_LOAD_INTERVAL", m_inputLoadInterval);
	fprintf(file, "%s\t%d\n", "INPUT_LOAD_START", m_inputLoadStart);
	fprintf(file, "%s\t%d\n", "STEP_COUNT_MAX", m_stepCountMax);
	fprintf(file, "%s\t%d\n", "SYNTHESIS_WINDOW_SIZE", m_synthesisWindowSize);
	DumpStringParam(file, "EXEMPLAR_NAME", m_exemplarName);
	fprintf(file, "%s\t%f\n", "SPRING_LENGTH", m_springLength);
	fprintf(file, "%s\t%f\n", "WEIGHT_INDEX_DIFFERENCE", m_wtIdxDiff);
	fprintf(file, "%s\t%f\n", "WEIGHT_PR_DIFFERENCE", m_wtPrDiff);
	fprintf(file, "%s\t%f\n", "WEIGHT_NEIGHBOR_PR_DIFFERENCE", m_wtNeighPrDiff);
	fprintf(file, "%s\t%f\n", "SHAPE_WEIGHT", m_shapeWt);
	fprintf(file, "%s\t%f\n", "SMOOTH_WEIGHT", m_smoothWt);
	fprintf(file, "%s\t%f\n", "ALIGN_WEIGHT", m_alignWt);
	fprintf(file, "%s\t%f\n", "TEMPORAL_WEIGHT", m_temporalWt);
	fprintf(file, "%s\t%f\n", "STRAND_RADIUS", m_strandRadius);
	fprintf(file, "%s\t%f\n", "STRAND_REPULSION_CONST", m_strandRepulsionConstant);
	fprintf(file, "%s\t%f\n", "NEIGHBOR_DISTANCE", m_neighDist);
	fprintf(file, "%s\t%f\n", "GAUSSIAN_SIGMA", m_gaussianSigma);
	fprintf(file, "%s\t%f\n", "TEMPORAL_SIGMA", m_temporalSigma);
	fprintf(file, "%s\t%f %f\n", "POISSON_2D_MIN", m_poisson2Dmin[0], m_poisson2Dmin[1]);
	fprintf(file, "%s\t%f %f\n", "POISSON_2D_MAX", m_poisson2Dmax[0], m_poisson2Dmax[1]);
}

void CMassSpringSynConfig::UpdateOutputPrefix()
{
	ostringstream oss;
	BOOL flag = CreateDirectory(m_outputPrefix.c_str(), NULL);
	const int numLength = 2;
	while ( flag == FALSE && m_outputPrefix.size() > 0 )
	{
		cout << "Failed to create the directory: " << m_outputPrefix.c_str() << endl;
		string subStr = m_outputPrefix.substr(m_outputPrefix.length()-numLength-1, m_outputPrefix.length()-2);
		int num = atoi(subStr.c_str()) + 1;
		char numChar[MAX_PATH];
		sprintf_s(numChar, "%02d", num);
		ostringstream oss;
		oss << m_outputPrefix.substr(0, m_outputPrefix.length()-numLength-1) << numChar << "\\";
		m_outputPrefix = oss.str();
		flag = CreateDirectory(m_outputPrefix.c_str(), NULL);
	}
	ostringstream oss2;
	oss2 << CMassSpringSynConfig::m_outputPrefix << "Snapshot";
	CreateDirectory(oss2.str().c_str(), NULL);
	ostringstream oss3;
	oss3 << CMassSpringSynConfig::m_outputPrefix << "Dumped";
	CreateDirectory(oss3.str().c_str(), NULL);
}
