
#include "SpaghettiSynConfig.h"

bool CSpaghettiSynConfig::m_flagTemporallyToroidal = true;
int CSpaghettiSynConfig::m_numOfMasses = 20;
int CSpaghettiSynConfig::m_spatialNeighSize = 1;
int CSpaghettiSynConfig::m_temporalWindowSize = 1;
int CSpaghettiSynConfig::m_iterationNum = 1;
int CSpaghettiSynConfig::m_inputLoadInterval = 1;
int CSpaghettiSynConfig::m_inputLoadStart = 1;
int CSpaghettiSynConfig::m_inputLoadEnd = 1000;
int CSpaghettiSynConfig::m_coarseLoadInterval = 1;
int CSpaghettiSynConfig::m_coarseLoadStart = 1;
int CSpaghettiSynConfig::m_coarseLoadEnd = 1000;
int CSpaghettiSynConfig::m_interleavedSteps = 0;
Flt CSpaghettiSynConfig::m_springLength = 0.05f;
Flt CSpaghettiSynConfig::m_springRadius = 0.02f;
Flt CSpaghettiSynConfig::m_springStiffness = 0.05f;
Flt CSpaghettiSynConfig::m_springFriction = 0.2f;
Flt CSpaghettiSynConfig::m_springRepulsion = 0.0f;
Flt CSpaghettiSynConfig::m_gravity = -9.81f;
Flt CSpaghettiSynConfig::m_airFriction = 0.02f;
Flt CSpaghettiSynConfig::m_groundHeight = -0.5f;
Flt CSpaghettiSynConfig::m_groundRepulsion = 100.0f;
Flt CSpaghettiSynConfig::m_groundFriction = 0.2f;
Flt CSpaghettiSynConfig::m_groundAbsorption = 20.0f;
Flt CSpaghettiSynConfig::m_spatialWt = 0.0f;
Flt CSpaghettiSynConfig::m_spatialSigma = 0.0f;
Flt CSpaghettiSynConfig::m_spatialNeighDist = 0.0f;
Flt CSpaghettiSynConfig::m_temporalWt = 0.0f;
Flt CSpaghettiSynConfig::m_temporalSigma = 0.0f;
Flt CSpaghettiSynConfig::m_compatibilityThresh = 0.01f;
Vec2f CSpaghettiSynConfig::m_poisson2Dmin = Vec2f(-0.1f,-0.1f);
Vec2f CSpaghettiSynConfig::m_poisson2Dmax = Vec2f( 0.1f, 0.1f);

CSpaghettiSynConfig::CSpaghettiSynConfig()
{
	LoadFromSpaghettiSynConfig();
	ResetConfig();
}

bool CSpaghettiSynConfig::ReloadConfigFromFile(string fileName)
{
	bool flag = LoadFromSpaghettiSynConfig(fileName);
	if ( flag == false ) return false;
	ResetConfig();
	return true;
}

bool CSpaghettiSynConfig::LoadFromSpaghettiSynConfig(string fileName /* = */ )
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
		if ( param == string("NUMBER_OF_MASS") )
		{
			fin >> m_numOfMasses;
		}
		else if ( param == string("SPATIAL_NEIGHBOR_SIZE") )
		{
			fin >> m_spatialNeighSize;
		}
		else if ( param == string("TEMPORAL_WINDOW_SIZE") )
		{
			fin >> m_temporalWindowSize;
		}
		else if ( param == string("ITERATION_NUMBER") )
		{
			fin >> m_iterationNum;
		}
		else if ( param == string("INPUT_LOAD_INTERVAL") )
		{
			fin >> m_inputLoadInterval;
		}
		else if ( param == string("INPUT_LOAD_START") )
		{
			fin >> m_inputLoadStart;
		}
		else if ( param == string("INPUT_LOAD_END") )
		{
			fin >> m_inputLoadEnd;
		}
		else if ( param == string("COARSE_LOAD_INTERVAL") )
		{
			fin >> m_coarseLoadInterval;
		}
		else if ( param == string("COARSE_LOAD_START") )
		{
			fin >> m_coarseLoadStart;
		}
		else if ( param == string("COARSE_LOAD_END") )
		{
			fin >> m_coarseLoadEnd;
		}
		else if ( param == string("INTERLEAVED_STEPS") )
		{
			fin >> m_interleavedSteps;
		}
		else if ( param == string("SPRING_LENGTH") )
		{
			fin >> m_springLength;
		}
		else if ( param == string("SPRING_RADIUS") )
		{
			fin >> m_springRadius;
		}
		else if ( param == string("SPRING_STIFFNESS") )
		{
			fin >> m_springStiffness;
		}
		else if ( param == string("SPRING_FRICTION") )
		{
			fin >> m_springFriction;
		}
		else if ( param == string("SPRING_REPULSION") )
		{
			fin >> m_springRepulsion;
		}
		else if ( param == string("GRAVITY") )
		{
			fin >> m_gravity;
		}
		else if ( param == string("AIR_FRICTION") )
		{
			fin >> m_airFriction;
		}
		else if ( param == string("GROUND_HEIGHT") )
		{
			fin >> m_groundHeight;
		}
		else if ( param == string("GROUND_REPULSION") )
		{
			fin >> m_groundRepulsion;
		}
		else if ( param == string("GROUND_FRICTION") )
		{
			fin >> m_groundFriction;
		}
		else if ( param == string("GROUND_ABSORPTION") )
		{
			fin >> m_groundAbsorption;
		}
		else if ( param == string("SPATIAL_WEIGHT") )
		{
			fin >> m_spatialWt;
		}
		else if ( param == string("SPATIAL_SIGMA") )
		{
			fin >> m_spatialSigma;
		}
		else if ( param == string("SPATIAL_NEIGHBOR_DISTANCE") )
		{
			fin >> m_spatialNeighDist;
		}
		else if ( param == string("TEMPORAL_WEIGHT") )
		{
			fin >> m_temporalWt;
		}
		else if ( param == string("TEMPORAL_SIGMA") )
		{
			fin >> m_temporalSigma;
		}
		else if ( param == string("COMPATIBILITY_THRESH") )
		{
			fin >> m_compatibilityThresh;
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
	return true;
}

bool CSpaghettiSynConfig::DumpToSpaghettiSynConfig()
{
	ostringstream oss;
	oss << m_outputPrefix << "SpaghettiSynConfig.txt";
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

void CSpaghettiSynConfig::ResetConfig()
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
	DumpToSpaghettiSynConfig();
}

void CSpaghettiSynConfig::DumpParameters(FILE* file)
{
	CSynConfigBase::DumpParameters(file);
	fprintf(file, "%s\t%d\n", "FLAG_TEMPORALLY_TOROIDAL", m_flagTemporallyToroidal);
	fprintf(file, "%s\t%d\n", "NUMBER_OF_MASS", m_numOfMasses);
	fprintf(file, "%s\t%d\n", "SPATIAL_NEIGHBOR_SIZE", m_spatialNeighSize);
	fprintf(file, "%s\t%d\n", "TEMPORAL_WINDOW_SIZE", m_temporalWindowSize);
	fprintf(file, "%s\t%d\n", "ITERATION_NUMBER", m_iterationNum);
	fprintf(file, "%s\t%d\n", "INPUT_LOAD_INTERVAL", m_inputLoadInterval);
	fprintf(file, "%s\t%d\n", "INPUT_LOAD_START", m_inputLoadStart);
	fprintf(file, "%s\t%d\n", "INPUT_LOAD_END", m_inputLoadEnd);
	fprintf(file, "%s\t%d\n", "COARSE_LOAD_INTERVAL", m_coarseLoadInterval);
	fprintf(file, "%s\t%d\n", "COARSE_LOAD_START", m_coarseLoadStart);
	fprintf(file, "%s\t%d\n", "COARSE_LOAD_END", m_coarseLoadEnd);
	fprintf(file, "%s\t%d\n", "INTERLEAVED_STEPS", m_interleavedSteps);
	fprintf(file, "%s\t%f\n", "SPRING_LENGTH", m_springLength);
	fprintf(file, "%s\t%f\n", "SPRING_RADIUS", m_springRadius);
	fprintf(file, "%s\t%f\n", "SPRING_STIFFNESS", m_springStiffness);
	fprintf(file, "%s\t%f\n", "SPRING_FRICTION", m_springFriction);
	fprintf(file, "%s\t%f\n", "SPRING_REPULSION", m_springRepulsion);
	fprintf(file, "%s\t%f\n", "GRAVITY", m_gravity);
	fprintf(file, "%s\t%f\n", "AIR_FRICTION", m_airFriction);
	fprintf(file, "%s\t%f\n", "GROUND_HEIGHT", m_groundHeight);
	fprintf(file, "%s\t%f\n", "GROUND_REPULSION", m_groundRepulsion);
	fprintf(file, "%s\t%f\n", "GROUND_FRICTION", m_groundFriction);
	fprintf(file, "%s\t%f\n", "GROUND_ABSORPTION", m_groundAbsorption);
	fprintf(file, "%s\t%f\n", "SPATIAL_WEIGHT", m_spatialWt);
	fprintf(file, "%s\t%f\n", "SPATIAL_SIGMA", m_spatialSigma);
	fprintf(file, "%s\t%f\n", "SPATIAL_NEIGHBOR_DISTANCE", m_spatialNeighDist);
	fprintf(file, "%s\t%f\n", "TEMPORAL_WEIGHT", m_temporalWt);
	fprintf(file, "%s\t%f\n", "TEMPORAL_SIGMA", m_temporalSigma);
	fprintf(file, "%s\t%f\n", "COMPATIBILITY_THRESH", m_compatibilityThresh);
	fprintf(file, "%s\t%f %f\n", "POISSON_2D_MIN", m_poisson2Dmin[0], m_poisson2Dmin[1]);
	fprintf(file, "%s\t%f %f\n", "POISSON_2D_MAX", m_poisson2Dmax[0], m_poisson2Dmax[1]);
}
