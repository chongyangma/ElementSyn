
#include "SynConfigBase.h"

bool CSynConfigBase::m_flagRandomness = false;
int CSynConfigBase::m_stepCountMax = 2000;
string CSynConfigBase::m_inputPrefix;
string CSynConfigBase::m_coarsePrefix;
string CSynConfigBase::m_outputPrefix;
Flt CSynConfigBase::m_timeStep = 0.0025f;

CSynConfigBase::CSynConfigBase()
{
	//LoadFromSynConfig();
	//ResetConfig();
}

bool CSynConfigBase::LoadFromSynConfig(const string& fileName)
{
	ifstream fin(fileName.c_str());
	if (fin.fail() == true)
	{
		cout << "Failed to load base class parameters from config file " << fileName << "!\n";
		return false;
	}
	string param;
	while (fin >> param)
	{
		if (param == string("FLAG_RANDOMNESS"))
		{
			fin >> m_flagRandomness;
		}
		else if (param == string("STEP_COUNT_MAX"))
		{
			fin >> m_stepCountMax;
		}
		else if (param == string("INPUT_PREFIX"))
		{
			fin >> m_inputPrefix;
		}
		else if (param == string("COARSE_PREFIX"))
		{
			fin >> m_coarsePrefix;
		}
		else if (param == string("OUTPUT_PREFIX"))
		{
			fin >> m_outputPrefix;
		}
		else if (param == string("TIME_STEP"))
		{
			fin >> m_timeStep;
		}
	}

	return true;
}

bool CSynConfigBase::DumpToSynConfig()
{
	ostringstream oss;
	oss << m_outputPrefix << "SynConfig.txt";
	string fileName = oss.str();
	FILE* file;
	file = fopen(fileName.c_str(), "w");
	if (!file)
	{
		cout << "Failed to dump parameters into config file " << fileName << "!\n";
		return false;
	}
	DumpTimeAndDate(file);
	DumpParameters(file);
	fclose(file);

	return true;
}

void CSynConfigBase::ResetConfig()
{
	if (m_flagRandomness == true)
	{
		srand((unsigned int)time(0));
	}
	if (m_outputPrefix.empty() == true)
	{
		return;
	}
	UpdateOutputPrefix();
	DumpToSynConfig();
}

void CSynConfigBase::DumpTimeAndDate(FILE* file)
{
	time_t myTime = time(NULL);
	tm* ptrTime = localtime(&myTime);
	fprintf(file, "%02d:%02d:%02d ", ptrTime->tm_hour, ptrTime->tm_min, ptrTime->tm_sec);
	fprintf(file, "%02d/%02d/%04d\n\n", ptrTime->tm_mon + 1, ptrTime->tm_mday, ptrTime->tm_year + 1900);
}

void CSynConfigBase::DumpParameters(FILE* file)
{
	fprintf(file, "%s\t%d\n", "FLAG_RANDOMNESS", m_flagRandomness);
	fprintf(file, "%s\t%d\n", "STEP_COUNT_MAX", m_stepCountMax);
	DumpStringParam(file, "INPUT_PREFIX", m_inputPrefix);
	DumpStringParam(file, "COARSE_PREFIX", m_coarsePrefix);
	DumpStringParam(file, "OUTPUT_PREFIX", m_outputPrefix);
	fprintf(file, "%s\t%f\n", "TIME_STEP", m_timeStep);
}

void CSynConfigBase::UpdateOutputPrefix()
{
	if (m_outputPrefix[m_outputPrefix.size() - 1] != '\\')
	{
#ifdef WIN32
		m_outputPrefix = m_outputPrefix + string("\\");
#else
		m_outputPrefix = m_outputPrefix + string("/");
#endif
	}
	ostringstream oss;
#ifdef WIN32
	BOOL flag = CreateDirectoryA(m_outputPrefix.c_str(), NULL);
	const int numLength = 2;
	while (flag == false && m_outputPrefix.size() >= 2)
	{
		string subStr = m_outputPrefix.substr(m_outputPrefix.length() - numLength - 1, 2);
		int num = atoi(subStr.c_str()) + 1;
		char numChar[MAX_PATH];
		sprintf(numChar, "%02d", num);
		ostringstream oss;
		oss << m_outputPrefix.substr(0, m_outputPrefix.length() - numLength - 1) << numChar << "\\";
		m_outputPrefix = oss.str();
		flag = CreateDirectoryA(m_outputPrefix.c_str(), NULL);
	}
	cout << "Dump results to the directory: " << m_outputPrefix.c_str() << endl;
	ostringstream oss2;
	oss2 << CSynConfigBase::m_outputPrefix << "Snapshot";
	CreateDirectory(oss2.str().c_str(), NULL);
	ostringstream oss3;
	oss3 << CSynConfigBase::m_outputPrefix << "Dumped";
	CreateDirectory(oss3.str().c_str(), NULL);
#else
	string outputFolder = m_outputPrefix.substr(0, m_outputPrefix.size() - 1);
	int flag = mkdir(outputFolder.c_str(), 0777);
	const int numLength = 2;
	while (flag == -1 && m_outputPrefix.size() >= 2)
	{
		string subStr = m_outputPrefix.substr(m_outputPrefix.length() - numLength - 1, 2);
		int num = atoi(subStr.c_str()) + 1;
		char numChar[MAX_PATH];
		sprintf(numChar, "%02d", num);
		ostringstream oss;
		oss << m_outputPrefix.substr(0, m_outputPrefix.length() - numLength - 1) << numChar << "/";
		m_outputPrefix = oss.str();
		outputFolder = m_outputPrefix.substr(0, m_outputPrefix.size() - 1);
		flag = mkdir(outputFolder.c_str(), 0777);
	}
	cout << "Dump results to the directory: " << m_outputPrefix.c_str() << endl;
	ostringstream oss2;
	oss2 << CSynConfigBase::m_outputPrefix << "Snapshot";
	mkdir(oss2.str().c_str(), 0777);
	ostringstream oss3;
	oss3 << CSynConfigBase::m_outputPrefix << "Dumped";
	mkdir(oss3.str().c_str(), 0777);
#endif
}

void CSynConfigBase::DumpStringParam(FILE* file, const char* param, const string str)
{
	if (str.empty() == true)
	{
		//fprintf(file, "%s\t%s\n", param, "NULL");
	}
	else
	{
		fprintf(file, "%s\t%s\n", param, str.c_str());
	}
}
