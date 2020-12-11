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
Vec2f CMassSpringSynConfig::m_poisson2Dmin = Vec2f(-0.1f, -0.1f);
Vec2f CMassSpringSynConfig::m_poisson2Dmax = Vec2f(0.1f, 0.1f);

CMassSpringSynConfig::CMassSpringSynConfig(const std::string& config_file_name)
{
    LoadFromMassSpringSynConfig(config_file_name);
    ResetConfig();
}

bool CMassSpringSynConfig::ReloadConfigFromFile(const std::string& config_file_name)
{
    bool flag = LoadFromMassSpringSynConfig(config_file_name);
    if (flag == false) return false;
    ResetConfig();
    return true;
}

bool CMassSpringSynConfig::LoadFromMassSpringSynConfig(const std::string& config_file_name)
{
    CSynConfigBase::LoadFromSynConfig(config_file_name.c_str());
    ifstream fin(config_file_name.c_str());
    if (fin.fail() == true)
    {
        cout << "Failed to load parameters from config file " << config_file_name << "!\n";
        return false;
    }
    string param;
    while (fin >> param)
    {
        if (param == string("FLAG_TEMPORALLY_TOROIDAL"))
        {
            fin >> m_flagTemporallyToroidal;
        }
        else if (param == string("NUM_OF_INPUT_FRAMES"))
        {
            fin >> m_numOfInputFrames;
        }
        else if (param == string("NUMBER_OF_MASS"))
        {
            fin >> m_numOfMasses;
        }
        else if (param == string("NUMBER_OF_STRAND"))
        {
            fin >> m_numOfStrands;
        }
        else if (param == string("NEIGHBOR_SIZE"))
        {
            fin >> m_neighSize;
        }
        else if (param == string("INPUT_LOAD_INTERVAL"))
        {
            fin >> m_inputLoadInterval;
        }
        else if (param == string("INPUT_LOAD_START"))
        {
            fin >> m_inputLoadStart;
        }
        else if (param == string("NUMBER_OF_INPUT_STRANDS"))
        {
            fin >> m_numOfInputStrands;
        }
        else if (param == string("SYNTHESIS_WINDOW_SIZE"))
        {
            fin >> m_synthesisWindowSize;
        }
        else if (param == string("EXEMPLAR_NAME"))
        {
            fin >> m_exemplarName;
        }
        else if (param == string("INPUT_WEIGHT"))
        {
            Flt inputWt;
            fin >> inputWt;
            m_vecInputWt.push_back(inputWt);
        }
        else if (param == string("SPRING_LENGTH"))
        {
            fin >> m_springLength;
        }
        else if (param == string("WEIGHT_INDEX_DIFFERENCE"))
        {
            fin >> m_wtIdxDiff;
        }
        else if (param == string("WEIGHT_PR_DIFFERENCE"))
        {
            fin >> m_wtPrDiff;
        }
        else if (param == string("WEIGHT_NEIGHBOR_PR_DIFFERENCE"))
        {
            fin >> m_wtNeighPrDiff;
        }
        else if (param == string("SHAPE_WEIGHT"))
        {
            fin >> m_shapeWt;
        }
        else if (param == string("SMOOTH_WEIGHT"))
        {
            fin >> m_smoothWt;
        }
        else if (param == string("ALIGN_WEIGHT"))
        {
            fin >> m_alignWt;
        }
        else if (param == string("TEMPORAL_WEIGHT"))
        {
            fin >> m_temporalWt;
        }
        else if (param == string("STRAND_RADIUS"))
        {
            fin >> m_strandRadius;
        }
        else if (param == string("STRAND_REPULSION_CONST"))
        {
            fin >> m_strandRepulsionConstant;
        }
        else if (param == string("NEIGHBOR_DISTANCE"))
        {
            fin >> m_neighDist;
        }
        else if (param == string("GAUSSIAN_SIGMA"))
        {
            fin >> m_gaussianSigma;
        }
        else if (param == string("TEMPORAL_SIGMA"))
        {
            fin >> m_temporalSigma;
        }
        else if (param == string("POISSON_2D_MIN"))
        {
            fin >> m_poisson2Dmin[0] >> m_poisson2Dmin[1];
        }
        else if (param == string("POISSON_2D_MAX"))
        {
            fin >> m_poisson2Dmax[0] >> m_poisson2Dmax[1];
        }
    }
    if ((m_vecInputPrefix.size() > 0))
    {
        if (m_vecInputWt.size() != m_vecInputPrefix.size())
        {
            m_vecInputWt.resize(m_vecInputPrefix.size(), 1.0f);
        }
        Flt inputWtSum = 0.0f;
        for (int i = 0; i < int(m_vecInputWt.size()); i++)
        {
            inputWtSum += m_vecInputWt[i];
        }
        m_inputWtSumInv = 1.0f / inputWtSum;
    }
    return true;
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
