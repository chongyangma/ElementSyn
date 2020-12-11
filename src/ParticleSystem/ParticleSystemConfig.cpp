#include "ParticleSystemConfig.h"

bool CParticleSystemConfig::m_flagHybridSolver = true;
bool CParticleSystemConfig::m_flagSmoothSynthesis = true;
int CParticleSystemConfig::m_initMethod = 1;
int CParticleSystemConfig::m_minimumPatchSize = 1;
int CParticleSystemConfig::m_iterationNum = 1;
int CParticleSystemConfig::m_temporalWindowSize = 1;
string CParticleSystemConfig::m_trajectoriesFileName = "trajectories.txt";
vector<string> CParticleSystemConfig::m_vecInputFileName;
Flt CParticleSystemConfig::m_shapeWt = 0.0f;
Flt CParticleSystemConfig::m_neighWt = 0.0f;
Flt CParticleSystemConfig::m_spatialWt = 0.0f;
Flt CParticleSystemConfig::m_boundaryWt = 0.0f;
Flt CParticleSystemConfig::m_temporalWt = 0.0f;
Flt CParticleSystemConfig::m_neighDist = 0.0f;
Flt CParticleSystemConfig::m_repulsionDist = 0.0f;
Flt CParticleSystemConfig::m_repulsionWt = 0.0f;
Flt CParticleSystemConfig::m_compatibilityThresh = 5.0f;
Flt CParticleSystemConfig::m_boundaryThresh = 0.0f;
Flt CParticleSystemConfig::m_gaussianSigma = 0.0f;
Vec3f CParticleSystemConfig::m_inputTrans = Vec3f(0.0f, 0.0f, 0.0f);
Vec3f CParticleSystemConfig::m_outputTrans = Vec3f(0.0f, 0.0f, 0.0f);
Vec3f CParticleSystemConfig::m_cubicPosMin = Vec3f(-5.f, -5.f, -5.f);
Vec3f CParticleSystemConfig::m_cubicPosMax = Vec3f(5.f, 5.f, 5.f);
Vec3f CParticleSystemConfig::m_patchSize = Vec3f(1.f, 1.f, 1.f);
Vec3f CParticleSystemConfig::m_particleAmbient = Vec3f(0.5f, 0.5f, 0.5f);
Vec3f CParticleSystemConfig::m_particleDiffuse = Vec3f(0.5f, 0.5f, 0.5f);
Vec3f CParticleSystemConfig::m_particleSpecular = Vec3f(0.1f, 0.1f, 0.1f);
vector<Vec3f> CParticleSystemConfig::m_vecDiffuseColor;

CParticleSystemConfig::CParticleSystemConfig(const std::string& fileName)
{
    LoadFromParticleSystemConfig(fileName);
    ResetConfig();
}

bool CParticleSystemConfig::ReloadConfigFromFile(const string& fileName)
{
    bool flag = LoadFromParticleSystemConfig(fileName);
    if (flag == false) return false;
    ResetConfig();
    return true;
}

bool CParticleSystemConfig::LoadFromParticleSystemConfig(const string& fileName)
{
    CSynConfigBase::LoadFromSynConfig(fileName);
    ifstream fin(fileName.c_str());
    if (fin.fail() == true)
    {
        cout << "Failed to load parameters from config file " << fileName << "!\n";
        return false;
    }
    m_vecInputFileName.clear();
    m_vecDiffuseColor.clear();
    string param;
    while (fin >> param)
    {
        if (param == string("FLAG_HYBRID_SOLVER"))
        {
            fin >> m_flagHybridSolver;
        }
        else if (param == string("FLAG_SMOOTH_SYNTHESIS"))
        {
            fin >> m_flagSmoothSynthesis;
        }
        else if (param == string("INIT_METHOD"))
        {
            fin >> m_initMethod;
        }
        else if (param == string("MINIMUM_PATCH_SIZE"))
        {
            fin >> m_minimumPatchSize;
        }
        else if (param == string("ITERATION_NUMBER"))
        {
            fin >> m_iterationNum;
        }
        else if (param == string("TEMPORAL_WINDOW_SIZE"))
        {
            fin >> m_temporalWindowSize;
        }
        else if (param == string("TRAJECTORIES_FILE_NAME"))
        {
            fin >> m_trajectoriesFileName;
        }
        else if (param == string("INPUT_FILE_NAME"))
        {
            string inputFileName;
            fin >> inputFileName;
            m_vecInputFileName.push_back(inputFileName);
        }
        else if (param == string("SHAPE_WEIGHT"))
        {
            fin >> m_shapeWt;
        }
        else if (param == string("NEIGH_WEIGHT"))
        {
            fin >> m_neighWt;
        }
        else if (param == string("SPATIAL_WEIGHT"))
        {
            fin >> m_spatialWt;
        }
        else if (param == string("BOUNDARY_WEIGHT"))
        {
            fin >> m_boundaryWt;
        }
        else if (param == string("TEMPORAL_WEIGHT"))
        {
            fin >> m_temporalWt;
        }
        else if (param == string("NEIGHBOR_DISTANCE"))
        {
            fin >> m_neighDist;
        }
        else if (param == string("REPULSION_DISTANCE"))
        {
            fin >> m_repulsionDist;
        }
        else if (param == string("REPULSION_WEIGHT"))
        {
            fin >> m_repulsionWt;
        }
        else if (param == string("COMPATIBILITY_THRESH"))
        {
            fin >> m_compatibilityThresh;
        }
        else if (param == string("BOUNDARY_THRESH"))
        {
            fin >> m_boundaryThresh;
        }
        else if (param == string("GAUSSIAN_SIGMA"))
        {
            fin >> m_gaussianSigma;
        }
        else if (param == string("INPUT_TRANSLATION"))
        {
            Flt px, py, pz;
            fin >> px >> py >> pz;
            m_inputTrans = Vec3f(px, py, pz);
        }
        else if (param == string("OUTPUT_TRANSLATION"))
        {
            Flt px, py, pz;
            fin >> px >> py >> pz;
            m_outputTrans = Vec3f(px, py, pz);
        }
        else if (param == string("CUBIC_POS_MIN"))
        {
            fin >> m_cubicPosMin[0] >> m_cubicPosMin[1] >> m_cubicPosMin[2];
        }
        else if (param == string("CUBIC_POS_MAX"))
        {
            fin >> m_cubicPosMax[0] >> m_cubicPosMax[1] >> m_cubicPosMax[2];
        }
        else if (param == string("PATCH_SIZE"))
        {
            fin >> m_patchSize[0] >> m_patchSize[1] >> m_patchSize[2];
        }
        else if (param == string("PARTICLE_AMBIENT"))
        {
            fin >> m_particleAmbient[0] >> m_particleAmbient[1] >> m_particleAmbient[2];
        }
        else if (param == string("PARTICLE_DIFFUSE"))
        {
            fin >> m_particleDiffuse[0] >> m_particleDiffuse[1] >> m_particleDiffuse[2];
        }
        else if (param == string("PARTICLE_SPECULAR"))
        {
            fin >> m_particleSpecular[0] >> m_particleSpecular[1] >> m_particleSpecular[2];
        }
        else if (param == string("DIFFUSE_COLOR"))
        {
            Vec3f diffuseColor;
            fin >> diffuseColor[0] >> diffuseColor[1] >> diffuseColor[2];
            m_vecDiffuseColor.push_back(diffuseColor);
        }
    }
    if (m_flagRandomness == true)
    {
        srand((unsigned int)time(0));
    }
    else
    {
        srand(0); // Reset rand seed!
    }

    return true;
}

void CParticleSystemConfig::DumpParameters(FILE* file)
{
    CSynConfigBase::DumpParameters(file);
    fprintf(file, "%s\t%d\n", "FLAG_HYBRID_SOLVER", m_flagHybridSolver);
    fprintf(file, "%s\t%d\n", "FLAG_SMOOTH_SYNTHESIS", m_flagSmoothSynthesis);
    fprintf(file, "%s\t%d\n", "INIT_METHOD", m_initMethod);
    fprintf(file, "%s\t%d\n", "MINIMUM_PATCH_SIZE", m_minimumPatchSize);
    fprintf(file, "%s\t%d\n", "ITERATION_NUMBER", m_iterationNum);
    fprintf(file, "%s\t%d\n", "TEMPORAL_WINDOW_SIZE", m_temporalWindowSize);
    DumpStringParam(file, "TRAJECTORIES_FILE_NAME", m_trajectoriesFileName);
    for (int i = 0; i < int(m_vecInputFileName.size()); i++)
    {
        DumpStringParam(file, "INPUT_FILE_NAME", m_vecInputFileName[i]);
    }
    fprintf(file, "%s\t%f\n", "SHAPE_WEIGHT", m_shapeWt);
    fprintf(file, "%s\t%f\n", "NEIGH_WEIGHT", m_neighWt);
    fprintf(file, "%s\t%f\n", "SPATIAL_WEIGHT", m_spatialWt);
    fprintf(file, "%s\t%f\n", "BOUNDARY_WEIGHT", m_boundaryWt);
    fprintf(file, "%s\t%f\n", "TEMPORAL_WEIGHT", m_temporalWt);
    fprintf(file, "%s\t%f\n", "NEIGHBOR_DISTANCE", m_neighDist);
    fprintf(file, "%s\t%f\n", "REPULSION_DISTANCE", m_repulsionDist);
    fprintf(file, "%s\t%f\n", "REPULSION_WEIGHT", m_repulsionWt);
    fprintf(file, "%s\t%f\n", "COMPATIBILITY_THRESH", m_compatibilityThresh);
    fprintf(file, "%s\t%f\n", "BOUNDARY_THRESH", m_boundaryThresh);
    fprintf(file, "%s\t%f\n", "GAUSSIAN_SIGMA", m_gaussianSigma);
    fprintf(file, "%s\t%f %f %f\n", "INPUT_TRANSLATION", m_inputTrans[0], m_inputTrans[1], m_inputTrans[2]);
    fprintf(file, "%s\t%f %f %f\n", "OUTPUT_TRANSLATION", m_outputTrans[0], m_outputTrans[1], m_outputTrans[2]);
    fprintf(file, "%s\t%f %f %f\n", "CUBIC_POS_MIN", m_cubicPosMin[0], m_cubicPosMin[1], m_cubicPosMin[2]);
    fprintf(file, "%s\t%f %f %f\n", "CUBIC_POS_MAX", m_cubicPosMax[0], m_cubicPosMax[1], m_cubicPosMax[2]);
    fprintf(file, "%s\t%f %f %f\n", "PATCH_SIZE", m_patchSize[0], m_patchSize[1], m_patchSize[2]);
    fprintf(file, "%s\t%f %f %f\n", "PARTICLE_AMBIENT", m_particleAmbient[0], m_particleAmbient[1], m_particleAmbient[2]);
    fprintf(file, "%s\t%f %f %f\n", "PARTICLE_DIFFUSE", m_particleDiffuse[0], m_particleDiffuse[1], m_particleDiffuse[2]);
    fprintf(file, "%s\t%f %f %f\n", "PARTICLE_SPECULAR", m_particleSpecular[0], m_particleSpecular[1], m_particleSpecular[2]);
    for (int i = 0; i < int(m_vecDiffuseColor.size()); i++)
    {
        Vec3f diffuseColor = m_vecDiffuseColor[i];
        fprintf(file, "%s\t%f %f %f\n", "DIFFUSE_COLOR", diffuseColor[0], diffuseColor[1], diffuseColor[2]);
    }
}
