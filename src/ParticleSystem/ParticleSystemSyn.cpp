#include "ParticleSystemSyn.h"
//#define MORPHABLE_PARTICLE_SYSTEM
//#define CONSTRAINED_SYNTHESIS
//#define FLOW_SYNTHESIS // 01/07/2013

CParticleSystemSyn::CParticleSystemSyn(const std::string& config_file_name)
{
    m_stepCount = 0;
    m_serialCount = 1;
    m_ptrSynConfig = new CParticleSystemConfig(config_file_name);
#ifdef FLOW_SYNTHESIS
    LoadInputDataNew();
#else
    LoadInputData();
#endif
    InitializeOutput();
}

CParticleSystemSyn::~CParticleSystemSyn()
{
    DELETE_OBJECT(m_ptrSynConfig);
}

void CParticleSystemSyn::SetInputNeighborhoods()
{
    for (int i = 0; i < int(m_inputSequence.size()); i++)
    {
        m_inputSequence[i].SetNeighboringSamples(CParticleSystemConfig::m_neighDist * 2.0f);
    }
    for (int i = 0; i < int(m_vecInputExemplar.size()); i++)
    {
        m_vecInputExemplar[i].SetNeighboringSamples(CParticleSystemConfig::m_neighDist * 2.0f);
    }
}

void CParticleSystemSyn::LoadInputData()
{
    m_inputSequence.clear();
    CParticleSystem inputGroup; // = LoadInputExemplarFromTXT(CParticleSystemConfig::m_trajectoriesFileName);
    bool flag = inputGroup.LoadParticleSystem(CParticleSystemConfig::m_trajectoriesFileName);
    if (flag == false)
    {
        exit(2);
    }
    char fileName[MAX_PATH];
    sprintf(fileName, "%sinput_exemplar.csv", CSynConfigBase::m_outputPrefix.c_str());
    inputGroup.SaveParticleSystemAsCSV(fileName);
    m_inputSequence.push_back(inputGroup);

    m_vecInputExemplar.clear();
    for (int i = 0; i < int(CParticleSystemConfig::m_vecInputFileName.size()); i++)
    {
        CParticleSystem inputExemplar = LoadInputExemplarFromTXT(CParticleSystemConfig::m_vecInputFileName[i]);
        m_vecInputExemplar.push_back(inputExemplar);
    }

    SetInputNeighborhoods();
    //LoadBoundaryConstraints();
}

void CParticleSystemSyn::LoadInputDataNew()
{
    m_inputSequence.clear();
    m_vecInputExemplar.clear();
    vector<CParticleSystem> vecInputExemplar;
    for (int n = 1; n <= 6; n++)
    {
        char fileName[MAX_PATH];
        sprintf(fileName, "test%02d.txt", n);
        ifstream fin(fileName);
        CParticleSystem particles;
        for (int i = 0; i < 60; i++)
        {
            Flt px, py;
            fin >> px >> py;
            Vec3f pi = Vec3f(px, py, 0.f);
            pi += Vec3f(-0.5f, -0.5f, 0.0f);
            pi *= 1.2f;
            vector<Vec3f> vertices;
            vertices.push_back(pi);
            CParticleData particle;
            particle.SetPos(pi);
            particle.SetVecSamplePos(vertices);
            particles.AddParticle(particle);
        }
        vecInputExemplar.push_back(particles);
        //m_inputSequence.push_back(particles);
    }
    for (int n = 0; n < int(vecInputExemplar.size()); n++)
    {
        int n1 = n;
        int n2 = (n + 1) % int(vecInputExemplar.size());
        CParticleSystem& particles1 = vecInputExemplar[n1];
        CParticleSystem& particles2 = vecInputExemplar[n2];
        vector<CParticleSystem> ps = InterpolateParticleSystem(particles1, particles2, 30);
        for (int i = 0; i < int(ps.size()); i++)
        {
            m_vecInputExemplar.push_back(ps[i]);
        }
    }
    SetSequencesVelocity(m_vecInputExemplar, false);
    //m_vecInputExemplar = m_inputSequence;
    m_inputSequence.push_back(m_vecInputExemplar[0]);
    SetInputNeighborhoods();
}

void CParticleSystemSyn::UpdateOutput()
{
#ifdef CONSTRAINED_SYNTHESIS
    int idx = min(m_stepCount, int(m_outputSequence.size()) - 1);
    m_outputGroup = m_outputSequence[idx];
    CParticleSystemConfig::m_particleDiffuse = m_diffuseColorSequence[idx];
#else
    #ifdef FLOW_SYNTHESIS
    AdvectOutputNew();
    m_inputSequence.clear();
    for (int i = 0; i < CParticleSystemConfig::m_temporalWindowSize; i++)
    {
        int inputIdx = (m_stepCount + i) % int(m_vecInputExemplar.size());
        m_inputSequence.push_back(m_vecInputExemplar[inputIdx]);
    }
    #else
    AdvectOutput();
    #endif
    for (int i = 0; i < CParticleSystemConfig::m_iterationNum; i++)
    {
        UpdateOutputGroupViaOptimization(m_outputGroup);
    }
#endif
    char fileName[MAX_PATH];
#ifdef WIN32
    sprintf(fileName, "%sDumped\\dumped_%04d.csv", CParticleSystemConfig::m_outputPrefix.c_str(), GetStepCount());
#else
    sprintf(fileName, "%sDumped/dumped_%04d.csv", CParticleSystemConfig::m_outputPrefix.c_str(), GetStepCount());
#endif
    m_outputGroup.SaveParticleSystem(fileName);
    m_stepCount++;
}

void CParticleSystemSyn::RenderInput()
{
    if (m_inputSequence.empty())
    {
        return;
    }
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    CParticleSystem& inputExemplar = m_inputSequence[0];
    int numOfParticles = inputExemplar.GetNumOfSoftBodies();
    glPushMatrix();
    glTranslatef(-2.25f, -1.25f, 0.0f); // Render the input at the lower left corner
    for (int i = 0; i < numOfParticles; i++)
    {
        Vec3f pi = inputExemplar.GetParticleData(i).GetPos();
        inputExemplar.GetParticleData(i).RenderSoftBody(pi);
    }
    glPopMatrix();
    glDisable(GL_LIGHTING);
}

void CParticleSystemSyn::RenderOutput()
{
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    // Render output...
    int numOfSoftBodies = m_outputGroup.GetNumOfSoftBodies();
    glPushMatrix();
    //glTranslatef(0.5f, 0.0f, 0.0f);
    glPointSize(3.0f);
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        Vec3f pi = m_outputGroup.GetParticleData(i).GetVecSamplePos()[0];
        m_outputGroup.GetParticleData(i).RenderSoftBody(pi);
    }
    glPopMatrix();
    glDisable(GL_LIGHTING);
}

void CParticleSystemSyn::RestartSynthesis()
{
    char fileName[MAX_PATH];
    sprintf(fileName, "ParticleSystemSyn_config%02d.txt", m_serialCount++);
    ResetOutput(fileName);
    m_stepCount = 0;
}

void CParticleSystemSyn::ResetOutput(const std::string& config_file_name)
{
    m_stepCount = 0;
    bool flag = m_ptrSynConfig->ReloadConfigFromFile(config_file_name);
    if (flag == false)
    {
        cout << "Have used all the configuration files!\n";
        exit(0);
    }
#ifdef FLOW_SYNTHESIS
    LoadInputDataNew();
#else
    LoadInputData();
#endif
    InitializeOutput();
}

void CParticleSystemSyn::InitializeOutput()
{
#ifndef CONSTRAINED_SYNTHESIS
    m_outputGroup.ClearParticleSystem();
    InitializeOutputTrajectoriesViaPatchCopy();
    for (int i = 0; i < 10; i++)
    {
        UpdateOutputGroupViaOptimization(m_outputGroup);
    }
    //m_outputGroup.GetParticleData(50).SetFlagFixed(true);
#else
    CBoundaryConstraint* ptrBoundaryConstraint0 = &(m_vecBoundaryConstraint[0]);
    CBoundaryConstraint* ptrBoundaryConstraint1 = NULL;
    m_outputSequence.clear();
    InitializeOutputTrajectoriesViaPatchCopy3(-1, ptrBoundaryConstraint0);
    for (int i = 0; i < 10; i++)
    {
        UpdateOutputGroupViaOptimization3(m_outputGroup, ptrBoundaryConstraint0);
    }
    int numOfParticles = m_outputGroup.GetNumOfSoftBodies();
    CParticleSystem outputGroup0 = m_outputGroup;
    for (int n = 1; n < int(m_vecBoundaryConstraint.size()); n++)
    {
        ptrBoundaryConstraint1 = &(m_vecBoundaryConstraint[n]);
        InitializeOutputTrajectoriesViaPatchCopy3(numOfParticles, ptrBoundaryConstraint1);
        for (int i = 0; i < 10; i++)
        {
            UpdateOutputGroupViaOptimization3(m_outputGroup, ptrBoundaryConstraint1);
        }
        int interpLength = ptrBoundaryConstraint1->GetFrameIdx() - ptrBoundaryConstraint0->GetFrameIdx();
        vector<CParticleSystem> outputSequence = InterpolateParticleSystem(outputGroup0, m_outputGroup, interpLength);
        for (int i = 0; i < int(outputSequence.size()); i++)
        {
            m_outputSequence.push_back(outputSequence[i]);
        }
        ptrBoundaryConstraint0 = ptrBoundaryConstraint1;
        outputGroup0 = m_outputGroup;
    }
    m_diffuseColorSequence.resize(m_outputSequence.size());
    vector<Vec3f>& vecDiffuseColor = CParticleSystemConfig::m_vecDiffuseColor;
    int n0 = 0;
    int n1 = 1;
    ptrBoundaryConstraint0 = &(m_vecBoundaryConstraint[n0]);
    ptrBoundaryConstraint1 = &(m_vecBoundaryConstraint[n1]);
    for (int j = 0; j < int(m_outputSequence.size()); j++)
    {
        int frameIdx1 = ptrBoundaryConstraint0->GetFrameIdx();
        int frameIdx2 = ptrBoundaryConstraint1->GetFrameIdx();
        Flt frameWt1 = (frameIdx2 - j) / Flt(frameIdx2 - frameIdx1);
        frameWt1 = min(frameWt1, 1.0f);
        frameWt1 = max(frameWt1, 0.0f);
        Flt frameWt2 = 1.0f - frameWt1;
        m_diffuseColorSequence[j] = frameWt1 * vecDiffuseColor[n0] + frameWt2 * vecDiffuseColor[n1];
        //m_diffuseColorSequence[j] = vecDiffuseColor[0];
        if (j == ptrBoundaryConstraint1->GetFrameIdx())
        {
            n0 = n1;
            n1 = (n1 + 1) % int(m_vecBoundaryConstraint.size());
            ptrBoundaryConstraint0 = ptrBoundaryConstraint1;
            ptrBoundaryConstraint1 = &(m_vecBoundaryConstraint[n1]);
        }
    }
    for (int i = 0; i < CParticleSystemConfig::m_iterationNum; i++)
    {
        int n0 = 0;
        int n1 = 1;
        //for ( int j=int(m_outputSequence.size())-1; j>=0; j-- )
        for (int j = 0; j < int(m_outputSequence.size()); j++)
        {
            ptrBoundaryConstraint0 = &(m_vecBoundaryConstraint[n0]);
            ptrBoundaryConstraint1 = &(m_vecBoundaryConstraint[n1]);
            UpdateOutputFrameViaBlending(j, ptrBoundaryConstraint0, ptrBoundaryConstraint1);
            if (j == ptrBoundaryConstraint1->GetFrameIdx())
            {
                n0 = n1;
                n1 = (n1 + 1) % int(m_vecBoundaryConstraint.size());
                ptrBoundaryConstraint0 = ptrBoundaryConstraint1;
                ptrBoundaryConstraint1 = &(m_vecBoundaryConstraint[n1]);
            }
        }
    }
    CParticleSystemConfig::m_stepCountMax = int(m_outputSequence.size()) + 30;
#endif
}

void CParticleSystemSyn::InitializeOutputTrajectoriesViaPatchCopy(int numOfParticles /* = -1 */, int boundaryType /* = 0 */)
{
    if (boundaryType == 1)
    {
        InitializeOutputTrajectoriesViaPatchCopy2(numOfParticles);
        return;
    }
    m_outputGroup.ClearParticleSystem();
    vector<CubicPatch> patches = GetCubicPatches(m_inputSequence);
    const int numOfPatches = int(patches.size());
    Vec3f patchSize = CParticleSystemConfig::m_patchSize;
    Vec3f cubicMin = CParticleSystemConfig::m_cubicPosMin;
    Vec3f cubicMax = CParticleSystemConfig::m_cubicPosMax;
    Vec3i patchNum;
    for (int i = 0; i < 3; i++)
    {
        patchNum[i] = floor((cubicMax[i] - cubicMin[i]) / patchSize[i] + 0.5f);
    }
    patchNum[2] = 1; // For 2D synthesis
    if (patchNum[2] == 1) cubicMin[2] = 0.0f;
    Vec3f patchCen;
    const int numOfInputFrames = int(m_inputSequence.size());
    for (int i = 0; i < patchNum[0]; i++)
    {
        patchCen[0] = cubicMin[0] + (i + 0.5f) * patchSize[0];
        for (int j = 0; j < patchNum[1]; j++)
        {
            patchCen[1] = cubicMin[1] + (j + 0.5f) * patchSize[1];
            for (int k = 0; k < patchNum[2]; k++)
            {
                patchCen[2] = cubicMin[2] + (k + 0.5f) * patchSize[2];
                int ri = int(rand() / Flt(RAND_MAX) * numOfPatches);
                ri = ri % numOfPatches;
                CubicPatch& patch = patches[ri];
                vector<int> indices = patch.m_patchIndices;
                int patchFrameIdx = patch.m_patchFrameIdx;
                Vec3f patchCenSrc = patch.m_patchCen;
                for (int n = 0; n < int(indices.size()); n++)
                {
                    int idxSrc = indices[n];
                    CParticleData softBodyData = m_inputSequence[patchFrameIdx].GetParticleData(idxSrc);
                    Vec3f center = softBodyData.GetPos() - patchCenSrc + patchCen;
                    softBodyData.SetPos(center);
                    softBodyData.TranslateSoftBody(patchCen - patchCenSrc);
                    if (numOfParticles > 0 && m_outputGroup.GetNumOfSoftBodies() >= numOfParticles)
                    {
                        continue;
                    }
                    m_outputGroup.AddParticle(softBodyData);
                } // End-For-n
            }     // End-For-k
        }         // End-For-j
    }             // End-For-i
    if (m_outputGroup.GetNumOfSoftBodies() < numOfParticles)
    {
        cout << "Try to initialize again...\n";
        InitializeOutputTrajectoriesViaPatchCopy(numOfParticles, boundaryType);
    }
    cout << "Have initialized " << m_outputGroup.GetNumOfSoftBodies() << " soft bodies via patch copy!\n";
}

void CParticleSystemSyn::InitializeOutputTrajectoriesViaPatchCopy2(int numOfParticles /* = -1 */)
{
    m_outputGroup.ClearParticleSystem();
    vector<CubicPatch> patches = GetCubicPatches(m_inputSequence);
    const int numOfPatches = int(patches.size());
    Vec3f patchSize = CParticleSystemConfig::m_patchSize;
    Vec3f cubicMin = CParticleSystemConfig::m_cubicPosMin;
    Vec3f cubicMax = CParticleSystemConfig::m_cubicPosMax;
    Vec3f circleCenter = (cubicMin + cubicMax) * 0.5f;
    Flt area = (cubicMax[0] - cubicMin[0]) * (cubicMax[1] - cubicMin[1]);
    Flt r = sqrt(area / M_PI);
    Vec3i patchNum;
    for (int i = 0; i < 3; i++)
    {
        cubicMin[i] = circleCenter[i] - r;
        cubicMax[i] = circleCenter[i] + r;
        patchNum[i] = floor((cubicMax[i] - cubicMin[i]) / patchSize[i] + 0.5f);
    }
    patchNum[2] = 1; // For 2D synthesis
    if (patchNum[2] == 1) cubicMin[2] = 0.0f;
    Vec3f patchCen;
    const int numOfInputFrames = int(m_inputSequence.size());
    vector<Vec3f> vecPatchCen;
    for (int i = 0; i < patchNum[0]; i++)
    {
        patchCen[0] = cubicMin[0] + (i + 0.5f) * patchSize[0];
        for (int j = 0; j < patchNum[1]; j++)
        {
            patchCen[1] = cubicMin[1] + (j + 0.5f) * patchSize[1];
            for (int k = 0; k < patchNum[2]; k++)
            {
                patchCen[2] = cubicMin[2] + (k + 0.5f) * patchSize[2];
                vecPatchCen.push_back(patchCen);
            }
        }
    }
    for (int i = 0; i < int(vecPatchCen.size()); i++)
    {
        patchCen = vecPatchCen[i];
        int ri = int(rand() / Flt(RAND_MAX) * numOfPatches);
        ri = ri % numOfPatches;
        CubicPatch& patch = patches[ri];
        vector<int> indices = patch.m_patchIndices;
        int patchFrameIdx = patch.m_patchFrameIdx;
        Vec3f patchCenSrc = patch.m_patchCen;
        for (int n = 0; n < int(indices.size()); n++)
        {
            int idxSrc = indices[n];
            CParticleData softBodyData = m_inputSequence[patchFrameIdx].GetParticleData(idxSrc);
            Vec3f center = softBodyData.GetPos() - patchCenSrc + patchCen;
            softBodyData.SetPos(center);
            softBodyData.TranslateSoftBody(patchCen - patchCenSrc);
            Vec3f pos = softBodyData.GetVecSamplePos()[0];
            if (dist2(pos, circleCenter) > r * r)
            {
                continue;
            }
            if (numOfParticles > 0 && m_outputGroup.GetNumOfSoftBodies() >= numOfParticles)
            {
                continue;
            }
            m_outputGroup.AddParticle(softBodyData);
        } // End-For-n
    }
    if (m_outputGroup.GetNumOfSoftBodies() < numOfParticles)
    {
        cout << "Try to initialize again...\n";
        InitializeOutputTrajectoriesViaPatchCopy2(numOfParticles);
    }
    cout << "Have initialized " << m_outputGroup.GetNumOfSoftBodies() << " soft bodies via patch copy within a circle!\n";
}

void CParticleSystemSyn::InitializeOutputTrajectoriesViaPatchCopy3(int numOfParticles /* = -1 */, CBoundaryConstraint* ptrBoundaryConstraint /* = NULL */)
{
    m_outputGroup.ClearParticleSystem();
    int inputIdx = ptrBoundaryConstraint->GetInputIdx();
    vector<CubicPatch> patches = GetCubicPatches(m_vecInputExemplar[inputIdx]); //GetCubicPatches(m_inputSequence); //
    const int numOfPatches = int(patches.size());
    Vec3f patchSize = CParticleSystemConfig::m_patchSize;
    Vec3f cubicMin = ptrBoundaryConstraint->GetPosMin();
    Vec3f cubicMax = ptrBoundaryConstraint->GetPosMax();
    Vec3i patchNum;
    for (int i = 0; i < 3; i++)
    {
        patchNum[i] = floor((cubicMax[i] - cubicMin[i]) / patchSize[i] + 0.5f);
    }
    patchNum[2] = 1; // For 2D synthesis
    if (patchNum[2] == 1) cubicMin[2] = 0.0f;
    Vec3f patchCen;
    //const int numOfInputFrames = int(m_inputSequence.size());
    vector<Vec3f> vecPatchCen;
    for (int i = 0; i < patchNum[0]; i++)
    {
        patchCen[0] = cubicMin[0] + (i + 0.5f) * patchSize[0];
        for (int j = 0; j < patchNum[1]; j++)
        {
            patchCen[1] = cubicMin[1] + (j + 0.5f) * patchSize[1];
            for (int k = 0; k < patchNum[2]; k++)
            {
                patchCen[2] = cubicMin[2] + (k + 0.5f) * patchSize[2];
                vecPatchCen.push_back(patchCen);
            }
        }
    }
    random_shuffle(vecPatchCen.begin(), vecPatchCen.end());
    for (int i = 0; i < int(vecPatchCen.size()); i++)
    {
        patchCen = vecPatchCen[i];
        int ri = int(rand() / Flt(RAND_MAX) * numOfPatches);
        ri = ri % numOfPatches;
        CubicPatch& patch = patches[ri];
        vector<int> indices = patch.m_patchIndices;
        int patchFrameIdx = patch.m_patchFrameIdx;
        Vec3f patchCenSrc = patch.m_patchCen;
        for (int n = 0; n < int(indices.size()); n++)
        {
            int idxSrc = indices[n];
            //CParticleData softBodyData = m_inputSequence[patchFrameIdx].GetParticleData(idxSrc);
            CParticleData softBodyData = m_vecInputExemplar[inputIdx].GetParticleData(idxSrc);
            Vec3f center = softBodyData.GetPos() - patchCenSrc + patchCen;
            softBodyData.SetPos(center);
            softBodyData.TranslateSoftBody(patchCen - patchCenSrc);
            Vec3f pos = softBodyData.GetVecSamplePos()[0];
            if (numOfParticles > 0 && m_outputGroup.GetNumOfSoftBodies() >= numOfParticles)
            {
                break;
            }
            if (ptrBoundaryConstraint->InsideBoundaryNew(pos) == false)
            {
                continue;
            }
            m_outputGroup.AddParticle(softBodyData);
        } // End-For-n
    }
    if (m_outputGroup.GetNumOfSoftBodies() < numOfParticles)
    {
        cout << "Try to initialize again...\n";
        InitializeOutputTrajectoriesViaPatchCopy3(numOfParticles, ptrBoundaryConstraint);
    }
    cout << "Have initialized " << m_outputGroup.GetNumOfSoftBodies() << " soft bodies via patch copy with a boundary constraint!\n";
}

void CParticleSystemSyn::AdvectOutput()
{
    Vec3f centerPos = Vec3f(0.0f, 0.0f, 0.0f);
    for (int i = 0; i < m_outputGroup.GetNumOfSoftBodies(); i++)
    {
        CParticleData& softBodyData = m_outputGroup.GetParticleData(i);
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        for (int j = 0; j < int(vecSamplePos.size()); j++)
        {
            Vec3f pj = vecSamplePos[j];
            Vec3f pd = pj - centerPos;
            if (mag2(pd) > 0.0f)
            {
                //normalize(pd);
                pd *= CParticleSystemConfig::m_timeStep;
                Vec3f pv = Vec3f(pd[1], -pd[0], 0.0f);
                vecSamplePos[j] = pj + pv;
            }
        }
    }
}

void CParticleSystemSyn::AdvectOutputNew()
{
    Vec3f centerPos = Vec3f(0.0f, 0.0f, 0.0f);
    for (int i = 0; i < m_outputGroup.GetNumOfSoftBodies(); i++)
    {
        CParticleData& softBodyData = m_outputGroup.GetParticleData(i);
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        Vec3f& vel = softBodyData.GetVel();
        for (int j = 0; j < int(vecSamplePos.size()); j++)
        {
            Vec3f pj = vecSamplePos[j];
            Vec3f pv = Vec3f(CParticleSystemConfig::m_timeStep, 0.0f, 0.0f);
            if (CParticleSystemConfig::m_iterationNum > 0)
            {
                pv += vel;
            }
            if (softBodyData.GetFlagFixed() == true)
            {
                //vecSamplePos[j] -= Vec3f(CParticleSystemConfig::m_timeStep, 0.0f, 0.0f);
                vecSamplePos[j] -= Vec3f(0.005f, 0.0f, 0.0f);
                if (vecSamplePos[j][0] < CParticleSystemConfig::m_cubicPosMin[0])
                {
                    vecSamplePos[j][0] += CParticleSystemConfig::m_cubicPosMax[0] - CParticleSystemConfig::m_cubicPosMin[0];
                }
                continue;
            }
            vecSamplePos[j] = pj + pv;
        }
    }
}

void CParticleSystemSyn::UpdateOutputGroupViaOptimization(CParticleSystem& outputGroup, int boundaryType /* = 0 */)
{
    ResetCoeffMatrix(outputGroup);
    const int numOfUnknownsTotal = outputGroup.GetNumOfSamples();
    const int numOfSoftBodies = outputGroup.GetNumOfSoftBodies();
    const int numOfSamplesPerData = outputGroup.GetParticleData(0).GetNumOfSamples();
    const Flt spatialWt = CParticleSystemConfig::m_spatialWt;
    const Flt shapeWt = CParticleSystemConfig::m_shapeWt;
    const Flt neighWt = CParticleSystemConfig::m_neighWt;
    const Flt sigma = CParticleSystemConfig::m_gaussianSigma / (CParticleSystemConfig::m_neighDist * CParticleSystemConfig::m_neighDist);
    if (neighWt > 0.0f || spatialWt > 0.0f)
    {
        outputGroup.SetNeighboringSamples(CParticleSystemConfig::m_neighDist);
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
    // Set the motion term...
    for (int n = 0; n < numOfUnknownsTotal; n++)
    {
        int softBodyIdx = n / numOfSamplesPerData;
        int sampleIdx = n % numOfSamplesPerData;
        CParticleData& softBodyData = outputGroup.GetParticleData(softBodyIdx);
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        Vec3f pos = vecSamplePos[sampleIdx];
        vecPx[n] += pos[0];
        vecPy[n] += pos[1];
        vecPz[n] += pos[2];
        m_vecCx[n] += pos[0];
        m_vecCy[n] += pos[1];
        m_vecCz[n] += pos[2];
        if (softBodyData.GetFlagBoundary() == true)
        {
            UpdateCoeffMatDiagonalVal(n, CParticleSystemConfig::m_boundaryWt, pos);
        }
    }
    FinishCoeffsAccordingToBoundary();
    SampleRepulsion();
    Flt distSum = 0.0f;
    Flt distMax = -1e10;
#pragma omp parallel for
    // Set the texture term...
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        CParticleData& softBodyData = outputGroup.GetParticleData(i);
        int startIdx = i * numOfSamplesPerData;
        CSampleMatchingResult matchingRst1;
        CSampleMatchingResult matchingRst2;
#ifdef MORPHABLE_PARTICLE_SYSTEM
        vector<CSampleMatchingResult> vecMatchingRst(2);
        Flt distMin = GetNearestInputSoftBodiesMorph(softBodyData, vecMatchingRst);
        distSum += distMin;
        vector<Vec3f>& vecPos = softBodyData.GetVecSamplePos();
        Flt weight1 = (vecPos[0][0] < 0.0f) ? 1.0f : 0.0f;
        Flt weight2 = 1.0f - weight1;
        CParticleData& inputData1 = m_vecInputExemplar[0].GetParticleData(vecMatchingRst[0].m_sampleIdx);
        CParticleData& inputData2 = m_vecInputExemplar[1].GetParticleData(vecMatchingRst[1].m_sampleIdx);
        matchingRst1 = vecMatchingRst[0];
        matchingRst2 = vecMatchingRst[1];
#else
        Flt distMin = GetNearestInputSoftBodies(softBodyData, i, matchingRst1, matchingRst2);
        distSum += distMin;
        distMax = max(distMax, distMin);
        Flt weight1 = matchingRst2.m_distSum;
        Flt weight2 = matchingRst1.m_distSum;
        if (distMin <= 0.0f)
        {
            weight1 = 1.0f;
            weight2 = 1.0f;
        }
        else
        {
            weight1 *= 1.0f / distMin;
            weight2 *= 1.0f / distMin;
        }
        CParticleData& inputData1 = m_inputSequence[matchingRst1.m_frameIdx].GetParticleData(matchingRst1.m_sampleIdx);
        CParticleData& inputData2 = m_inputSequence[matchingRst2.m_frameIdx].GetParticleData(matchingRst2.m_sampleIdx);
#endif
        // sample pairs between different elements...
        vector<CNeighboringSample>& vecNeighOutput = outputGroup.GetParticleData(i).GetNeighboringSamples();
        vector<CNeighboringSample>& vecNeighInput1 = inputData1.GetNeighboringSamples();
        vector<CNeighboringSample>& vecNeighInput2 = inputData2.GetNeighboringSamples();
        if (neighWt > 0.0f)
        {
            for (int j = 0; j < int(vecNeighOutput.size()); j++)
            {
                int idxi = startIdx + vecNeighOutput[j].m_sampleIdx;
                CNeighboringSample& neighInput1 = vecNeighInput1[matchingRst1.m_matchingIndices[j]];
                CNeighboringSample& neighInput2 = vecNeighInput2[matchingRst2.m_matchingIndices[j]];
                int idxj = vecNeighOutput[j].m_neighBodyIdx * numOfSamplesPerData + vecNeighOutput[j].m_neighSampleIdx;
                Vec3f prSrc1 = neighInput1.m_neighPr;
                Vec3f prSrc2 = neighInput2.m_neighPr;
                Flt neighWtNew = neighWt * exp(sigma * mag2(vecNeighOutput[j].m_neighPr));
                Flt wt1 = neighWt * exp(sigma * mag2(prSrc1));
                Flt wt2 = neighWt * exp(sigma * mag2(prSrc2));
                UpdateCoeffMatPairVals(idxi, idxj, (neighWtNew + wt1) * 0.5f * weight1, prSrc1);
                UpdateCoeffMatPairVals(idxi, idxj, (neighWtNew + wt2) * 0.5f * weight2, prSrc2);
            }
        }
    }
    cout << "distSum = " << distSum << endl;
    //cout << "distMax = " << distMax << endl;
    // Update via least squares...
    vector<Flt> vecPxNew = machy_math::GetSolution(m_vecPtrCoeffMatrix[0], m_vecCx);
    vector<Flt> vecPyNew = machy_math::GetSolution(m_vecPtrCoeffMatrix[1], m_vecCy);
    vector<Flt> vecPzNew = machy_math::GetSolution(m_vecPtrCoeffMatrix[2], m_vecCz);
    for (int n = 0; n < numOfUnknownsTotal; n++)
    {
        Vec3f pos = Vec3f(vecPxNew[n], vecPyNew[n], vecPzNew[n]);
        int softBodyIdx = n / numOfSamplesPerData;
        int sampleIdx = n % numOfSamplesPerData;
        CParticleData& softBodyData = outputGroup.GetParticleData(softBodyIdx);
        if (softBodyData.GetFlagFixed() == true)
        {
            continue;
        }
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        vecSamplePos[sampleIdx] = pos;
    }
}

void CParticleSystemSyn::UpdateOutputGroupViaOptimization3(CParticleSystem& outputGroup, CBoundaryConstraint* ptrBoundaryConstraint)
{
    int inputIdx = ptrBoundaryConstraint->GetInputIdx();
    const int numOfUnknownsTotal = outputGroup.GetNumOfSamples();
    const int numOfSoftBodies = outputGroup.GetNumOfSoftBodies();
    const int numOfSamplesPerData = outputGroup.GetParticleData(0).GetNumOfSamples();
    const Flt spatialWt = CParticleSystemConfig::m_spatialWt;
    const Flt shapeWt = CParticleSystemConfig::m_shapeWt;
    const Flt neighWt = CParticleSystemConfig::m_neighWt;
    const Flt boundaryWt = CParticleSystemConfig::m_boundaryWt;
    const Flt sigma = CParticleSystemConfig::m_gaussianSigma / (CParticleSystemConfig::m_neighDist * CParticleSystemConfig::m_neighDist);
    if (neighWt > 0.0f || spatialWt > 0.0f)
    {
        outputGroup.SetNeighboringSamples(CParticleSystemConfig::m_neighDist);
    }
    vector<Vec3f> vecPosSum(numOfUnknownsTotal, Vec3f(0.0f, 0.0f, 0.0f));
    vector<Flt> vecWtSum(numOfUnknownsTotal, 0.0f);
    // Set the motion term...
    for (int n = 0; n < numOfUnknownsTotal; n++)
    {
        int softBodyIdx = n / numOfSamplesPerData;
        int sampleIdx = n % numOfSamplesPerData;
        CParticleData& softBodyData = outputGroup.GetParticleData(softBodyIdx);
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        Vec3f pos = vecSamplePos[sampleIdx];
        vecPosSum[n] += pos;
        vecWtSum[n] += 1.0f;
    }
    Flt distSum = 0.0f;
    Flt distMax = -1e10;
#pragma omp parallel for
    // Set the texture term...
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        CParticleData& softBodyData = outputGroup.GetParticleData(i);
        int startIdx = i * numOfSamplesPerData;
        CSampleMatchingResult matchingRst1;
        CSampleMatchingResult matchingRst2;
#ifdef MORPHABLE_PARTICLE_SYSTEM
        vector<CSampleMatchingResult> vecMatchingRst(2);
        Flt distMin = GetNearestInputSoftBodiesMorph(softBodyData, vecMatchingRst);
        distSum += distMin;
        vector<Vec3f>& vecPos = softBodyData.GetVecSamplePos();
        Flt weight1 = (vecPos[0][0] < 0.0f) ? 1.0f : 0.0f;
        Flt weight2 = 1.0f - weight1;
        CParticleData& inputData1 = m_vecInputExemplar[0].GetParticleData(vecMatchingRst[0].m_sampleIdx);
        CParticleData& inputData2 = m_vecInputExemplar[1].GetParticleData(vecMatchingRst[1].m_sampleIdx);
        matchingRst1 = vecMatchingRst[0];
        matchingRst2 = vecMatchingRst[1];
#else
        Flt distMin = GetNearestInputSoftBodies(softBodyData, inputIdx, i, matchingRst1, matchingRst2);
        distSum += distMin;
        distMax = max(distMax, distMin);
        Flt weight1 = matchingRst2.m_distSum;
        Flt weight2 = matchingRst1.m_distSum;
        if (distMin <= 0.0f)
        {
            weight1 = 1.0f;
            weight2 = 1.0f;
        }
        else
        {
            weight1 *= 1.0f / distMin;
            weight2 *= 1.0f / distMin;
        }
        CParticleData& inputData1 = m_vecInputExemplar[matchingRst1.m_frameIdx].GetParticleData(matchingRst1.m_sampleIdx);
        CParticleData& inputData2 = m_vecInputExemplar[matchingRst2.m_frameIdx].GetParticleData(matchingRst2.m_sampleIdx);
#endif
        // sample pairs between different elements...
        vector<CNeighboringSample>& vecNeighOutput = outputGroup.GetParticleData(i).GetNeighboringSamples();
        vector<CNeighboringSample>& vecNeighInput1 = inputData1.GetNeighboringSamples();
        vector<CNeighboringSample>& vecNeighInput2 = inputData2.GetNeighboringSamples();
        if (neighWt > 0.0f)
        {
            for (int j = 0; j < int(vecNeighOutput.size()); j++)
            {
                int idxi = startIdx + vecNeighOutput[j].m_sampleIdx;
                CNeighboringSample& neighInput1 = vecNeighInput1[matchingRst1.m_matchingIndices[j]];
                CNeighboringSample& neighInput2 = vecNeighInput2[matchingRst2.m_matchingIndices[j]];
                int idxj = vecNeighOutput[j].m_neighBodyIdx * numOfSamplesPerData + vecNeighOutput[j].m_neighSampleIdx;
                Vec3f prSrc1 = neighInput1.m_neighPr;
                Vec3f prSrc2 = neighInput2.m_neighPr;
                Flt neighWtNew = neighWt * exp(sigma * mag2(vecNeighOutput[j].m_neighPr));
                Flt wt1 = neighWt * exp(sigma * mag2(prSrc1));
                Flt wt2 = neighWt * exp(sigma * mag2(prSrc2));
                Vec3f posPredict1 = outputGroup.GetParticleData(idxj).GetVecSamplePos()[0] + prSrc1;
                Vec3f posPredict2 = outputGroup.GetParticleData(idxi).GetVecSamplePos()[0] - prSrc1;
                Flt w1 = (neighWtNew + wt1) * 0.5f * weight1;
                Flt w2 = (neighWtNew + wt2) * 0.5f * weight2;
                if (boundaryWt > 0.0f)
                {
                    Vec3f posPredict1new = posPredict1;
                    Vec3f posPredict2new = posPredict2;
                    w1 = ptrBoundaryConstraint->ApplyBoundaryConstraintNew(posPredict1, posPredict1new) ? (w1 * boundaryWt) : w1;
                    w2 = ptrBoundaryConstraint->ApplyBoundaryConstraintNew(posPredict2, posPredict2new) ? (w2 * boundaryWt) : w2;
                    posPredict1 = posPredict1new;
                    posPredict2 = posPredict2new;
                }
                vecPosSum[idxi] += posPredict1 * w1;
                vecWtSum[idxi] += w1;
                vecPosSum[idxj] += posPredict2 * w2;
                vecWtSum[idxj] += w2;
            }
        }
    }
    cout << "distSum = " << distSum << endl;
    //cout << "distMax = " << distMax << endl;
    for (int n = 0; n < numOfUnknownsTotal; n++)
    {
        Vec3f pos = vecPosSum[n] / vecWtSum[n];
        int softBodyIdx = n / numOfSamplesPerData;
        int sampleIdx = n % numOfSamplesPerData;
        CParticleData& softBodyData = outputGroup.GetParticleData(softBodyIdx);
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        vecSamplePos[sampleIdx] = pos;
    }
}

void CParticleSystemSyn::ResetCoeffMatrix(CParticleSystem& outputGroup)
{
    const int numOfUnknownsTotal = outputGroup.GetNumOfSamples();
    const int numOfSamplesPerData = outputGroup.GetParticleData(0).GetNumOfSamples();
    if (m_vecPtrCoeffMatrix.empty() == true)
    {
        m_vecPtrCoeffMatrix.resize(3);
    }
    for (int n = 0; n < 3; n++)
    {
        m_vecPtrCoeffMatrix[n] = Eigen::MatrixXf::Identity(numOfUnknownsTotal, numOfUnknownsTotal);
    }
}

void CParticleSystemSyn::FinishCoeffsAccordingToBoundary()
{
    const Flt boundaryWt = CParticleSystemConfig::m_boundaryWt;
    const Flt boundaryThresh = CParticleSystemConfig::m_boundaryThresh;
    if (boundaryWt <= 0.0f) return;
    Vec3f boundaryThreshCubic = Vec3f(boundaryThresh, boundaryThresh, boundaryThresh);
    const Vec3f cubicPosMin = CParticleSystemConfig::m_cubicPosMin - boundaryThreshCubic;
    const Vec3f cubicPosMax = CParticleSystemConfig::m_cubicPosMax + boundaryThreshCubic;
    vector<int> vecUpdateCount(3, 0);
    for (int i = 0; i < m_outputGroup.GetNumOfSoftBodies(); i++)
    {
        CParticleData& softBodyData = m_outputGroup.GetParticleData(i);
        int startIdx = i * softBodyData.GetNumOfSamples();
        for (int j = 0; j < softBodyData.GetNumOfSamples(); j++)
        {
            Vec3f pos = softBodyData.GetVecSamplePos()[j];
            Vec3f posBoundary = pos;
            vector<bool> updateFlags(3, false);
            for (int k = 0; k < 3; k++)
            {
                if (pos[k] < cubicPosMin[k])
                {
                    posBoundary[k] = 2.0f * cubicPosMin[k] - pos[k];
                    updateFlags[k] = true;
                    vecUpdateCount[k]++;
                }
                else if (pos[k] > cubicPosMax[k])
                {
                    posBoundary[k] = 2.0f * cubicPosMax[k] - pos[k];
                    updateFlags[k] = true;
                    vecUpdateCount[k]++;
                }
            }
            for (int k = 0; k < 3; k++)
            {
                int idx = startIdx + j;
                if (updateFlags[k] == true)
                {
                    UpdateCoeffMatDiagonalVal(idx, boundaryWt, posBoundary, k);
                }
            }
        }
    }
    cout << "Boundary update count: " << vecUpdateCount[0] << ", " << vecUpdateCount[1] << ", " << vecUpdateCount[2] << endl;
}

void CParticleSystemSyn::UpdateOutputFrameViaBlending(int frameIdx, CBoundaryConstraint* ptrBoundaryConstraint1, CBoundaryConstraint* ptrBoundaryConstraint2)
{
    int frameIdx1 = ptrBoundaryConstraint1->GetFrameIdx();
    int frameIdx2 = ptrBoundaryConstraint2->GetFrameIdx();
    Flt frameWt1 = (frameIdx2 - frameIdx) / Flt(frameIdx2 - frameIdx1);
    frameWt1 = min(frameWt1, 1.0f);
    frameWt1 = max(frameWt1, 0.0f);
    Flt frameWt2 = 1.0f - frameWt1;
    int inputIdx1 = ptrBoundaryConstraint1->GetInputIdx();
    int inputIdx2 = ptrBoundaryConstraint2->GetInputIdx();
    vector<Flt> vecFrameWt;
    vecFrameWt.push_back(frameWt1);
    vecFrameWt.push_back(frameWt2);
    vector<int> vecInputIdx;
    vecInputIdx.push_back(inputIdx1);
    vecInputIdx.push_back(inputIdx2);
    CParticleSystem& outputGroup = m_outputSequence[frameIdx];
    const int numOfUnknownsTotal = outputGroup.GetNumOfSamples();
    const int numOfSoftBodies = outputGroup.GetNumOfSoftBodies();
    const int numOfSamplesPerData = outputGroup.GetParticleData(0).GetNumOfSamples();
    const Flt spatialWt = CParticleSystemConfig::m_spatialWt;
    const Flt shapeWt = CParticleSystemConfig::m_shapeWt;
    const Flt neighWt = CParticleSystemConfig::m_neighWt;
    const Flt boundaryWt = CParticleSystemConfig::m_boundaryWt;
    const Flt sigma = CParticleSystemConfig::m_gaussianSigma / (CParticleSystemConfig::m_neighDist * CParticleSystemConfig::m_neighDist);
    if (neighWt > 0.0f || spatialWt > 0.0f)
    {
        outputGroup.SetNeighboringSamples(CParticleSystemConfig::m_neighDist);
    }
    CBoundaryConstraint* ptrBoundaryConstraint = NULL;
    for (int i = 0; i < int(m_vecBoundaryConstraint.size()); i++)
    {
        if (m_vecBoundaryConstraint[i].GetFrameIdx() == frameIdx)
        {
            ptrBoundaryConstraint = &(m_vecBoundaryConstraint[i]);
            break;
        }
    }
    int temporalWindowSize = CParticleSystemConfig::m_temporalWindowSize;
    temporalWindowSize = min(temporalWindowSize, frameIdx);
    temporalWindowSize = min(temporalWindowSize, int(m_outputSequence.size()) - frameIdx - 1);
    vector<Vec3f> vecPosSum(numOfUnknownsTotal, Vec3f(0.0f, 0.0f, 0.0f));
    vector<Flt> vecWtSum(numOfUnknownsTotal, 0.0f);
#pragma omp parallel for
    // Set the motion term...
    for (int n = 0; n < numOfUnknownsTotal; n++)
    {
        for (int j = -temporalWindowSize; j <= temporalWindowSize; j++)
        {
            CParticleSystem& outputGroupTmp = m_outputSequence[frameIdx + j];
            Vec3f& posTmp = outputGroupTmp.GetParticleData(n).GetVecSamplePos()[0];
            Flt temporalWt = CParticleSystemConfig::m_temporalWt;
            vecPosSum[n] += posTmp * temporalWt;
            vecWtSum[n] += temporalWt;
        }
    }
    Flt distSum = 0.0f;
    Flt distMax = -1e10;
#pragma omp parallel for
    // Set the texture term...
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        CParticleData& softBodyData = outputGroup.GetParticleData(i);
        int startIdx = i * numOfSamplesPerData;
        for (int n = 0; n < int(vecInputIdx.size()); n++)
        {
            Flt frameWt = vecFrameWt[n];
            int inputIdx = vecInputIdx[n];
            CSampleMatchingResult matchingRst1;
            CSampleMatchingResult matchingRst2;
            Flt distMin = GetNearestInputSoftBodies(softBodyData, inputIdx, i, matchingRst1, matchingRst2);
            distSum += distMin;
            distMax = max(distMax, distMin);
            Flt weight1 = matchingRst2.m_distSum;
            Flt weight2 = matchingRst1.m_distSum;
            if (distMin <= 0.0f)
            {
                weight1 = 1.0f;
                weight2 = 1.0f;
            }
            else
            {
                weight1 *= 1.0f / distMin;
                weight2 *= 1.0f / distMin;
            }
            CParticleData& inputData1 = m_vecInputExemplar[matchingRst1.m_frameIdx].GetParticleData(matchingRst1.m_sampleIdx);
            CParticleData& inputData2 = m_vecInputExemplar[matchingRst2.m_frameIdx].GetParticleData(matchingRst2.m_sampleIdx);
            // sample pairs between different elements...
            vector<CNeighboringSample>& vecNeighOutput = outputGroup.GetParticleData(i).GetNeighboringSamples();
            vector<CNeighboringSample>& vecNeighInput1 = inputData1.GetNeighboringSamples();
            vector<CNeighboringSample>& vecNeighInput2 = inputData2.GetNeighboringSamples();
            if (neighWt > 0.0f)
            {
                for (int j = 0; j < int(vecNeighOutput.size()); j++)
                {
                    int idxi = startIdx + vecNeighOutput[j].m_sampleIdx;
                    CNeighboringSample& neighInput1 = vecNeighInput1[matchingRst1.m_matchingIndices[j]];
                    CNeighboringSample& neighInput2 = vecNeighInput2[matchingRst2.m_matchingIndices[j]];
                    int idxj = vecNeighOutput[j].m_neighBodyIdx * numOfSamplesPerData + vecNeighOutput[j].m_neighSampleIdx;
                    Vec3f prSrc1 = neighInput1.m_neighPr;
                    Vec3f prSrc2 = neighInput2.m_neighPr;
                    Flt neighWtNew = neighWt * exp(sigma * mag2(vecNeighOutput[j].m_neighPr));
                    Flt wt1 = neighWt * exp(sigma * mag2(prSrc1));
                    Flt wt2 = neighWt * exp(sigma * mag2(prSrc2));
                    Vec3f posPredict1 = outputGroup.GetParticleData(idxj).GetVecSamplePos()[0] + prSrc1;
                    Vec3f posPredict2 = outputGroup.GetParticleData(idxi).GetVecSamplePos()[0] - prSrc1;
                    Flt w1 = (neighWtNew + wt1) * 0.5f * weight1;
                    Flt w2 = (neighWtNew + wt2) * 0.5f * weight2;
                    if (ptrBoundaryConstraint != NULL)
                    {
                        Vec3f posPredict1new = posPredict1;
                        Vec3f posPredict2new = posPredict2;
                        w1 = ptrBoundaryConstraint->ApplyBoundaryConstraintNew(posPredict1, posPredict1new) ? (w1 * boundaryWt) : w1;
                        w2 = ptrBoundaryConstraint->ApplyBoundaryConstraintNew(posPredict1, posPredict1new) ? (w2 * boundaryWt) : w2;
                        posPredict1 = posPredict1new;
                        posPredict2 = posPredict2new;
                    }
                    wt1 *= frameWt;
                    wt2 *= frameWt;
                    vecPosSum[idxi] += posPredict1 * w1;
                    vecWtSum[idxi] += w1;
                    vecPosSum[idxj] += posPredict2 * w2;
                    vecWtSum[idxj] += w2;
                }
            } // End-If
        }
    }
    cout << "distSum = " << distSum << endl;
    for (int n = 0; n < numOfUnknownsTotal; n++)
    {
        Vec3f pos = vecPosSum[n] / vecWtSum[n];
        int softBodyIdx = n / numOfSamplesPerData;
        int sampleIdx = n % numOfSamplesPerData;
        CParticleData& softBodyData = outputGroup.GetParticleData(softBodyIdx);
        vector<Vec3f>& vecSamplePos = softBodyData.GetVecSamplePos();
        vecSamplePos[sampleIdx] = pos;
    }
}

bool CParticleSystemSyn::ApplyBoundaryCondition(Vec3f& posOld, Vec3f& posNew)
{
    const Flt boundaryWt = CParticleSystemConfig::m_boundaryWt;
    const Flt boundaryThresh = CParticleSystemConfig::m_boundaryThresh;
    if (boundaryWt <= 0.0f) return false;
    Vec3f boundaryThreshCubic = Vec3f(boundaryThresh, boundaryThresh, boundaryThresh);
#ifdef FLOW_SYNTHESIS // 01/07/2013
    boundaryThreshCubic[0] += CParticleSystemConfig::m_timeStep * m_stepCount;
#endif
    const Vec3f cubicPosMin = CParticleSystemConfig::m_cubicPosMin - boundaryThreshCubic;
    const Vec3f cubicPosMax = CParticleSystemConfig::m_cubicPosMax + boundaryThreshCubic;
    Vec3f pos = posOld;
    Vec3f posBoundary = pos;
    bool flag = false;
    for (int k = 0; k < 3; k++)
    {
        if (pos[k] < cubicPosMin[k])
        {
            posBoundary[k] = 2.0f * cubicPosMin[k] - pos[k];
            flag = true;
        }
        else if (pos[k] > cubicPosMax[k])
        {
            posBoundary[k] = 2.0f * cubicPosMax[k] - pos[k];
            flag = true;
        }
    }
    if (flag == true)
    {
        posNew = posBoundary;
    }
    return flag;
}

bool CParticleSystemSyn::ApplyBoundaryCondition2(Vec3f& posOld, Vec3f& posNew)
{
    const Flt boundaryWt = CParticleSystemConfig::m_boundaryWt;
    const Flt boundaryThresh = CParticleSystemConfig::m_boundaryThresh;
    if (boundaryWt <= 0.0f) return false;
    Vec3f boundaryThreshCubic = Vec3f(boundaryThresh, boundaryThresh, boundaryThresh);
    const Vec3f cubicPosMin = CParticleSystemConfig::m_cubicPosMin - boundaryThreshCubic;
    const Vec3f cubicPosMax = CParticleSystemConfig::m_cubicPosMax + boundaryThreshCubic;
    Vec3f center = (cubicPosMin + cubicPosMax) * 0.5f;
    Flt area = (cubicPosMax[0] - cubicPosMin[0]) * (cubicPosMax[1] - cubicPosMin[1]);
    Flt r = sqrt(area / M_PI);
    r -= 0.05f;
    Vec3f pos = posOld;
    bool flag = false;
    Flt distSqr = dist2(pos, center);
    if (distSqr > r * r)
    {
        flag = true;
        Vec3f pr = pos - center;
        normalize(pr);
        pr *= r;
        posNew = center + pr;
    }
    return flag;
}

Flt CParticleSystemSyn::GetNearestInputSoftBodies(CParticleData& softBodyOut, int outputIdx, CSampleMatchingResult& matchingRst1, CSampleMatchingResult& matchingRst2)
{
    Flt distMin1 = 1e10;
    Flt distMin2 = 1e10;
    for (int i = 0; i < int(m_inputSequence.size()); i++)
    {
        CParticleSystem& inputGroup = m_inputSequence[i];
        for (int j = 0; j < inputGroup.GetNumOfSoftBodies(); j++)
        {
            CParticleData& softBodyInput = inputGroup.GetParticleData(j);
            vector<int> matchIndicesTmp;
            Flt distTmp = SoftBodyNeighDistance(softBodyOut, softBodyInput, matchIndicesTmp, CParticleSystemConfig::m_neighWt);
            if (distTmp < distMin1)
            {
                distMin2 = distMin1;
                matchingRst2 = matchingRst1;
                distMin1 = distTmp;
                matchingRst1.m_frameIdx = i;
                matchingRst1.m_sampleIdx = j;
                matchingRst1.m_matchingIndices = matchIndicesTmp;
                matchingRst1.m_distSum = distTmp;
            }
            else if (distTmp < distMin2)
            {
                distMin2 = distTmp;
                matchingRst2.m_frameIdx = i;
                matchingRst2.m_sampleIdx = j;
                matchingRst2.m_matchingIndices = matchIndicesTmp;
                matchingRst2.m_distSum = distTmp;
            }
        }
    }
    if (CParticleSystemConfig::m_flagSmoothSynthesis == false)
    {
        matchingRst1.m_weight = 1.0f;
        matchingRst2.m_weight = 0.0f;
        return distMin1 + distMin2;
    }
    Flt a = 0.0f;
    Flt b = 0.0f;
    vector<CNeighboringSample>& vecNeighOut = softBodyOut.GetNeighboringSamples();
    vector<CNeighboringSample>& vecNeigh1 = m_inputSequence[matchingRst1.m_frameIdx].GetParticleData(matchingRst1.m_sampleIdx).GetNeighboringSamples();
    vector<CNeighboringSample>& vecNeigh2 = m_inputSequence[matchingRst2.m_frameIdx].GetParticleData(matchingRst2.m_sampleIdx).GetNeighboringSamples();
    for (int i = 0; i < int(vecNeighOut.size()); i++)
    {
        Vec3f pr0 = vecNeighOut[i].m_neighPr;
        Vec3f pr1 = vecNeigh1[i].m_neighPr;
        Vec3f pr2 = vecNeigh2[i].m_neighPr;
        a += dist2(pr1, pr2);
        b += dot(pr1 - pr2, pr0 - pr2);
    }
    Flt weight1 = b / a;
    weight1 = min(weight1, 1.0f);
    weight1 = max(weight1, 0.0f);
    matchingRst1.m_weight = weight1;
    matchingRst2.m_weight = 1.0f - weight1;
    return distMin1 + distMin2;
}

Flt CParticleSystemSyn::GetNearestInputSoftBodies(CParticleData& softBodyOut, int inputIdx, int outputIdx, CSampleMatchingResult& matchingRst1, CSampleMatchingResult& matchingRst2)
{
    Flt distMin1 = 1e10;
    Flt distMin2 = 1e10;
    //for ( int i=0; i<int(m_inputSequence.size()); i++ )
    {
        CParticleSystem& inputGroup = m_vecInputExemplar[inputIdx];
        for (int j = 0; j < inputGroup.GetNumOfSoftBodies(); j++)
        {
            CParticleData& softBodyInput = inputGroup.GetParticleData(j);
            vector<int> matchIndicesTmp;
            Flt distTmp = SoftBodyNeighDistance(softBodyOut, softBodyInput, matchIndicesTmp, CParticleSystemConfig::m_neighWt);
            if (distTmp < distMin1)
            {
                distMin2 = distMin1;
                matchingRst2 = matchingRst1;
                distMin1 = distTmp;
                matchingRst1.m_frameIdx = inputIdx;
                matchingRst1.m_sampleIdx = j;
                matchingRst1.m_matchingIndices = matchIndicesTmp;
                matchingRst1.m_distSum = distTmp;
            }
            else if (distTmp < distMin2)
            {
                distMin2 = distTmp;
                matchingRst2.m_frameIdx = inputIdx;
                matchingRst2.m_sampleIdx = j;
                matchingRst2.m_matchingIndices = matchIndicesTmp;
                matchingRst2.m_distSum = distTmp;
            }
        }
    }
    if (CParticleSystemConfig::m_flagSmoothSynthesis == false)
    {
        matchingRst1.m_weight = 1.0f;
        matchingRst2.m_weight = 0.0f;
        return distMin1 + distMin2;
    }
    Flt a = 0.0f;
    Flt b = 0.0f;
    vector<CNeighboringSample>& vecNeighOut = softBodyOut.GetNeighboringSamples();
    vector<CNeighboringSample>& vecNeigh1 = m_vecInputExemplar[matchingRst1.m_frameIdx].GetParticleData(matchingRst1.m_sampleIdx).GetNeighboringSamples();
    vector<CNeighboringSample>& vecNeigh2 = m_vecInputExemplar[matchingRst2.m_frameIdx].GetParticleData(matchingRst2.m_sampleIdx).GetNeighboringSamples();
    for (int i = 0; i < int(vecNeighOut.size()); i++)
    {
        Vec3f pr0 = vecNeighOut[i].m_neighPr;
        Vec3f pr1 = vecNeigh1[i].m_neighPr;
        Vec3f pr2 = vecNeigh2[i].m_neighPr;
        a += dist2(pr1, pr2);
        b += dot(pr1 - pr2, pr0 - pr2);
    }
    Flt weight1 = b / a;
    weight1 = min(weight1, 1.0f);
    weight1 = max(weight1, 0.0f);
    matchingRst1.m_weight = weight1;
    matchingRst2.m_weight = 1.0f - weight1;
    return distMin1 + distMin2;
}

Flt CParticleSystemSyn::GetNearestInputSoftBodiesMorph(CParticleData& softBodyOut, vector<CSampleMatchingResult>& vecMatchingRst)
{
    Flt distSum = 0.0f;
    for (int i = 0; i < int(m_vecInputExemplar.size()); i++)
    {
        Flt distMin = 1e10;
        CParticleSystem& inputExemplar = m_vecInputExemplar[i];
        CSampleMatchingResult& matchingRst = vecMatchingRst[i];
        for (int j = 0; j < inputExemplar.GetNumOfSoftBodies(); j++)
        {
            CParticleData& softBodyInput = inputExemplar.GetParticleData(j);
            vector<int> matchIndicesTmp;
            Flt distTmp = SoftBodyNeighDistance(softBodyOut, softBodyInput, matchIndicesTmp, CParticleSystemConfig::m_neighWt);
            if (distTmp < distMin)
            {
                distMin = distTmp;
                matchingRst.m_frameIdx = i;
                matchingRst.m_sampleIdx = j;
                matchingRst.m_matchingIndices = matchIndicesTmp;
                matchingRst.m_distSum = distTmp;
            }
        }
        distSum += distMin;
    }
    return distSum;
}

Flt CParticleSystemSyn::SoftBodyNeighDistance(CParticleData& softBody1, CParticleData& softBody2, vector<int>& matchIndices, Flt weight)
{
    if (weight <= 0.0f)
    {
        return 0.0f;
    }
    const Flt sigma = CParticleSystemConfig::m_gaussianSigma / (CParticleSystemConfig::m_neighDist * CParticleSystemConfig::m_neighDist);
    vector<CNeighboringSample>& vecNeigh1 = softBody1.GetNeighboringSamples();
    vector<CNeighboringSample>& vecNeigh2 = softBody2.GetNeighboringSamples();
    int size1 = int(vecNeigh1.size());
    int size2 = int(vecNeigh2.size());
    if (size1 > size2)
    {
        return 1e10;
    }
    matchIndices.resize(size1);
    vector<double> distMatrix(size1 * size2);
    int cnt = 0;
    for (int j = 0; j < size2; j++)
    {
        for (int i = 0; i < size1; i++)
        {
            Vec3f pr1 = vecNeigh1[i].m_neighPr;
            Vec3f pr2 = vecNeigh2[j].m_neighPr;
            Flt wt1 = exp(sigma * mag2(pr1));
            Flt wt2 = exp(sigma * mag2(pr2));
            Vec3f prr = pr1 - pr2;
            Flt distTmp = mag2(prr);
            distTmp *= (wt1 + wt2) * 0.5f;
            distMatrix[cnt++] = distTmp;
        }
    }
    vector<int> assignment(size1);
    double cost = 0.0f;
    CHungarianAlgorithm ha;
    ha.FindOptimalAssignment(assignment, &cost, distMatrix, size1, size2);
    for (int i = 0; i < size1; i++)
    {
        matchIndices[i] = assignment[i] - 1;
    }
    return (cost * weight);
}

void CParticleSystemSyn::UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos)
{
    for (int n = 0; n < 3; n++)
    {
        m_vecPtrCoeffMatrix[n](idx, idx) += wt;
    }
    m_vecCx[idx] += wt * pos[0];
    m_vecCy[idx] += wt * pos[1];
    m_vecCz[idx] += wt * pos[2];
}

void CParticleSystemSyn::UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos, int dim)
{
    m_vecPtrCoeffMatrix[dim](idx, idx) += wt;
    switch (dim)
    {
    case 0:
        m_vecCx[idx] += wt * pos[0];
        break;
    case 1:
        m_vecCy[idx] += wt * pos[1];
        break;
    case 2:
        m_vecCz[idx] += wt * pos[2];
        break;
    default:
        break;
    }
}

void CParticleSystemSyn::UpdateCoeffMatPairVals(int idxi, int idxj, Flt wt, Vec3f pr)
{
    for (int n = 0; n < 3; n++)
    {
        m_vecPtrCoeffMatrix[n](idxi, idxi) += wt;
        m_vecPtrCoeffMatrix[n](idxj, idxj) += wt;
        m_vecPtrCoeffMatrix[n](idxi, idxj) -= wt;
        m_vecPtrCoeffMatrix[n](idxj, idxi) -= wt;
    }
    m_vecCx[idxi] += wt * pr[0];
    m_vecCy[idxi] += wt * pr[1];
    m_vecCz[idxi] += wt * pr[2];
    m_vecCx[idxj] -= wt * pr[0];
    m_vecCy[idxj] -= wt * pr[1];
    m_vecCz[idxj] -= wt * pr[2];
}

CParticleSystem CParticleSystemSyn::LoadInputExemplarFromTXT(const string& fileName)
{
    CParticleSystem inputGroup;
    ifstream fin(fileName.c_str());
    int numOfElementTypes; // useless
    fin >> numOfElementTypes;

    int numOfElements;
    fin >> numOfElements;

    for (int i = 0; i < numOfElements; i++)
    {
        int elementType;
        fin >> elementType;
        // The matrix (orientation, position) of the element...
        float mat[16];
        for (int j = 0; j < 16; j++)
        {
            fin >> mat[j];
        }
        Vec3f pos = Vec3f(mat[12], mat[13], 0.0f);
        CParticleData softBodyData;
        vector<Vec3f> vertices(1);
        vertices[0] = pos;
        softBodyData.SetVecSamplePos(vertices);
        softBodyData.SetPos(pos);
        inputGroup.AddParticle(softBodyData);
    }
    return inputGroup;
}

vector<CubicPatch> CParticleSystemSyn::GetCubicPatches(vector<CParticleSystem>& groupSequence)
{
    vector<CubicPatch> cubicPatches;
    int nonEmptyPatchCount = 0;
    for (int n = 0; n < int(groupSequence.size()); n++)
    {
        vector<CubicPatch> cubicPatchesTmp = GetCubicPatches(groupSequence[n]);
        for (int i = 0; i < int(cubicPatchesTmp.size()); i++)
        {
            if (int(cubicPatchesTmp[i].m_patchIndices.size()) > CParticleSystemConfig::m_minimumPatchSize)
            {
                if (cubicPatchesTmp[i].m_patchIndices.empty() == false) nonEmptyPatchCount++;
            }
            else
            {
                continue; // Skip patches smaller than the thresh!
            }
            cubicPatchesTmp[i].m_patchFrameIdx = n;
            cubicPatches.push_back(cubicPatchesTmp[i]);
        }
    }
    cout << nonEmptyPatchCount << " of " << cubicPatches.size() << " patches are non-empty!\n";
    return cubicPatches;
}

vector<CubicPatch> CParticleSystemSyn::GetCubicPatches(CParticleSystem& softBodyGroup)
{
    vector<CubicPatch> cubicPatches;
    Vec3f patchSize = CParticleSystemConfig::m_patchSize;
    Vec3f cubicMin = Vec3f(-1.0f, -1.0f, -1.0f);
    Vec3f cubicMax = Vec3f(1.0f, 1.0f, 1.0f);
    //cubicMin = Vec3f(-0.9f, -0.9f, -0.9f);
    cubicMin = Vec3f(-1.1f, -1.1f, -1.1f);
    Vec3i patchNum;
    for (int i = 0; i < 3; i++)
    {
        patchNum[i] = floor((cubicMax[i] - cubicMin[i]) / patchSize[i] + 0.5f);
    }
    Vec3f patchCen;
    for (int i = 0; i < patchNum[0]; i++)
    {
        patchCen[0] = cubicMin[0] + (i + 0.5f) * patchSize[0];
        for (int j = 0; j < patchNum[1]; j++)
        {
            patchCen[1] = cubicMin[1] + (j + 0.5f) * patchSize[1];
            for (int k = 0; k < patchNum[2]; k++)
            {
                patchCen[2] = cubicMin[2] + (k + 0.5f) * patchSize[2];
                vector<int> indices = GetCubicPatchIndices(softBodyGroup, patchCen, patchSize);
                CubicPatch patch;
                patch.m_patchCen = patchCen;
                patch.m_patchIndices = indices;
                patch.m_patchFrequency = 0;
                if (int(patch.m_patchIndices.size()) <= CParticleSystemConfig::m_minimumPatchSize)
                {
                    continue;
                }
                cubicPatches.push_back(patch);
            }
        }
    }
    return cubicPatches;
}

vector<int> CParticleSystemSyn::GetCubicPatchIndices(CParticleSystem& softBodyGroup, Vec3f patchCen, Vec3f patchSize)
{
    vector<int> indices;
    Vec3f patchSizeHalf = 0.5f * patchSize;
    for (int i = 0; i < softBodyGroup.GetNumOfSoftBodies(); i++)
    {
        Vec3f pi = softBodyGroup.GetParticleData(i).GetPos();
        if (abs(pi[0] - patchCen[0]) > patchSizeHalf[0])
        {
            continue;
        }
        if (abs(pi[1] - patchCen[1]) > patchSizeHalf[1])
        {
            continue;
        }
        if (abs(pi[2] - patchCen[2]) > patchSizeHalf[2])
        {
            continue;
        }
        indices.push_back(i);
    }
    return indices;
}

vector<CParticleSystem> CParticleSystemSyn::InterpolateParticleSystem(CParticleSystem& softBodyGroup1, CParticleSystem& softBodyGroup2, int sequenceLength)
{
    int numOfSoftBodies = softBodyGroup1.GetNumOfSoftBodies();
    vector<int> matchIndices(numOfSoftBodies);
    vector<double> distMatrix(numOfSoftBodies * numOfSoftBodies);
    int cnt = 0;
    for (int j = 0; j < numOfSoftBodies; j++)
    {
        for (int i = 0; i < numOfSoftBodies; i++)
        {
            Vec3f p1 = softBodyGroup1.GetParticleData(i).GetVecSamplePos()[0];
            Vec3f p2 = softBodyGroup2.GetParticleData(j).GetVecSamplePos()[0];
            Vec3f pr = p1 - Vec3f(p2[0], p2[1], p2[2]);
            Flt distTmp = mag2(pr);
            distMatrix[cnt++] = distTmp;
        }
    }
    vector<int> assignment(numOfSoftBodies);
    double cost = 0.0f;
    CHungarianAlgorithm ha;
    ha.FindOptimalAssignment(assignment, &cost, distMatrix, numOfSoftBodies, numOfSoftBodies);
    for (int i = 0; i < numOfSoftBodies; i++)
    {
        matchIndices[i] = assignment[i] - 1;
    }
    vector<CParticleSystem> sequence;
    for (int i = 0; i < sequenceLength; i++)
    {
        Flt wt2 = Flt(i) / Flt(sequenceLength - 1);
        Flt wt1 = 1.0f - wt2;
        CParticleSystem softBodyGroup = softBodyGroup1;
        for (int j = 0; j < numOfSoftBodies; j++)
        {
            Vec3f& p1 = softBodyGroup.GetParticleData(j).GetVecSamplePos()[0];
            Vec3f p2 = softBodyGroup2.GetParticleData(matchIndices[j]).GetVecSamplePos()[0];
            p1 = wt1 * p1 + wt2 * p2;
            softBodyGroup.GetParticleData(j).SetPos(p1);
        }
        sequence.push_back(softBodyGroup);
    }
    softBodyGroup2 = sequence.back();
    return sequence;
}

void CParticleSystemSyn::LoadBoundaryConstraints()
{
    m_vecBoundaryConstraint.clear();
    CBoundaryConstraint bc;
    bc.LoadBoundaryConstraintFromImage("S03.png");
    bc.SetPosMin(CParticleSystemConfig::m_cubicPosMin);
    bc.SetPosMax(CParticleSystemConfig::m_cubicPosMax);
    bc.SetFrameIdx(0);
    bc.SetInputIdx(0);
    m_vecBoundaryConstraint.push_back(bc);
    //CParticleSystemConfig::m_cubicPosMin[0] += 1.0f;
    //CParticleSystemConfig::m_cubicPosMax[0] += 1.0f;
    CBoundaryConstraint bc2;
    bc2.LoadBoundaryConstraintFromImage("G04.png");
    bc2.SetPosMin(CParticleSystemConfig::m_cubicPosMin);
    bc2.SetPosMax(CParticleSystemConfig::m_cubicPosMax);
    bc2.SetFrameIdx(CParticleSystemConfig::m_stepCountMax);
    bc2.SetInputIdx(1);
    m_vecBoundaryConstraint.push_back(bc2);
    CBoundaryConstraint bc3;
    bc3.LoadBoundaryConstraintFromImage("R02.png");
    bc3.SetPosMin(CParticleSystemConfig::m_cubicPosMin);
    bc3.SetPosMax(CParticleSystemConfig::m_cubicPosMax);
    bc3.SetFrameIdx(CParticleSystemConfig::m_stepCountMax * 2);
    bc3.SetInputIdx(2);
    m_vecBoundaryConstraint.push_back(bc3);
}

void CParticleSystemSyn::SetSequencesVelocity(vector<CParticleSystem>& vecSequence, bool flagTemporallyToroidal /* = false */)
{
    const int numOfFrames = int(vecSequence.size());
    for (int t = 1; t < numOfFrames; t++)
    {
        CParticleSystem& sequence1 = vecSequence[t];
        CParticleSystem& sequence0 = vecSequence[t - 1];
        for (int i = 0; i < sequence1.GetNumOfSoftBodies(); i++)
        {
            CParticleData& particle1 = sequence1.GetParticleData(i);
            CParticleData& particle0 = sequence0.GetParticleData(i);
            vector<Vec3f>& vecSample1 = particle1.GetVecSamplePos();
            vector<Vec3f>& vecSample0 = particle0.GetVecSamplePos();
            Vec3f vel = vecSample1[0] - vecSample0[0];
            particle1.SetVel(vel);
        }
    }
    // For the first frame...
    if (vecSequence.size() <= 1) return;
    if (flagTemporallyToroidal == true)
    {
        CParticleSystem& sequence1 = vecSequence[0];
        CParticleSystem& sequence0 = vecSequence[vecSequence.size() - 1];
        for (int i = 0; i < sequence1.GetNumOfSoftBodies(); i++)
        {
            CParticleData& particle1 = sequence1.GetParticleData(i);
            CParticleData& particle0 = sequence0.GetParticleData(i);
            vector<Vec3f>& vecSample1 = particle1.GetVecSamplePos();
            vector<Vec3f>& vecSample0 = particle0.GetVecSamplePos();
            Vec3f vel = vecSample1[0] - vecSample0[0];
            particle1.SetVel(vel);
        }
        return;
    }
    CParticleSystem& sequence1 = vecSequence[1];
    CParticleSystem& sequence0 = vecSequence[0];
    for (int i = 0; i < sequence1.GetNumOfSoftBodies(); i++)
    {
        CParticleData& particle1 = sequence1.GetParticleData(i);
        CParticleData& particle0 = sequence0.GetParticleData(i);
        vector<Vec3f>& vecSample1 = particle1.GetVecSamplePos();
        vector<Vec3f>& vecSample0 = particle0.GetVecSamplePos();
        Vec3f vel = vecSample1[0] - vecSample0[0];
        particle0.SetVel(vel);
    }
}

void CParticleSystemSyn::SampleRepulsion()
{
    const Flt repulsionWt = CParticleSystemConfig::m_repulsionWt;
    const Flt repulsionDist = CParticleSystemConfig::m_repulsionDist * 2.0f;
    const Flt repulsionDistSqr = repulsionDist * repulsionDist;
    if (repulsionWt <= 0.0f || repulsionDist <= 0.0f) return;
    const int numOfUnknownsTotal = m_outputGroup.GetNumOfSamples();
    const int numOfSoftBodies = m_outputGroup.GetNumOfSoftBodies();
    const int numOfSamplesPerData = m_outputGroup.GetParticleData(0).GetNumOfSamples();
    int repulsionCount = 0;
#pragma omp parallel for
    for (int i1 = 0; i1 < numOfSoftBodies; i1++)
    {
        CParticleData& softBody1 = m_outputGroup.GetParticleData(i1);
        vector<Vec3f>& vecSamplePos1 = softBody1.GetVecSamplePos();
        for (int i2 = i1 + 1; i2 < numOfSoftBodies; i2++)
        {
            CParticleData& softBody2 = m_outputGroup.GetParticleData(i2);
            vector<Vec3f>& vecSamplePos2 = softBody2.GetVecSamplePos();
            for (int j1 = 0; j1 < numOfSamplesPerData; j1++)
            {
                for (int j2 = 0; j2 < numOfSamplesPerData; j2++)
                {
                    Vec3f p1 = vecSamplePos1[j1];
                    Vec3f p2 = vecSamplePos2[j2];
                    if (dist2(p1, p2) < repulsionDistSqr)
                    {
                        Vec3f pr = p2 - p1;
                        normalize(pr);
                        pr *= repulsionDist;
                        int idx1 = i1 * numOfSamplesPerData + j1;
                        int idx2 = i2 * numOfSamplesPerData + j2;
                        UpdateCoeffMatPairVals(idx1, idx2, repulsionWt, -pr);
                        repulsionCount++;
                    }
                }
            }
        }
    }
    if (repulsionCount > 0)
    {
        cout << "Have repulsed " << repulsionCount << " sample pairs!\n";
    }
}
