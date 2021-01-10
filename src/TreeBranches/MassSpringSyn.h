#ifndef MASSSPRINGSYN_H
#define MASSSPRINGSYN_H

#include "../DynamicElement/MathFunc.h"
#include "MassSpringData.h"
#include "MassSpringSynConfig.h"

class CMassSpringSyn
{
public:
    CMassSpringSyn(const std::string& config_file_name);

    ~CMassSpringSyn();

    void ResetOutput(const std::string& config_file_name);

    void LoadInputData();

    void GatherInputNeighborsExtended();

    void InitializeOutput();

    void ResetCoeffMatrix();

    void UpdateOutput();

    void UpdateStrandsExtended();

    void CollisionResponse(Eigen::VectorXf& vecCx, Eigen::VectorXf& vecCy, Eigen::VectorXf& vecCz);

    void RenderInput();

    void RenderOutput();

    void RestartSynthesis();

    bool LoadMassSequenceFromCSV(vector<MassSpringSequence>& vecSequence, const std::string& fileName);

    bool SaveMassSequenceAsCSV(const vector<MassSpringSequence>& vecSequence, const std::string& fileName);

    inline int GetStepCount() { return m_stepCount; }

private:
    MassNeighborhoodExtended GetMassNeighborhoodExtended(vector<CMassSpringData>& strands, int strandIdx, int massIdx, Flt neighDist);

    MassNeighborhoodExtended GetMassNeighborhoodExtended(int frameIdx, int strandIdx, int massIdx, Flt neighDist); // For input

    MassNeighborhoodExtended GetMassNeighborhoodExtended(int strandIdx, int massIdx, Flt neighDist); // For output

    Flt NeighborhoodMetricExtended(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices);

    Flt GetNearestInputNeighborExtended(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices);

    Flt NeighborhoodMetricExtendedNew(MassNeighborhoodExtended& neigh1, MassNeighborhoodExtended& neigh2, vector<int>& matchIndices);

    Flt GetNearestInputNeighborExtendedNew(MassNeighborhoodExtended& neighOut, int& nearesNeighIdx, vector<int>& matchIndices);

    vector<SegmentPatch> GetSegmentPatch(CMassSpringData& msData);

    Flt SegmentPatchCompatibility(SegmentPatch* ptrPatch1, SegmentPatch* ptrPatch2);

    vector<Vec3f> GetConnectedInputPatch(int patchLength);

    void UpdateCoeffMatDiagonalVal(int idx, Flt wt, Vec3f pos, bool flagUpdateMat = true);

    void UpdateCoeffMatPairVals(int idxi, int idxj, Flt wt, Vec3f pr, bool flagUpdateMat = true);

    vector<CMassSpringData> LoadStrandsFromTXT(const string& fileName);

    vector<CMassSpringData> InitializeStrandsViaPoissonDisk(const string& fileName);

    int GetNumOfUnknownsTotal(vector<CMassSpringData>& strands);

    void SetSequencesVelocity(vector<MassSpringSequence>& vecSequence, bool flagTemporallyToroidal = false);

    inline Flt VecPrDifference(vector<Vec3f>& vecPr1, vector<Vec3f>& vecPr2)
    {
        Flt diff = 0.0f;
        int sizeDiff = int(vecPr2.size() - vecPr1.size());
        for (int i = 0; i < int(vecPr1.size()); i++)
        {
            diff += dist2(vecPr1[i], vecPr2[i + sizeDiff]) * exp(CMassSpringSynConfig::m_gaussianSigma * i * i);
        }
        return diff;
    }

    inline Flt VecPrDifferenceTemporal(vector<Vec3f>& vecPr1, vector<Vec3f>& vecPr2)
    {
        if (vecPr2.size() < vecPr1.size()) return 1e10;
        Flt diff = 0.0f;
        int sizeDiff = int(vecPr2.size() - vecPr1.size());
        for (int i = 0; i < int(vecPr1.size()); i++)
        {
            Flt i2 = (i + 1) * (i + 1);
            diff += dist2(vecPr1[i], vecPr2[i]) * exp(CMassSpringSynConfig::m_temporalSigma * i2);
        }
        return diff;
    }

    int m_stepCount;
    int m_serialCount;
    CMassSpringSynConfig* m_ptrSynConfig;
    CCrossList* m_ptrCoeffMatrix;
    Eigen::VectorXf m_vecCx;
    Eigen::VectorXf m_vecCy;
    Eigen::VectorXf m_vecCz;
    vector<MassNeighborhoodExtended> m_vecInputNeighExtended;
    vector<SegmentPatch> m_vecInputSegmentPatch;
    vector<MassSpringSequence> m_vecInputSequence;
    vector<CMassSpringData> m_vecOutputStrand;
    vector<MassSpringSequence> m_vecOutputSequence;
};

#endif // MASSSPRINGSYN_H
