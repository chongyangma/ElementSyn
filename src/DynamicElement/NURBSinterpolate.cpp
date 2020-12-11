#include "NURBSinterpolate.h"

vector<float> Chase1(CoeffMatrix matA, vector<float> vecD)
{
    int n = int(vecD.size());
    vector<float> vecP(n), vecQ(n - 1), vecZ(n), vecX(n);
    vecP[0] = matA.vecB[0];
    for (int i = 0; i < n - 1; i++)
    {
        vecQ[i] = matA.vecC[i] / vecP[i];
        vecP[i + 1] = matA.vecB[i + 1] - matA.vecA[i] * vecQ[i];
    }
    vecZ[0] = vecD[0] / vecP[0];
    for (int i = 1; i < n; i++)
    {
        vecZ[i] = (vecD[i] - matA.vecA[i - 1] * vecZ[i - 1]) / vecP[i];
    }
    vecX[n - 1] = vecZ[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        vecX[i] = vecZ[i] - vecQ[i] * vecX[i + 1];
    }
    return vecX;
}

vector<float> CNURBSinterpolate::InterpolateUniformOpen(const vector<float>& vecShapePts)
{
    int numOfShapePts = int(vecShapePts.size());
    int numOfCtrlPts = numOfShapePts + 2;
    int order = numOfCtrlPts + ORDER;
    m_vecKnot.resize(order);
    for (int i = 0; i < order; i++)
    {
        m_vecKnot[i] = float(i);
    }
    SetMatrixUniformOpen(numOfCtrlPts);
    vector<float> vecB;
    vecB.push_back(0.0f);
    for (int i = 0; i < int(vecShapePts.size()); i++)
    {
        vecB.push_back(6.0f * vecShapePts[i]);
    }
    vecB.push_back(0.0f);
    vector<float> vecCtrlPts = Chase1(m_MatrixA, vecB);
    return vecCtrlPts;
}

vector<float> CNURBSinterpolate::InterpolateUniformOpen(const vector<float>& vecShapePts, int interval)
{
    int sz = int(vecShapePts.size() - 1) * interval + 1;
    vector<float> vecU(sz);
    for (int i = 0; i < sz; i++)
    {
        vecU[i] = 3.0f + i / float(interval);
    }
    vector<float> vecInterpPts = InterpolateUniformOpen(vecShapePts, vecU);
    return vecInterpPts;
}

vector<float> CNURBSinterpolate::InterpolateUniformOpen(const vector<float>& vecShapePts, vector<float>& vecU)
{
    vector<float> vecCtrlPts = InterpolateUniformOpen(vecShapePts);
    vector<float> vecInterpPts(vecU.size());
    for (int i = 0; i < int(vecU.size()); i++)
    {
        float u = vecU[i];
        vecInterpPts[i] = NURBSvalue(u, vecCtrlPts);
    }
    return vecInterpPts;
}

void CNURBSinterpolate::SetMatrixUniformOpen(int numOfCtrlPts)
{
    vector<float> CoeffA((numOfCtrlPts - 1), 1.0f);
    vector<float> CoeffB((numOfCtrlPts - 0), 4.0f);
    vector<float> CoeffC((numOfCtrlPts - 1), 1.0f);

    CoeffB[0] = 6;
    CoeffC[0] = -6;
    CoeffA[numOfCtrlPts - 2] = -6;
    CoeffB[numOfCtrlPts - 1] = 6;

    m_MatrixA.vecA = CoeffA;
    m_MatrixA.vecB = CoeffB;
    m_MatrixA.vecC = CoeffC;
}

float CNURBSinterpolate::NURBSvalue(float u, const vector<float>& vecCtrlPts)
{
    float sum1 = 0.0f;
    float sum2 = 0.0f;
    for (int i = 0; i < int(vecCtrlPts.size()); i++)
    {
        float b = BsplineBasis(i, ORDER - 1, u);
        sum1 += b * vecCtrlPts[i];
        sum2 += b;
    }
    return (sum1 / sum2);
}

float CNURBSinterpolate::BsplineBasis(int i, int k, float u)
{
    // Implement according to \url{http://web.cs.wpi.edu/~matt/courses/cs563/talks/nurbs.html}
    if (k == 0)
    {
        if (u >= m_vecKnot[i] && u < m_vecKnot[i + 1])
        {
            return 1.0f;
        }
        else
        {
            return 0.0f;
        }
    }
    else
    {
        float w1 = (u - m_vecKnot[i]) / (m_vecKnot[i + k] - m_vecKnot[i]);
        float w2 = (m_vecKnot[i + k + 1] - u) / (m_vecKnot[i + k + 1] - m_vecKnot[i + 1]);
        return (w1 * BsplineBasis(i, k - 1, u) + w2 * BsplineBasis(i + 1, k - 1, u));
    }
}
