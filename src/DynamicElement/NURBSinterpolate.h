#ifndef NURBSINTERPOLATE_H
#define NURBSINTERPOLATE_H

#include <iostream>
#include <vector>

using namespace std;

#define ORDER 4

// struct definition of triple diagonal matrix
typedef struct CoeffMatrix
{
    vector<float> vecA; // Lower minor secondary elements
    vector<float> vecB; // Diagonal elements
    vector<float> vecC; // Upper minor secondary elements
} CoeffMatrix;

class CNURBSinterpolate
{
public:
    vector<float> InterpolateUniformOpen(const vector<float>& vecShapePts);

    vector<float> InterpolateUniformOpen(const vector<float>& vecShapePts, int interval);

    vector<float> InterpolateUniformOpen(const vector<float>& vecShapePts, vector<float>& vecU);

    void SetMatrixUniformOpen(int numOfCtrlPts);

    float NURBSvalue(float u, const vector<float>& vecCtrlPts);

    float BsplineBasis(int i, int k, float u);

private:
    CoeffMatrix m_MatrixA;
    vector<float> m_vecKnot;
};

#endif // NURBSINTERPOLATE_H
