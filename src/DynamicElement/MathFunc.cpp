#include "MathFunc.h"

namespace machy_math
{
Vec4f QuatConjugate(const Vec4f& q0)
{
    Vec4f q;
    q[0] = q0[0];
    q[1] = -q0[1];
    q[2] = -q0[2];
    q[3] = -q0[3];

    return q;
}

Vec3f QuatMultiplyVec(const Vec4f& q0, const Vec3f& v0)
{
    // From Wiki (Quaternions and spatial rotation)
    Vec3f v;

    Flt t2 = q0[0] * q0[1];
    Flt t3 = q0[0] * q0[2];
    Flt t4 = q0[0] * q0[3];
    Flt t5 = -q0[1] * q0[1];
    Flt t6 = q0[1] * q0[2];
    Flt t7 = q0[1] * q0[3];
    Flt t8 = -q0[2] * q0[2];
    Flt t9 = q0[2] * q0[3];
    Flt t10 = -q0[3] * q0[3];

    v[0] = 2.f * ((t8 + t10) * v0[0] + (t6 - t4) * v0[1] + (t3 + t7) * v0[2]) + v0[0];
    v[1] = 2.f * ((t4 + t6) * v0[0] + (t5 + t10) * v0[1] + (t9 - t2) * v0[2]) + v0[1];
    v[2] = 2.f * ((t7 - t3) * v0[0] + (t2 + t9) * v0[1] + (t5 + t8) * v0[2]) + v0[2];

    return v;
}

Vec4f QuatSlerp(Vec4f q1, Vec4f q2, Flt wt)
{
    Flt dotVal = dot(q1, q2);
    if (dotVal < 0.0f) q2 = -q2; // To avoid the longer path
    Flt dotabs = abs(dotVal);
    if (dotabs >= 1.0f)
    {
        return q1;
    }
    Flt theta = acos(dotabs);
    Flt sv = sin(theta);
    Flt sv1 = sin(wt * theta);
    Flt sv2 = sin(theta - wt * theta);
    Vec4f q = (sv1 / sv) * q1 + (sv2 / sv) * q2;
    normalize(q);
    return q;
}

void SetMatrixFromQuat(Flt* ptrMat, Vec4f quat)
{
    for (int i = 0; i < 16; i++)
    {
        ptrMat[i] = 0.f;
    }
    ptrMat[15] = 1.f;
    normalize(quat);
    Flt s = quat[0];
    Flt vx = quat[1];
    Flt vy = quat[2];
    Flt vz = quat[3];
    ptrMat[0] = 1. - 2. * vy * vy - 2. * vz * vz;
    ptrMat[4] = 2. * vx * vy - 2. * s * vz;
    ptrMat[8] = 2. * vx * vz + 2. * s * vy;
    ptrMat[1] = 2. * vx * vy + 2. * s * vz;
    ptrMat[5] = 1. - 2. * vx * vx - 2. * vz * vz;
    ptrMat[9] = 2. * vy * vz - 2. * s * vx;
    ptrMat[2] = 2. * vx * vz - 2. * s * vy;
    ptrMat[6] = 2. * vy * vz + 2. * s * vx;
    ptrMat[10] = 1. - 2. * vx * vx - 2. * vy * vy;
}

void SetQuatFromMatrix(Flt* ptrMat, Vec4f& quat)
{
    Flt tr, s;
    tr = ptrMat[0] + ptrMat[5] + ptrMat[10];

    if (tr >= 0)
    {
        s = sqrt(tr + Flt(1));
        quat[0] = 0.5 * s;
        s = 0.5 / s;
        quat[1] = (ptrMat[6] - ptrMat[9]) * s;
        quat[2] = (ptrMat[8] - ptrMat[2]) * s;
        quat[3] = (ptrMat[1] - ptrMat[4]) * s;
    }
    else
    {
        int i = 0;
        if (ptrMat[5] > ptrMat[0])
        {
            i = 1;
        }
        if (ptrMat[10] > ptrMat[5 * i])
        {
            i = 2;
        }

        switch (i)
        {
        case 0:
            s = sqrt((ptrMat[0] - (ptrMat[5] + ptrMat[10])) + Flt(1));
            quat[1] = 0.5 * s;
            s = 0.5 / s;
            quat[2] = (ptrMat[4] + ptrMat[1]) * s;
            quat[3] = (ptrMat[2] + ptrMat[8]) * s;
            quat[0] = (ptrMat[6] - ptrMat[9]) * s;
            break;
        case 1:
            s = sqrt((ptrMat[5] - (ptrMat[10] + ptrMat[0])) + Flt(1));
            quat[2] = 0.5 * s;
            s = 0.5 / s;
            quat[3] = (ptrMat[9] + ptrMat[6]) * s;
            quat[1] = (ptrMat[4] + ptrMat[1]) * s;
            quat[0] = (ptrMat[8] - ptrMat[2]) * s;
            break;
        case 2:
            s = sqrt((ptrMat[10] - (ptrMat[0] + ptrMat[5])) + Flt(1));
            quat[3] = 0.5 * s;
            s = 0.5 / s;
            quat[1] = (ptrMat[2] + ptrMat[8]) * s;
            quat[2] = (ptrMat[9] + ptrMat[6]) * s;
            quat[0] = (ptrMat[1] - ptrMat[4]) * s;
            break;
        }
    }
    normalize(quat);
}

Vec3f GetTriFaceNormal(const Vec3f& v1, const Vec3f& v2, const Vec3f& v3)
{
    Vec3f e1 = v2 - v1;
    Vec3f e2 = v3 - v1;
    Vec3f n = cross(e1, e2);
    normalize(n);
    return n;
}

Flt GetVertexTriDist(const Vec3f& v1, const Vec3f& v2, const Vec3f& v3, const Vec3f& p0, const Vec3f& direc)
{
    Vec3f n = GetTriFaceNormal(v1, v2, v3);
    Flt vd = dot(n, direc);
    if (vd >= 0.f)
    {
        //return 1.f; // Do NOT return for double-sided mesh!
    }
    Flt dCoeff = -dot(n, v1);
    Flt t = -(dot(n, p0) + dCoeff) / vd;
    if (t > 0.f)
    {
        //return 1.f;
    }
    return t;
}

Flt GetVertexTriDistNew(const Vec3f& v1, const Vec3f& v2, const Vec3f& v3, const Vec3f& p0, const Vec3f& direc)
{
    Vec3f d1 = v1 - p0;
    Vec3f d2 = v2 - p0;
    Vec3f d3 = v3 - p0;
    Flt t1 = dot(d1, direc);
    Flt t2 = dot(d2, direc);
    Flt t3 = dot(d3, direc);
    return min(t1, min(t2, t3));
}

Flt GetTriPairDist(const Vec3f& v11, const Vec3f& v12, const Vec3f& v13, const Vec3f& v21, const Vec3f& v22, const Vec3f& v23, const Vec3f& direc)
{
    Flt d0 = GetVertexTriDist(v21, v22, v23, v11, direc);
    Flt d1 = GetVertexTriDist(v21, v22, v23, v12, direc);
    Flt d2 = GetVertexTriDist(v21, v22, v23, v13, direc);
    return min(min(d0, d1), d2);
}

Flt GetTriPairDistNew(const Vec3f& v11, const Vec3f& v12, const Vec3f& v13, const Vec3f& v21, const Vec3f& v22, const Vec3f& v23, const Vec3f& direc)
{
    Flt d0 = GetVertexTriDistNew(v21, v22, v23, v11, direc);
    Flt d1 = GetVertexTriDistNew(v21, v22, v23, v12, direc);
    Flt d2 = GetVertexTriDistNew(v21, v22, v23, v13, direc);
    return min(min(d0, d1), d2);
}

vector<Flt> GetSolution(CCrossList* ptrCoeffMatrix, const vector<Flt>& vecB)
{
    typedef Eigen::SparseMatrix<float> SparseMatrix;
    typedef Eigen::Triplet<float> Triplet;

    int nCols = ptrCoeffMatrix->GetColNum();
    int nRows = ptrCoeffMatrix->GetRowNum();

    Eigen::SparseMatrix<float> M(nRows, nCols);
    std::vector<Triplet> pending;

    OLink* ptrRhead = ptrCoeffMatrix->GetPtrRhead();
    for (int i = 1; i <= nRows; i++)
    {
        for (OLNode* q = ptrRhead[i]; q != NULL; q = q->m_ptrR)
        {
            pending.push_back(Triplet(q->m_x - 1, q->m_y - 1, q->m_value));
        }
    }
    M.setFromTriplets(pending.begin(), pending.end());
    M.makeCompressed();

    Eigen::VectorXf b(nRows);
    for (int i = 0; i < nRows; i++)
    {
        b[i] = vecB[i];
    }
    Eigen::VectorXf x(nRows);
    vector<Flt> vecX(nRows);

    Eigen::ConjugateGradient<SparseMatrix> cg;
    cg.compute(M);
    x = cg.solve(b);

    for (int i = 0; i < nRows; i++)
    {
        vecX[i] = x[i];
    }
    return vecX;
}

} // namespace machy_math
