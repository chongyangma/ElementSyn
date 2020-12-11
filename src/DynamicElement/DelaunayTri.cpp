#include "DelaunayTri.h"

CDelaunayTri::CDelaunayTri()
{
    m_pMesh = NULL;
}

CDelaunayTri::~CDelaunayTri()
{
    if (m_pMesh != NULL)
    {
        TRIANGLE_PTR triPtr = m_pMesh->pTriArr;
        if (triPtr != NULL)
        {
            delete[] triPtr;
        }
        VERTEX2D_PTR verPtr = m_pMesh->pVerArr;
        if (verPtr != NULL)
        {
            delete[] verPtr;
        }
        delete m_pMesh;
        m_pMesh = NULL;
    }
}

void CDelaunayTri::InitMesh(vector<Vec2f> vecPoint)
{
    int ver_num = int(vecPoint.size());

    m_pMesh = new MESH;
    m_pMesh->pTriArr = NULL;
    m_pMesh->pVerArr = new VERTEX2D[ver_num + 3];
    m_pMesh->vertex_num = ver_num;

    VERTEX2D_PTR pVertices = m_pMesh->pVerArr;
    for (int i = 3; i < ver_num + 3; i++)
    {
        pVertices[i].x = vecPoint[i - 3][0];
        pVertices[i].y = vecPoint[i - 3][1];
    }
}

REAL CDelaunayTri::CounterClockWise(VERTEX2D_PTR pa, VERTEX2D_PTR pb, VERTEX2D_PTR pc)
{
    return ((pb->x - pa->x) * (pc->y - pb->y) - (pc->x - pb->x) * (pb->y - pa->y));
}

REAL CDelaunayTri::CounterClockWise(VERTEX2D va, VERTEX2D vb, VERTEX2D vc)
{
    return ((vb.x - va.x) * (vc.y - vb.y) - (vc.x - vb.x) * (vb.y - va.y));
}

TRIANGLE_PTR CDelaunayTri::AddTriangleNode(TRIANGLE_PTR pPrevTri, int i1, int i2, int i3)
{
    // test if 3 vertices are co-linear
    if (CounterClockWise(m_pMesh->pVerArr[i1], m_pMesh->pVerArr[i2], m_pMesh->pVerArr[i3]) == 0)
    {
        return NULL;
    }

    TRIANGLE_PTR pNewTestTri = new TRIANGLE;

    pNewTestTri->i1 = i1;
    pNewTestTri->i2 = i2;
    pNewTestTri->i3 = i3;

    // insert after prev triangle
    if (pPrevTri == NULL) // add root
    {
        m_pMesh->pTriArr = pNewTestTri;
        pNewTestTri->pNext = NULL;
        pNewTestTri->pPrev = NULL;
    }
    else
    {
        pNewTestTri->pNext = pPrevTri->pNext;
        pNewTestTri->pPrev = pPrevTri;
        if (pPrevTri->pNext != NULL)
        {
            pPrevTri->pNext->pPrev = pNewTestTri;
        }
        pPrevTri->pNext = pNewTestTri;
    }

    return pNewTestTri;
}

void CDelaunayTri::AddBoundingBox()
{
    REAL max_c = 0;
    REAL max_x = 0;
    REAL max_y = 0;
    REAL t;

    for (int i = 3; i < m_pMesh->vertex_num + 3; i++)
    {
        t = abs(m_pMesh->pVerArr[i].x);
        if (max_x < t)
        {
            max_x = t;
        }

        t = abs(m_pMesh->pVerArr[i].y);
        if (max_y < t)
        {
            max_y = t;
        }
    }

    max_c = max_x > max_y ? max_x : max_y;

    VERTEX2D v1 = { 0, 4 * max_c };
    VERTEX2D v2 = { -4 * max_c, -4 * max_c };
    VERTEX2D v3 = { 4 * max_c, 0 };

    // Assign to Vertex array
    *(m_pMesh->pVerArr) = v1;
    *(m_pMesh->pVerArr + 1) = v2;
    *(m_pMesh->pVerArr + 2) = v3;

    // add the Triangle bounding-box
    AddTriangleNode(NULL, 0, 1, 2);
}

void CDelaunayTri::RemoveBoundingBox()
{
    int statify[3] = { 0, 0, 0 };
    int vertex_index;
    int* pi;
    int k = 1;

    // Remove the first triangle bounding-box
    //pMesh->pTriArr = pMesh->pTriArr->pNext;
    //pMesh->pTriArr->pPrev = NULL; // as head

    TRIANGLE_PTR pTri = m_pMesh->pTriArr;
    TRIANGLE_PTR pNext = NULL;
    while (pTri != NULL)
    {
        pNext = pTri->pNext;

        statify[0] = 0;
        statify[1] = 0;
        statify[2] = 0;

        pi = &(pTri->i1);
        for (int j = 0, k = 1; j < 3; j++, k *= 2)
        {
            vertex_index = *pi++;

            if (vertex_index == 0 || vertex_index == 1 || vertex_index == 2) // bounding box vertex
            {
                statify[j] = k;
            }
        }

        switch (statify[0] | statify[1] | statify[2])
        {
        case 0: // no statify
            break;
        case 1:
        case 2:
        case 4: // 1 statify, remove 1 triangle, 1 vertex
            RemoveTriangleNode(pTri);
            break;
        case 3:
        case 5:
        case 6: // 2 statify, remove 1 triangle, 2 vertices
            RemoveTriangleNode(pTri);
            break;
        case 7: // 3 statify, remove 1 triangle, 3 vertices
            RemoveTriangleNode(pTri);
            break;
        default:
            break;
        }

        // go to next item
        pTri = pNext;
    }
}

void CDelaunayTri::RemoveTriangleNode(TRIANGLE_PTR pTri)
{
    if (pTri == NULL)
    {
        return;
    }

    // remove from the triangle list
    if (pTri->pPrev != NULL)
    {
        pTri->pPrev->pNext = pTri->pNext;
    }
    else // remove the head, need to reset the root node
    {
        m_pMesh->pTriArr = pTri->pNext;
    }

    if (pTri->pNext != NULL)
    {
        pTri->pNext->pPrev = pTri->pPrev;
    }

    // deallocate memory
    delete pTri;
}

REAL CDelaunayTri::InCircle(VERTEX2D_PTR pa, VERTEX2D_PTR pb, VERTEX2D_PTR pp, VERTEX2D_PTR pd)
{
    REAL det;
    REAL alift, blift, plift, bdxpdy, pdxbdy, pdxady, adxpdy, adxbdy, bdxady;

    REAL adx = pa->x - pd->x;
    REAL ady = pa->y - pd->y;

    REAL bdx = pb->x - pd->x;
    REAL bdy = pb->y - pd->y;

    REAL pdx = pp->x - pd->x;
    REAL pdy = pp->y - pd->y;

    bdxpdy = bdx * pdy;
    pdxbdy = pdx * bdy;
    alift = adx * adx + ady * ady;

    pdxady = pdx * ady;
    adxpdy = adx * pdy;
    blift = bdx * bdx + bdy * bdy;

    adxbdy = adx * bdy;
    bdxady = bdx * ady;
    plift = pdx * pdx + pdy * pdy;

    det = alift * (bdxpdy - pdxbdy) + blift * (pdxady - adxpdy) + plift * (adxbdy - bdxady);

    return -det;
}

REAL CDelaunayTri::InTriangle(VERTEX2D_PTR pVer, TRIANGLE_PTR pTri)
{
    int vertex_index;
    VERTEX2D_PTR pV1, pV2, pV3;

    vertex_index = pTri->i1;
    pV1 = (VERTEX2D_PTR)(m_pMesh->pVerArr + vertex_index);
    vertex_index = pTri->i2;
    pV2 = (VERTEX2D_PTR)(m_pMesh->pVerArr + vertex_index);
    vertex_index = pTri->i3;
    pV3 = (VERTEX2D_PTR)(m_pMesh->pVerArr + vertex_index);

    REAL ccw1 = CounterClockWise(pV1, pV2, pVer);
    REAL ccw2 = CounterClockWise(pV2, pV3, pVer);
    REAL ccw3 = CounterClockWise(pV3, pV1, pVer);

    REAL r = -1;
    if (ccw1 > 0 && ccw2 > 0 && ccw3 > 0)
    {
        r = 1;
    }
    else if (ccw1 * ccw2 * ccw3 == 0 && (ccw1 * ccw2 > 0 || ccw1 * ccw3 > 0 || ccw2 * ccw3 > 0))
    {
        r = 0;
    }

    return r;
}

bool CDelaunayTri::FlipTest(TRIANGLE_PTR pTestTri)
{
    bool flipped = false;

    int index_a = pTestTri->i1;
    int index_b = pTestTri->i2;
    int index_p = pTestTri->i3;

    int statify[3] = { 0, 0, 0 };
    int vertex_index;
    int* pi;
    int k = 1;

    // find the triangle which has edge consists of start and end
    TRIANGLE_PTR pTri = m_pMesh->pTriArr;

    int index_d = -1;
    while (pTri != NULL)
    {
        statify[0] = 0;
        statify[1] = 0;
        statify[2] = 0;

        pi = &(pTri->i1);
        for (int j = 0, k = 1; j < 3; j++, k *= 2)
        {
            vertex_index = *pi++;
            if (vertex_index == index_a || vertex_index == index_b)
            {
                statify[j] = k;
            }
        }

        switch (statify[0] | statify[1] | statify[2])
        {
        case 3:
            if (CounterClockWise((VERTEX2D_PTR)(m_pMesh->pVerArr + index_a), (VERTEX2D_PTR)(m_pMesh->pVerArr + index_b), (VERTEX2D_PTR)(m_pMesh->pVerArr + pTri->i3)) < 0)
            {
                index_d = pTri->i3;
            }
            break;
        case 5:
            if (CounterClockWise((VERTEX2D_PTR)(m_pMesh->pVerArr + index_a), (VERTEX2D_PTR)(m_pMesh->pVerArr + index_b), (VERTEX2D_PTR)(m_pMesh->pVerArr + pTri->i2)) < 0)
            {
                index_d = pTri->i2;
            }
            break;
        case 6:
            if (CounterClockWise((VERTEX2D_PTR)(m_pMesh->pVerArr + index_a), (VERTEX2D_PTR)(m_pMesh->pVerArr + index_b), (VERTEX2D_PTR)(m_pMesh->pVerArr + pTri->i1)) < 0)
            {
                index_d = pTri->i1;
            }
            break;
        default:
            break;
        }

        if (index_d != -1)
        {
            VERTEX2D_PTR pa = (VERTEX2D_PTR)(m_pMesh->pVerArr + index_a);
            VERTEX2D_PTR pb = (VERTEX2D_PTR)(m_pMesh->pVerArr + index_b);
            VERTEX2D_PTR pd = (VERTEX2D_PTR)(m_pMesh->pVerArr + index_d);
            VERTEX2D_PTR pp = (VERTEX2D_PTR)(m_pMesh->pVerArr + index_p);

            if (InCircle(pa, pb, pp, pd) < 0) // not local Delaunay
            {
                flipped = true;

                // add new triangle adp,  dbp, remove abp, abd.
                // allocate memory for adp
                TRIANGLE_PTR pT1 = AddTriangleNode(pTestTri, pTestTri->i1, index_d, pTestTri->i3);
                // allocate memory for dbp
                TRIANGLE_PTR pT2 = AddTriangleNode(pT1, index_d, pTestTri->i2, index_p);
                // remove abp
                RemoveTriangleNode(pTestTri);
                // remove abd
                RemoveTriangleNode(pTri);

                FlipTest(pT1); // pNewTestTri satisfies CCW order
                FlipTest(pT2); // pNewTestTri2  satisfies CCW order

                break;
            }
        }

        // go to next item
        pTri = pTri->pNext;
    }

    return flipped;
}

void CDelaunayTri::InsertInTriangle(TRIANGLE_PTR pTargetTri, int ver_index)
{
    int index_a, index_b, index_c;
    TRIANGLE_PTR pTri = NULL;
    TRIANGLE_PTR pNewTri = NULL;

    pTri = pTargetTri;
    if (pTri == NULL)
    {
        return;
    }

    // Inset p into target triangle
    index_a = pTri->i1;
    index_b = pTri->i2;
    index_c = pTri->i3;

    // Insert edge pa, pb, pc
    for (int i = 0; i < 3; i++)
    {
        // allocate memory
        if (i == 0)
        {
            pNewTri = AddTriangleNode(pTri, index_a, index_b, ver_index);
        }
        else if (i == 1)
        {
            pNewTri = AddTriangleNode(pTri, index_b, index_c, ver_index);
        }
        else
        {
            pNewTri = AddTriangleNode(pTri, index_c, index_a, ver_index);
        }

        // go to next item
        if (pNewTri != NULL)
        {
            pTri = pNewTri;
        }
        else
        {
            pTri = pTri;
        }
    }

    // Get the three sub-triangles
    pTri = pTargetTri;
    TRIANGLE_PTR pTestTri[3];
    for (int i = 0; i < 3; i++)
    {
        pTestTri[i] = pTri->pNext;

        pTri = pTri->pNext;
    }

    // remove the Target Triangle
    RemoveTriangleNode(pTargetTri);

    for (int i = 0; i < 3; i++)
    {
        // Flip test
        FlipTest(pTestTri[i]);
    }
}

void CDelaunayTri::InsertOnEdge(TRIANGLE_PTR pTargetTri, int ver_index)
{
    int index_a, index_b, index_c;
    TRIANGLE_PTR pTri = NULL;
    TRIANGLE_PTR pNewTri = NULL;

    pTri = pTargetTri;
    if (pTri == NULL)
    {
        return;
    }

    // Inset p into target triangle
    index_a = pTri->i1;
    index_b = pTri->i2;
    index_c = pTri->i3;

    // Insert edge pa, pb, pc
    for (int i = 0; i < 3; i++)
    {
        // allocate memory
        if (i == 0)
        {
            pNewTri = AddTriangleNode(pTri, index_a, index_b, ver_index);
        }
        else if (i == 1)
        {
            pNewTri = AddTriangleNode(pTri, index_b, index_c, ver_index);
        }
        else
        {
            pNewTri = AddTriangleNode(pTri, index_c, index_a, ver_index);
        }

        // go to next item
        if (pNewTri != NULL)
        {
            pTri = pNewTri;
        }
        else
        {
            pTri = pTri;
        }
    }

    // Get the two sub-triangles
    pTri = pTargetTri;
    TRIANGLE_PTR pTestTri[2];
    for (int i = 0; i < 2; i++)
    {
        pTestTri[i] = pTri->pNext;
        pTri = pTri->pNext;
    }

    // remove the Target Triangle
    RemoveTriangleNode(pTargetTri);

    for (int i = 0; i < 2; i++)
    {
        // Flip test
        FlipTest(pTestTri[i]);
    }
}

void CDelaunayTri::Insert(int ver_index)
{
    VERTEX2D_PTR pVer = (VERTEX2D_PTR)(m_pMesh->pVerArr + ver_index);
    TRIANGLE_PTR pTargetTri = NULL;
    TRIANGLE_PTR pEqualTri1 = NULL;
    TRIANGLE_PTR pEqualTri2 = NULL;

    int j = 0;
    TRIANGLE_PTR pTri = m_pMesh->pTriArr;
    while (pTri != NULL)
    {
        REAL r = InTriangle(pVer, pTri);
        if (r > 0) // should be in triangle
        {
            pTargetTri = pTri;
        }
        else if (r == 0) // should be on edge
        {
            if (j == 0)
            {
                pEqualTri1 = pTri;
                j++;
            }
            else
            {
                pEqualTri2 = pTri;
            }
        }

        pTri = pTri->pNext;
    }

    if (pEqualTri1 != NULL && pEqualTri2 != NULL)
    {
        InsertOnEdge(pEqualTri1, ver_index);
        InsertOnEdge(pEqualTri2, ver_index);
    }
    else
    {
        InsertInTriangle(pTargetTri, ver_index);
    }
}

void CDelaunayTri::IncrementalDelaunay()
{
    // Add a appropriate triangle bounding-box to contain V
    AddBoundingBox();

    // Get a vertex/point vi from V and Insert(vi)
    for (int i = 3; i < m_pMesh->vertex_num + 3; i++)
    {
        Insert(i);
    }

    // Remove the bounding box
    RemoveBoundingBox();
}
