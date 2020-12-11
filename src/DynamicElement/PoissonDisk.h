#ifndef POISSONDISK_H
#define POISSONDISK_H

#include "main.h"
#include "vec.h"

template <unsigned int N, class T>
class CPoissonDisk
{
public:
    CPoissonDisk<N, T>();

    vector<Vec<N, T>> GeneratePoints(Vec<N, T> pMin, Vec<N, T> pMax, T distSqrMin);

private:
};

template <unsigned int N, class T>
CPoissonDisk<N, T>::CPoissonDisk()
{
}

template <unsigned int N, class T>
vector<Vec<N, T>> CPoissonDisk<N, T>::GeneratePoints(Vec<N, T> pMin, Vec<N, T> pMax, T distSqrMin)
{
    vector<Vec<N, T>> vecPoint;
    int unsuccessfulTimes = 0;
    Vec<N, T> pMean = (pMin + pMax) * 0.5;
    Vec<N, T> pd = pMax - pMin;
    while (unsuccessfulTimes < 2000)
    {
        Vec<N, T> p;
        for (unsigned int i = 0; i < N; i++)
        {
            p[i] = (rand() / T(RAND_MAX) - 0.5) * pd[i] + pMean[i];
        }
        bool flag = false;
        for (int i = 0; i < int(vecPoint.size()); i++)
        {
            Flt distSqr = dist2(p, vecPoint[i]);
            if (distSqr < distSqrMin)
            {
                flag = true;
                break;
            }
        }
        if (flag == false)
        {
            unsuccessfulTimes = 0;
            vecPoint.push_back(p);
        }
        else
        {
            unsuccessfulTimes++;
        }
    }

    return vecPoint;
}

typedef CPoissonDisk<2, Flt> PoissonDisk2f;
typedef CPoissonDisk<3, Flt> PoissonDisk3f;

#endif // POISSONDISK_H
