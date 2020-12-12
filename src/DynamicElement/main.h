#ifndef MAIN_H
#define MAIN_H

#pragma warning(disable : 4244) // disable warning of data conversion.
#pragma warning(disable : 4305) // disable warning of data conversion.
#pragma warning(disable : 4996) // disable warning of security.

#include <string>
using std::string;
#include <iostream>
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <algorithm>
using std::copy;
using std::fill;
using std::max;
using std::min;
using std::sort;
using std::swap;
#include <vector>
using std::vector;
#include <list>
using std::list;
#include <set>
using std::set;
#include <map>
using std::map;
#include <deque>
using std::deque;
#include <cmath>
using std::abs;
using std::atan;
using std::atan2;
using std::ceil;
using std::cos;
using std::exp;
using std::floor;
using std::log;
using std::sin;
using std::sqrt;
#include <ctime>
using std::clock;
using std::difftime;
using std::time;
#include <cassert>
#include <cstdlib>
#include <fstream>
using std::ios;
using std::ostringstream;
#include <cstdio>
#include <iomanip>
#include <sstream>
#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif
#ifdef WIN32
    #include <windows.h>
#endif

//------------------------------CONSTANTS---------------------------------------------
// User-specified data type
typedef float Flt;
typedef unsigned int uint;

#ifndef INFINITY
    #define INFINITY HUGE_VAL
#endif

const float HALF_PI = 2.0 * atan(1.0);
//#ifdef WIN32
//const float M_PI = 4.0 * atan(1.0);
//#endif
const float TWO_PI = 8.0 * atan(1.0);
const float INV_PI = 1.0 / M_PI;
const float INV_TWO_PI = 1.0 / TWO_PI;
const float PI_INV180 = M_PI / 180.;
const float INV_PI180 = INV_PI * 180;
const float ONE_THIRD = 1. / 3.;
const float INV_255 = 1. / 255;
const float SQRT_TWO = sqrt(2.);
const float INV_SQRT_TWO = 1. / sqrt(2.);

#define EPS 1e-5
#define IS_ZERO(x) (x < EPS && x > -EPS)

#define DELETE_OBJECT(p) \
    {                    \
        if (NULL != p)   \
        {                \
            delete p;    \
            p = NULL;    \
        }                \
    }
#define DELETE_ARRAY(p) \
    {                   \
        if (NULL != p)  \
        {               \
            delete[] p; \
            p = NULL;   \
        }               \
    }

inline float Radians(const float& deg)
{
    return PI_INV180 * deg;
}
inline float Degrees(const float& rad)
{
    return INV_PI180 * rad;
}

#endif // MAIN_H
