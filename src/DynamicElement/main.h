
#ifndef MAIN_H
#define MAIN_H

#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glut32.lib")
#pragma warning( disable: 4244 )	// disable warning of data conversion.
#pragma warning( disable: 4996 )	// disable warning of security.

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::ostream;
using std::istream;
using std::ifstream;
using std::ofstream;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;
using std::fill;
using std::copy;
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
using std::sqrt;
using std::abs;
using std::exp;
using std::cos;
using std::sin;
using std::atan;
using std::atan2;
using std::log;
using std::exp;
using std::ceil;
using std::floor;
#include <ctime>
using std::clock;
using std::time;
using std::difftime;
#include <cstdlib>
#include <cassert>
#include <fstream>
using std::ostringstream;
using std::ios;
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <GL/glut.h>
#include <atlimage.h> // To use CImage class

//------------------------------CONSTANTS---------------------------------------------
// User-specified data type
typedef float Flt;
typedef unsigned int uint;

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

const float HALF_PI = 2.0 * atan(1.0);
const float M_PI = 4.0 * atan(1.0);
const float TWO_PI = 8.0 * atan(1.0);
const float INV_PI = 1.0/M_PI;
const float INV_TWO_PI = 1.0/TWO_PI;
const float PI_INV180 = M_PI   / 180.;
const float INV_PI180 = INV_PI * 180;
const float ONE_THIRD = 1./3.;
const float INV_255 = 1. / 255;
const float SQRT_TWO = sqrt(2.);
const float	INV_SQRT_TWO = 1. / sqrt(2.);

#define EPS					1e-5
#define IS_ZERO(x)			(x < EPS && x > -EPS)

#define DELETE_OBJECT(p)	{if(NULL!=p){delete p;p = NULL;}}
#define DELETE_ARRAY(p)		{if(NULL!=p){delete[] p;p = NULL;}}

inline float Radians(const float &deg) {return PI_INV180 * deg; }
inline float Degrees(const float &rad) {return INV_PI180 * rad; }

#endif MAIN_H
