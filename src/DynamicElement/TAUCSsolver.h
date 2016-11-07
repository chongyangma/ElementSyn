
#ifndef TAUCSSOLVER_H
#define TAUCSSOLVER_H

/********************************************************************
created:	 9:11:2009   9:04
filename: 	 TAUCSsolver.h
author:		 Liangliang Nan
contact:     liangliang.nan@gmail.com
purpose:	 Wrapper of Taucs Solver for easy use 

modified by Chongyang Ma, 07/30/2012
*********************************************************************/

#include "taucs_matrix.h"
#include "MathTypes.h"
#include <string>

#pragma comment(lib, "libatlas.lib")
#pragma comment(lib, "libcblas.lib")
#pragma comment(lib, "libf77blas.lib")
#pragma comment(lib, "liblapack.lib")
#pragma comment(lib, "libmetis.lib")
#pragma comment(lib, "libtaucs.lib")
//#pragma comment(lib, "libtstatlas.lib")
#pragma comment(lib, "vcf2c.lib")

class CTAUCSsolver
{
public:
	std::string title() { return "[CTAUCSsolver]: "; }

	// solve for "A*x=b"
	// A: the coefficient symmetry matrix, 
	// b: the right side column vector
	// x: the result
	bool solve_symmetry(const Taucs_matrix<double>& A, const std::vector<double>& b, std::vector<double>& x);

	// solve for "A*x=b"
	// A: the coefficient asymmetry matrix, 
	// b: the right side column vector
	// x: the result
	bool solve_unsymmetry(const Taucs_matrix<double>& A, const std::vector<double>& b, std::vector<double>& x);

	// solve for "A*x=b" in least square sence
	// A: the coefficient m * n matrix (m >= n)
	// b: the right side column vector
	// x: the result
	bool solve_linear_least_square(const Taucs_matrix<double>& A, const std::vector<double>& b, std::vector<double>& x);

	vector<Flt> GetSolution(CDenseMatrix* ptrCoeffMatrix, vector<Flt>& vecB);

	vector<Flt> GetSolution(CCrossList* ptrCoeffMatrix, vector<Flt>& vecB);
};

#endif TAUCSSOLVER_H
