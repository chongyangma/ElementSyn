
#ifndef MKLSOLVER_H
#define MKLSOLVER_H

#include "MathTypes.h"
#include "mkl_dss.h"
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#pragma comment(lib,"mkl_intel_c_dll.lib")
#pragma comment(lib,"mkl_solver.lib")
#pragma comment(lib,"mkl_intel_thread_dll.lib")

class CMKLsolver
{
public:
	CMKLsolver();

	~CMKLsolver();

	void AllocatePtrs(int nRows, int nCols, int nNonZeros);

	void DeallocatePtrs();

	void SetCoeffsFromCoeffMatrix(CDenseMatrix* ptrMatD);

	void SetCoeffsFromCoeffMatrix(CCrossList* ptrMatD);

	void SetCoeffsFromCSparseMatCRS(CSparseMatrixCRS& matS);

	vector<Flt> GetSolutionViaDCG(vector<Flt>& vecC, vector<Flt>& initGuess);

private:
	int m_nRows, m_nCols, m_nNonZeros;
	int m_nRhs;
	int* m_ptrRowIndex;
	int* m_ptrColumns;
	_DOUBLE_PRECISION_t* m_ptrValues;
	_DOUBLE_PRECISION_t* m_ptrSolValues;
	_DOUBLE_PRECISION_t* m_ptrRhs;
	_MKL_DSS_HANDLE_t m_handle;
};

#endif MKLSOLVER_H
