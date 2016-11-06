#ifndef _TAUCS_ADDON_
#define _TAUCS_ADDON_

/****************************************************************/
// Sivan's library!!!
#ifndef WIN32
#include <unistd.h>
#include <pthread.h>
#endif

#include <vector>
#include <map>


#define TAUCS_CORE_DOUBLE

extern "C" {
#include "taucs.h"
}


//////////////////////////////////////////////////////////////////////////
typedef  double  TaucsType;
/************************************************************************/


class TaucsAddOn
{
public:

	// Assuming nothing about the result (the result is NOT stored symmetric).
	static taucs_ccs_matrix* Mul2NonSymmetricMatrices(
		const taucs_ccs_matrix* matA,
		const taucs_ccs_matrix* matB);

	// For usage when it's known that the result is symmetric, like A^T * A.
	static taucs_ccs_matrix* Mul2NonSymmMatSymmResult(
		const taucs_ccs_matrix* matA,
		const taucs_ccs_matrix* matB);

	// Computes the transpose of a matrix.
	static taucs_ccs_matrix* MatrixTranspose(const taucs_ccs_matrix* mat);

	static taucs_ccs_matrix* CreateTaucsMatrixFromColumns(
		const std::vector< std::map<int, TaucsType> >& cols, 
		int nRows,
		int flags);

	// Multiplies matA by x and stores the result in b. Assumes all memory has 
	// been allocated and the sizes match; assumes matA is not symmetric!!
	static void MulNonSymmMatrixVector(
		const taucs_ccs_matrix* matA,
		const TaucsType* x,
		TaucsType* b);

	// Adds two vectors vecA and VecB and stores the result in vecResult.
	// Assumes all memory has been allocated and the sizes match!
	static void Add2Vectors(
		int n, 
		const TaucsType* vecA, 
		const TaucsType* vecB, 
		TaucsType* vecResult);

	// Copy mat to a new matrix, memory will be allocated during the copy process.
	static taucs_ccs_matrix* MatrixCopy(const taucs_ccs_matrix* mat);

};
#endif // _TAUCS_ADDON_