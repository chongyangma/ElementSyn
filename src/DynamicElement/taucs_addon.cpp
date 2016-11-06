
#include <cassert>
#include "taucs_addon.h"

using namespace std;

/************************************************************************/


// Assuming nothing about the result (the result is NOT stored symmetric).
taucs_ccs_matrix* TaucsAddOn::Mul2NonSymmetricMatrices(const taucs_ccs_matrix* matA,
													   const taucs_ccs_matrix* matB) 
{
	// Compatibility of dimensions        
	if (matA->m != matB->n)
		return NULL;

	if ((matA->flags & TAUCS_SYMMETRIC) ||
		(matB->flags & TAUCS_LOWER) )
		return NULL;

	// (m x n)*(n x k) = (m x k)
	int m=matA->m;
	int n=matA->n;
	int k=matB->n;

	TaucsType biv, valA;
	int rowInd, rowA;
	std::vector<std::map<int, TaucsType> > rowsC(k);
	for (int i=0; i<k; ++i) {
		// creating column i of C
		std::map<int, TaucsType> & mapRow2Val = rowsC[i];
		// travel on bi
		for (int rowptrBi = matB->colptr[i];rowptrBi < matB->colptr[i+1];++rowptrBi) {
			rowInd = matB->rowind[rowptrBi];
			biv = matB->taucs_values[rowptrBi];
			// make biv*a_{rowInd} and insert into mapRow2Val
			for (int rowptrA=matA->colptr[rowInd];rowptrA<matA->colptr[rowInd+1];++rowptrA) {
				rowA=matA->rowind[rowptrA];
				valA=matA->taucs_values[rowptrA];
				// insert valA*biv into map
				std::map<int, TaucsType>::iterator it = mapRow2Val.find(rowA);
				if (it == mapRow2Val.end()) {
					// first time
					mapRow2Val[rowA] = valA*biv;
				}
				else {
					it->second = it->second + valA*biv;
				}
			}
		}
		// now column i is created
	}

	return CreateTaucsMatrixFromColumns(rowsC,m,TAUCS_DOUBLE);
}

// For usage when it's known that the result is symmetric, like A^T * A.
taucs_ccs_matrix* TaucsAddOn::Mul2NonSymmMatSymmResult(const taucs_ccs_matrix* matA,
													   const taucs_ccs_matrix* matB)
{
	// Compatibility of dimensions        
	if ((matA->m != matB->n) || (matA->n != matB->m))
		return NULL;

	if ((matA->flags & TAUCS_SYMMETRIC) ||
		(matB->flags & TAUCS_LOWER) )
		return NULL;

	// (m x n)*(n x m) = (m x m)
	int m=matA->m;
	int n=matA->n;

	TaucsType biv, valA;
	int rowInd, rowA;
	std::vector<std::map<int, TaucsType> > rowsC(m);
	for (int i=0; i<m; ++i) {
		// creating column i of C
		std::map<int, TaucsType> & mapRow2Val = rowsC[i];
		// travel on bi
		for (int rowptrBi = matB->colptr[i];rowptrBi < matB->colptr[i+1];++rowptrBi) {
			rowInd = matB->rowind[rowptrBi];
			biv = matB->taucs_values[rowptrBi];
			// make biv*a_{rowInd} and insert into mapRow2Val
			// Ignore anything above the diagonal!!
			for (int rowptrA=matA->colptr[rowInd];rowptrA<matA->colptr[rowInd+1];++rowptrA) {
				rowA=matA->rowind[rowptrA];
				if (rowA >= i) {
					valA=matA->taucs_values[rowptrA];
					// insert valA*biv into map
					std::map<int, TaucsType>::iterator it = mapRow2Val.find(rowA);
					if (it == mapRow2Val.end()) {
						// first time
						mapRow2Val[rowA] = valA*biv;
					}
					else {
						it->second = it->second + valA*biv;
					}
				}
			}
		}
		// now column i is created
	}

	return CreateTaucsMatrixFromColumns(rowsC,m,TAUCS_DOUBLE|TAUCS_SYMMETRIC|TAUCS_LOWER);
}


/// Computes the transpose of a matrix.
taucs_ccs_matrix* TaucsAddOn::MatrixTranspose(const taucs_ccs_matrix* mat)
{
	taucs_ccs_matrix* ret;
	ret = taucs_ccs_create(mat->n, mat->m, mat->colptr[mat->n], mat->flags);
	if (! ret)
		return NULL;


	if (mat->flags & TAUCS_SYMMETRIC) {
		// symmetric - just copy the matrix

		memcpy(ret->colptr, mat->colptr, sizeof(int) * (mat->n + 1));
		memcpy(ret->rowind, mat->rowind, sizeof(int) * (mat->colptr[mat->n]));
		memcpy(ret->taucs_values, mat->taucs_values, sizeof(TaucsType) * (mat->colptr[mat->n]));

		return ret;
	}

	// non-symmetric matrix -> need to build data structure.
	// we'll go over the columns and build the rows
	std::vector< std::vector<int> >       rows(mat->m);
	std::vector< std::vector<TaucsType> > values(mat->m);
	for (int c = 0; c < mat->n; ++c) {
		for (int rowi = mat->colptr[c]; rowi < mat->colptr[c+1]; ++rowi) {
			rows[mat->rowind[rowi]].push_back(c);
			values[mat->rowind[rowi]].push_back(mat->taucs_values[rowi]);
		}
	}

	// copying the rows as columns in ret
	int cind = 0;
	for (int r = 0; r < mat->m; ++r) {
		ret->colptr[r] = cind;

		for (int j = 0; j < (int)rows[r].size(); ++j) {
			ret->rowind[cind] = rows[r][j];
			ret->taucs_values[cind] = values[r][j];
			cind++;
		}
	}
	ret->colptr[mat->m] = cind;

	assert(cind == mat->colptr[mat->n]);

	return ret;
}


taucs_ccs_matrix* TaucsAddOn::CreateTaucsMatrixFromColumns(const vector< map<int, TaucsType> >& cols, 
														   int nRows,
														   int flags)
{
	// count nnz:
	int nCols = (int)cols.size();

	int nnz = 0;
	for (int counter=0; counter < nCols; ++counter) {
		nnz += (int)cols[counter].size();
	}

	taucs_ccs_matrix *matC = taucs_ccs_create(nRows,nCols,nnz,flags);
	if (! matC)
		return NULL;

	// copy cols into matC
	std::map<int,TaucsType>::const_iterator rit;
	int rowptrC = 0;
	for (int c=0;c<nCols;++c) {
		matC->colptr[c] = rowptrC;
		for (rit = cols[c].begin();rit!= cols[c].end();++rit) {
			matC->rowind[rowptrC]=rit->first;
			int ind = rit->first;
			matC->taucs_values[rowptrC]=rit->second;
			double val = rit->second;
			++rowptrC;
		}
	}
	matC->colptr[nCols]=nnz;
	return matC;
}


// Multiplies matA by x and stores the result in b. Assumes all memory has 
// been allocated and the sizes match; assumes matA is not symmetric!!
void TaucsAddOn::MulNonSymmMatrixVector(const taucs_ccs_matrix* matA,
									   const TaucsType* x,
									   TaucsType* b)
{
	// make b all zero
	memset(b, 0, matA->m * sizeof(TaucsType));

	for (int col = 0; col < matA->n; ++col) {
		// going over column col of matA, multiplying it by x[col] and 
		// setting the appropriate values of vector b
		for (int p = matA->colptr[col]; p < matA->colptr[col+1]; ++p) {
			b[matA->rowind[p]] += x[col]*matA->taucs_values[p];
		}
	}
}


//////////////////////////////////////////////////////////////////////////
// Adds two vectors vecA and VecB and stores the result in vecResult.
// Assumes all memory has been allocated and the sizes match!
void TaucsAddOn::Add2Vectors(int n, 
							 const TaucsType* vecA, 
							 const TaucsType* vecB, 
							 TaucsType* vecResult)
{
	for (int i=0; i<n; i++)
		vecResult[i] = vecA[i] + vecB[i];
}

// Copy mat to a new matrix, memory will be allocated during the copy process.
taucs_ccs_matrix* TaucsAddOn::MatrixCopy(const taucs_ccs_matrix *mat) 
{
	taucs_ccs_matrix* ret;
	ret = taucs_ccs_create(mat->m, mat->n, mat->colptr[mat->n], mat->flags);
	if (! ret)
		return NULL;

	memcpy(ret->colptr, mat->colptr, sizeof(int) * (mat->n + 1));
	memcpy(ret->rowind, mat->rowind, sizeof(int) * (mat->colptr[mat->n]));
	memcpy(ret->taucs_values, mat->taucs_values, sizeof(TaucsType) * (mat->colptr[mat->n]));

	return ret;
}
