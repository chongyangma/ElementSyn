
#include "MKLsolver.h"

CMKLsolver::CMKLsolver()
{
	m_nRows = 0;
	m_nCols = 0;
	m_nNonZeros = 0;
	m_nRhs = 1;
	m_ptrRowIndex = NULL;
	m_ptrColumns = NULL;
	m_ptrValues = NULL;
	m_ptrSolValues = NULL;
	m_ptrRhs = NULL;
}

CMKLsolver::~CMKLsolver()
{
	DeallocatePtrs();
}

void CMKLsolver::AllocatePtrs(int nRows, int nCols, int nNonZeros)
{
	if ( m_nRows == nRows && m_nCols == nCols && m_nNonZeros == nNonZeros )
	{
		return;
	}
	DeallocatePtrs();
	m_nRows = nRows;
	m_nCols = nCols;
	m_nNonZeros = nNonZeros;
	m_ptrRowIndex = new int[m_nRows + 1];
	m_ptrColumns = new int[m_nNonZeros];
	m_ptrValues = new _DOUBLE_PRECISION_t[m_nNonZeros];
	m_ptrSolValues = new _DOUBLE_PRECISION_t[m_nCols];
	m_ptrRhs = new _DOUBLE_PRECISION_t[m_nRows];
}

void CMKLsolver::DeallocatePtrs()
{
	DELETE_ARRAY(m_ptrRhs);
	DELETE_ARRAY(m_ptrRowIndex);
	DELETE_ARRAY(m_ptrColumns);
	DELETE_ARRAY(m_ptrValues);
	DELETE_ARRAY(m_ptrSolValues);
}

void CMKLsolver::SetCoeffsFromCoeffMatrix(CDenseMatrix* ptrMatD)
{
	CSparseMatrixCRS matS;
	matS.ConvertFromDenseMatrix(ptrMatD);
	SetCoeffsFromCSparseMatCRS(matS);
}

void CMKLsolver::SetCoeffsFromCoeffMatrix(CCrossList* ptrMatD)
{
	CSparseMatrixCRS matS;
	matS.ConvertFromCrossList(ptrMatD);
	SetCoeffsFromCSparseMatCRS(matS);
}

void CMKLsolver::SetCoeffsFromCSparseMatCRS(CSparseMatrixCRS& matS)
{
	int nRows = matS.GetRowNum();
	int nCols = matS.GetColNum();
	int nNonZeros = matS.GetNonZerosNum();
	AllocatePtrs(nRows, nCols, nNonZeros);

	int nRhs = 1;
	vector<int> vecRowIndex = matS.GetVecRowIndex();
	vector<int> vecColumns = matS.GetVecColumn();
	vector<Flt> vecVal = matS.GetVecVal();
	for ( int i=0; i<(nRows + 1); i++ )
	{
		m_ptrRowIndex[i] = vecRowIndex[i];
	}
	for ( int i=0; i<nNonZeros; i++ )
	{
		m_ptrColumns[i] = vecColumns[i];
		m_ptrValues[i] = vecVal[i];
	}
}

vector<Flt> CMKLsolver::GetSolutionViaDCG(vector<Flt>& vecC, vector<Flt>& initGuess)
{
	int rci_request, itercount;
	int ipar[128];
	double dpar[128];
	double* tmp = new double[4 * m_nRows];
	char tr = 'u';
	double eone = -1.;
	int ione = 1;

	// Initialize the right hand side...
	for ( int i=0; i<m_nRows; i++ )
	{
		m_ptrRhs[i] = vecC[i];
	}

	// Initialize the initial guess...
	for ( int i=0; i<m_nCols; i++ )
	{
		m_ptrSolValues[i] = initGuess[i];
	}

	// Initialize the solver
	dcg_init(&m_nRows, m_ptrSolValues, m_ptrRhs, &rci_request, ipar, dpar, tmp);

	// Set the desired parameters...
	// LOGICAL parameters:
	// do residual stopping test
	// do not request for the user defined stopping test
	// DOUBLE parameters
	// set the relative tolerance to 1.0D-5 instead of default value 1.0D-6
	ipar[8]=1;
	ipar[9]=0;
	dpar[0]=1e-5;

	// Check the correctness and consistency of the newly set parameters...
	dcg_check(&m_nRows, m_ptrSolValues, m_ptrRhs, &rci_request, ipar, dpar, tmp);

	// Compute the solution by RCI (P)CG solver without preconditioning...
	rci_request = 1;
	while ( rci_request == 1 )
	{
		dcg(&m_nRows, m_ptrSolValues, m_ptrRhs, &rci_request, ipar, dpar, tmp);
		mkl_dcsrsymv(&tr, &m_nRows, m_ptrValues, m_ptrRowIndex, m_ptrColumns, tmp, &tmp[m_nRows]);
	}

	// Get the current iteration number into itercount
	dcg_get(&m_nRows, m_ptrSolValues, m_ptrRhs, &rci_request, ipar, dpar, tmp, &itercount);
	delete [] tmp;

	//cout << "Number of iterations: " << itercount << endl;

	vector<Flt> solution(m_nCols);
	for ( int i=0; i<m_nCols; i++ )
	{
		solution[i] = m_ptrSolValues[i];
	}

	return solution;
}
