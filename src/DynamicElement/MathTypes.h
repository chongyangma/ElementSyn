
#ifndef MATHTYPES_H
#define MATHTYPES_H

#include "main.h"

class CCrossList;

// type definition for an element (in a dense matrix)
typedef struct MatElem
{
	int m_i, m_j;
	Flt m_val;
} MatElem;

// base class definition for matrix
class CMatrixBase
{
public:
	CMatrixBase();

	int GetRowNum() { return m_row; }

	int GetColNum() { return m_col; }

protected:
	int m_row, m_col;
};

// class definition for dense matrix
class CDenseMatrix : public CMatrixBase
{
public:
	CDenseMatrix();

	CDenseMatrix(int row, int col);

	void ResizeMatrix(int row, int col);

	void PrintMatrix();

	void UpdateDiagonalVal(int i, Flt dv);

	void UpdatePairVals(int i, int j, Flt dv);

	void UpdateMatVal(int i, int j, Flt dv);

	inline vector<vector<Flt> >& GetVals() { return m_vals; }

	inline Flt GetVal(int i, int j) { return m_vals[i][j]; }

	inline void SetVal(int i, int j, Flt val) { m_vals[i][j] = val; }

private:
	vector<vector<Flt> > m_vals;
};

// class definition for sparse matrix (compressed row storage)
class CSparseMatrixCRS : public CMatrixBase
{
public:
	CSparseMatrixCRS();

	void ConvertFromDenseMatrix(CDenseMatrix* ptrMatD);

	void ConvertFromCrossList(CCrossList* ptrMatD);

	vector<int>& GetVecRowIndex() { return m_vecRowIndex; }

	vector<int>& GetVecColumn() { return m_vecColumn; }

	vector<Flt>& GetVecVal() { return m_vecVal; }

	int GetNonZerosNum() { return int(m_vecVal.size()); }

private:
	vector<int> m_vecRowIndex;
	vector<int> m_vecColumn;
	vector<Flt> m_vecVal; // non-zero values
};

// type definition for a node (element in a cross-list)
typedef struct OLNode
{
	int m_x, m_y;
	Flt m_value;
	struct OLNode *m_ptrR, *m_ptrD;
} OLNode;

typedef OLNode* OLink;

// class definition for sparse matrix (cross list storage)
class CCrossList : public CMatrixBase
{
public:
	CCrossList(int row, int col);

	~CCrossList();

	void UpdateDiagonalVal(int i, Flt dv);

	void UpdatePairVals(int i, int j, Flt dv);

	void UpdateMatVal(int i, int j, Flt dv);

	void PrintCrossList();

	OLink* GetPtrRhead() { return m_ptrRhead; }

private:
	int m_Ne;
	OLink *m_ptrRhead, *m_ptrChead;
};

#endif // MATHTYPES_H
