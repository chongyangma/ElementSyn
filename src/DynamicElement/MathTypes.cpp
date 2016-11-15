
#include "MathTypes.h"

// Implementation of CMatrixBase
CMatrixBase::CMatrixBase()
{
	m_row = 0;
	m_col = 0;
}

// Implementation of CDenseMatrix
CDenseMatrix::CDenseMatrix()
{
}

CDenseMatrix::CDenseMatrix(int row, int col)
{
	ResizeMatrix(row, col);
}

void CDenseMatrix::ResizeMatrix(int row, int col)
{
	m_row = row;
	m_col = col;
	m_vals.clear();
	m_vals.resize(row);
	for ( int i=0; i<m_row; i++ )
	{
		m_vals[i].resize(m_col);
	}
	for ( int i=0; i<m_row; i++ )
	{
		for ( int j=0; j<m_col; j++ )
		{
			m_vals[i][j] = 0.f;
		}
	}
}

void CDenseMatrix::PrintMatrix()
{
	for ( int i=0; i<m_row; i++ )
	{
		for ( int j=0; j<m_col; j++ )
		{
			cout << m_vals[i][j] << " ";
		}
		cout << endl;
	}
}

void CDenseMatrix::UpdateDiagonalVal(int i, Flt dv)
{
	m_vals[i][i] += dv;
}

void CDenseMatrix::UpdatePairVals(int i, int j, Flt dv)
{
	m_vals[i][i] += dv;
	m_vals[j][j] += dv;
	m_vals[i][j] -= dv;
	m_vals[j][i] -= dv;
}

void CDenseMatrix::UpdateMatVal(int i, int j, Flt dv)
{
	m_vals[i][j] += dv;
}

// Implementation of CSparseMatrixCRS
CSparseMatrixCRS::CSparseMatrixCRS()
{
}

void CSparseMatrixCRS::ConvertFromDenseMatrix(CDenseMatrix* ptrMatD)
{
	m_row = ptrMatD->GetRowNum();
	m_col = ptrMatD->GetColNum();
	m_vecRowIndex.resize(m_row + 1);
	m_vecColumn.clear();
	m_vecVal.clear();
	int nonZeroCnt = 1;
	vector<vector<Flt> >& vals = ptrMatD->GetVals();
	for ( int i=0; i<m_row; i++ )
	{
		m_vecRowIndex[i] = nonZeroCnt;
		for ( int j=i; j<m_col; j++ )
		{
			if ( vals[i][j] != 0.f )
			{
				m_vecVal.push_back(vals[i][j]);
				m_vecColumn.push_back(j + 1);
				nonZeroCnt ++;
			}
		}
	}
	m_vecRowIndex[m_row] = nonZeroCnt;
}

void CSparseMatrixCRS::ConvertFromCrossList(CCrossList* ptrMatD)
{
	m_row = ptrMatD->GetRowNum();
	m_col = ptrMatD->GetColNum();
	m_vecRowIndex.resize(m_row + 1);
	m_vecColumn.clear();
	m_vecVal.clear();
	int nonZeroCnt = 1;
	for ( int i=1; i<=m_row; i++ )
	{
		m_vecRowIndex[i-1] = nonZeroCnt;
		for ( OLNode* q=ptrMatD->GetPtrRhead()[i]; q!=NULL; q=q->m_ptrR )
		{
			int j = q->m_y;
			m_vecVal.push_back(q->m_value);
			m_vecColumn.push_back(j);
			nonZeroCnt ++;
		}
	}
	m_vecRowIndex[m_row] = nonZeroCnt;
}

// Implementation of CCrossList
CCrossList::CCrossList(int row, int col)
{
	m_row = row;
	m_col = col;
	m_Ne = 0;
	m_ptrRhead = new OLink[m_row + 1];
	if ( m_ptrRhead == NULL )
	{
		printf("Memory allocation error!\n");
	}
	m_ptrChead = new OLink[m_col + 1];
	if ( m_ptrChead == NULL )
	{
		printf("Memory allocation error!\n");
	}
	for ( int i=0; i<m_row+1; i++ )
	{
		m_ptrRhead[i] = NULL;
	}
	for ( int i=0; i<m_col+1; i++ )
	{
		m_ptrChead[i] = NULL;
	}
}

CCrossList::~CCrossList()
{
	for ( int i=1; i<=m_row; i++ )
	{
		OLNode *p, *q;
		p = m_ptrRhead[i];
		while ( p != NULL )
		{
			q = p->m_ptrR;
			DELETE_OBJECT(p);
			p = q;
		}
	}
	DELETE_ARRAY(m_ptrRhead);
	DELETE_ARRAY(m_ptrChead);
}

void CCrossList::UpdateDiagonalVal(int i, Flt dv)
{
	UpdateMatVal(i, i, dv);
}

void CCrossList::UpdatePairVals(int i, int j, Flt dv)
{
	UpdateMatVal(i, i, dv);
	UpdateMatVal(j, j, dv);
	UpdateMatVal(i, j,-dv);
	UpdateMatVal(j, i,-dv);
}

void CCrossList::UpdateMatVal(int i, int j, Flt dv)
{
	i = i + 1;
	j = j + 1;
	OLink pa = m_ptrRhead[i];
	OLink pre = NULL;
	OLink *hl = new OLink[m_col+1];
	for ( int n=1; n<=m_col; n++ )
	{
		hl[n] = m_ptrChead[n];
	}
	while ( 1 )
	{
		if ( pa == NULL || pa->m_y > j )
		{
			// Insert a new node...
			m_Ne ++;
			OLNode* p = new OLNode;
			p->m_x = i;
			p->m_y = j;
			p->m_value = dv;
			if ( pre == NULL )
			{
				m_ptrRhead[p->m_x] = p;
			}
			else
			{
				pre->m_ptrR = p;
			}
			p->m_ptrR = pa;
			pre = p;
			if ( m_ptrChead[p->m_y] == NULL )
			{
				m_ptrChead[p->m_y] = p;
				p->m_ptrD = NULL;
			}
			else
			{
				p->m_ptrD = hl[p->m_y]->m_ptrD;
				hl[p->m_y]->m_ptrD = p;
			}
			hl[p->m_y] = p;
			break;
		}
		else if ( pa != NULL && pa->m_y < j )
		{
			// Move pa to the next node in this row
			pre = pa;
			pa = pa->m_ptrR;
		}
		else if ( pa->m_y == j )
		{
			// Add the current element of M
			pa->m_value += dv;
			if ( pa->m_value == 0 )
			{ // The sum is 0, delete the node...
				m_Ne --;
				if ( pre == NULL )
				{
					m_ptrRhead[pa->m_x] = pa->m_ptrR;
				}
				else
				{
					pre->m_ptrR = pa->m_ptrR;
				}
				OLink p = pa;
				pa = pa->m_ptrR;
				if ( m_ptrChead[p->m_y] == p )
				{
					m_ptrChead[p->m_y] = p->m_ptrD;
					hl[p->m_y] = p->m_ptrD;
				}
				else
				{
					hl[p->m_y]->m_ptrD = p->m_ptrD;
				}
				delete p;
			}
			else
			{
				pre = pa;
				pa = pa->m_ptrR;
			}
			break;
		} // End-If
	} // End-While
	DELETE_ARRAY(hl);
}

void CCrossList::PrintCrossList()
{
	printf("Dimensions: (%d, %d)\n", m_row, m_col);
	printf("# of non-zero elements:%d\n", m_Ne);
	for ( int i=1; i<=m_row; i++ )
	{
		int j=1;
		int y;
		for ( OLNode* q=m_ptrRhead[i]; q!=NULL; q=q->m_ptrR )
		{
			y = q->m_y;
			for ( int n=j; n<y; n++ )
			{
				printf("0 ");
			}
			j = y + 1;
			printf("%f ", q->m_value);
		}
		if ( j<=m_col )
		{
			for ( int n=j; n<=m_col; n++ )
			{
				printf("0 ");
			}
		}
		printf("\n");
	}
	printf("\n");
}
