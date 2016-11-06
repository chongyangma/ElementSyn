
#include "HungarianAlgorithm.h"

void CHungarianAlgorithm::FindOptimalAssignment(vector<int>& assignment, double* cost, vector<double>& distMatrixIn, const int nOfRows, const int nOfColumns)
{
	double distMin = 1e10;
	double distMax =-1e10;
	for ( int i=0; i<int(distMatrixIn.size()); i++ )
	{
		distMin = min(distMin, distMatrixIn[i]);
		distMax = max(distMax, distMatrixIn[i]);
	}
	int n = nOfColumns; // nOfColumns >= nOfRows
	vector<vector<double>> costMat(n);
	for ( int i=0; i<n; i++ )
	{
		vector<double> costVec(n);
		for ( int j=0; j<n; j++ )
		{
			costVec[j] = (j < nOfRows) ? (distMax - distMatrixIn[i*nOfRows+j]) : -(distMax - distMin);
		}
		costMat[i] = costVec;
	}
	SetMatCost(costMat);
	InitLabel(n);
	Augment();
	*cost = 0.0;
	for ( int i=0; i<nOfColumns; i++ )
	{
		int idx = GetXY()[i];
		if ( idx < nOfRows )
		{
			assignment[idx] = i + 1;
			*cost += distMax - costMat[i][GetXY()[i]];
		}
	}
}

void CHungarianAlgorithm::FindOptimalAssignmentOld(vector<int>& assignment, double* cost, vector<double>& distMatrixIn, const int nOfRows, const int nOfColumns)
{
	bool allSinglyValidated, singleValidationFound;
	int n, row, col, tmpRow, tmpCol;
	double value, minValue;

	const double inf = haGetInf();

	// Make working copy of distance matrix
	int nOfElements = nOfRows * nOfColumns;
	vector<double> distMatrix(nOfElements);
	for(n=0; n<nOfElements; n++)
		distMatrix[n] = distMatrixIn[n];

	// Initialization
	*cost = 0;
#ifdef ONE_INDEXING
	for(row=0; row<nOfRows; row++)
		assignment[row] = 0;
#else
	for(row=0; row<nOfRows; row++)
		assignment[row] = -1;
#endif

	// Allocate memory
	vector<int> nOfValidObservations(nOfRows);
	vector<int> nOfValidTracks(nOfColumns);

	// Compute number of validations
	bool infiniteValueFound = false;
	bool finiteValueFound  = false;
	for(row=0; row<nOfRows; row++)
		for(col=0; col<nOfColumns; col++)
			if(haIsFinite(distMatrix[row + nOfRows*col]))
			{
				nOfValidTracks[col]       += 1;
				nOfValidObservations[row] += 1;
				finiteValueFound = true;
			}
			else
				infiniteValueFound = true;

	if(infiniteValueFound)
	{
		if(!finiteValueFound)
			return;

		bool repeatSteps = true;

		while(repeatSteps)
		{
			repeatSteps = false;

			// Step 1: Reject assignments of multiply validated tracks to singly validated observations
			for(col=0; col<nOfColumns; col++)
			{
				singleValidationFound = false;
				for(row=0; row<nOfRows; row++)
					if(haIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidObservations[row] == 1))
					{
						singleValidationFound = true;
						break;
					}

					if(singleValidationFound)
					{
						for(row=0; row<nOfRows; row++)
							if((nOfValidObservations[row] > 1) && haIsFinite(distMatrix[row + nOfRows*col]))
							{
								distMatrix[row + nOfRows*col] = inf;
								nOfValidObservations[row] -= 1;                        
								nOfValidTracks[col]       -= 1;
								repeatSteps = true;            
							}
					}
			}

			// Step 2: Reject assignments of multiply validated observations to singly validated tracks
			if(nOfColumns > 1)          
			{  
				for(row=0; row<nOfRows; row++)
				{
					singleValidationFound = false;
					for(col=0; col<nOfColumns; col++)
						if(haIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidTracks[col] == 1))
						{
							singleValidationFound = true;
							break;
						}

						if(singleValidationFound)
						{
							for(col=0; col<nOfColumns; col++)
								if((nOfValidTracks[col] > 1) && haIsFinite(distMatrix[row + nOfRows*col]))
								{
									distMatrix[row + nOfRows*col] = inf;
									nOfValidObservations[row] -= 1;
									nOfValidTracks[col]       -= 1;
									repeatSteps = true;                            
								}
						}
				}
			}
		} // while(repeatSteps)

		// For each multiply validated track that validates only with singly validated
		// observations, choose the observation with minimum distance
		for(row=0; row<nOfRows; row++)
		{
			if(nOfValidObservations[row] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(col=0; col<nOfColumns; col++)
				{
					value = distMatrix[row + nOfRows*col];
					if(haIsFinite(value))
					{
						if(nOfValidTracks[col] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidTracks[col] == 1) && (value < minValue))
						{
							tmpCol   = col;
							minValue = value;
						}
					}
				}

				if(allSinglyValidated)
				{
#ifdef ONE_INDEXING
					assignment[row] = tmpCol + 1;
#else
					assignment[row] = tmpCol;
#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*tmpCol] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[row + nOfRows*n] = inf;
				}
			}
		}

		// For each multiply validated observation that validates only with singly validated
		// track, choose the track with minimum distance
		for(col=0; col<nOfColumns; col++)
		{
			if(nOfValidTracks[col] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(row=0; row<nOfRows; row++)
				{
					value = distMatrix[row + nOfRows*col];
					if(haIsFinite(value))
					{
						if(nOfValidObservations[row] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidObservations[row] == 1) && (value < minValue))
						{
							tmpRow   = row;
							minValue = value;
						}
					}
				}

				if(allSinglyValidated)
				{
#ifdef ONE_INDEXING
					assignment[tmpRow] = col + 1;
#else
					assignment[tmpRow] = col;
#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*col] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[tmpRow + nOfRows*n] = inf;
				}
			}
		}  
	} // if(infiniteValueFound)

	// Now, recursively search for the minimum element and do the assignment
	while(true)
	{
		// Find minimum distance observation-to-track pair
		minValue = inf;
		for(row=0; row<nOfRows; row++)
			for(col=0; col<nOfColumns; col++)
			{
				value = distMatrix[row + nOfRows*col];
				if(haIsFinite(value) && (value < minValue))
				{
					minValue = value;
					tmpRow   = row;
					tmpCol   = col;
				}
			}

			if(haIsFinite(minValue))
			{
#ifdef ONE_INDEXING
				assignment[tmpRow] = tmpCol+ 1;
#else
				assignment[tmpRow] = tmpCol;
#endif
				*cost += minValue;
				for(n=0; n<nOfRows; n++)
					distMatrix[n + nOfRows*tmpCol] = inf;
				for(n=0; n<nOfColumns; n++)
					distMatrix[tmpRow + nOfRows*n] = inf;
			}
			else
				break;

	} // while(true)
}

void CHungarianAlgorithm::SetMatCost(vector<vector<double>>& matCost)
{
	m_matCost = matCost;
}

void CHungarianAlgorithm::InitLabel(int n)
{
	m_N = n;
	m_maxMatch = 0;

	m_vecLx.resize(m_N);
	m_vecLy.resize(m_N);
	m_vecXY.resize(m_N);
	m_vecYX.resize(m_N);
	m_vecS.resize(m_N);
	m_vecT.resize(m_N);
	m_vecPrev.resize(m_N);
	m_vecSlack.resize(m_N);
	m_vecSlackx.resize(m_N);
	for ( int i=0; i<m_N; i++ )
	{
		m_vecLx[i] = 0;
		m_vecLy[i] = 0;
		m_vecXY[i] = -1;
		m_vecYX[i] = -1;
	}
	for ( int i=0; i<m_N; i++ )
	{
		for ( int j=0; j<m_N; j++ )
		{
			m_vecLx[i] = max(m_vecLx[i], m_matCost[i][j]);
		}
	}
}

void CHungarianAlgorithm::UpdateLabel()
{
	double delta = INF;
	for ( int i=0; i<m_N; i++ )
	{
		if ( m_vecT[i] == false )
		{
			delta = min(delta, m_vecSlack[i]);
		}
	}
	for ( int i=0; i<m_N; i++ )
	{
		if ( m_vecS[i] == true )
		{
			m_vecLx[i] -= delta;
		}
	}
	for ( int i=0; i<m_N; i++ )
	{
		if ( m_vecT[i] == true )
		{
			m_vecLy[i] += delta;
		}
	}
	for ( int i=0; i<m_N; i++ )
	{
		if ( m_vecT[i] == false )
		{
			m_vecSlack[i] -= delta;
		}
	}
}

void CHungarianAlgorithm::AddToTree(int i, int previ)
{
	m_vecS[i] = true; // add x to S
	m_vecPrev[i] = previ; // need this for augmenting
	for ( int j=0; j<m_N; j++ )
	{
		if ( m_vecLx[i] + m_vecLy[j] - m_matCost[i][j] < m_vecSlack[j] )
		{
			m_vecSlack[j] = m_vecLx[i] + m_vecLy[j] - m_matCost[i][j];
			m_vecSlackx[j] = i;
		}
	}
}

void CHungarianAlgorithm::Augment()
{
	if ( m_maxMatch == m_N )
	{
		// matching is already perfect
		return;
	}

	int x, y, root; // counters and root vertex
	vector<int> vecQ;
	vecQ.resize(m_N);
	int wr = 0;
	int rd = 0;

	for ( x=0; x<m_N; x++ )
	{
		m_vecS[x] = false;
		m_vecT[x] = false;
		m_vecPrev[x] = -1;
	}

	// find root of the tree...
	for ( int i=0; i<m_N; i++ )
	{
		if ( m_vecXY[i] == -1 )
		{
			vecQ[wr++] = root = i;
			m_vecPrev[i] = -2;
			m_vecS[i] = true;
			break;
		}
	}

	// initialize slack array...
	for ( y=0; y<m_N; y++ )
	{
		m_vecSlack[y] = m_vecLx[root] + m_vecLy[y] - m_matCost[root][y];
		m_vecSlackx[y] = root;
	}

	while(true)
	{
		while(rd < wr)
		{
			x = vecQ[rd++];
			for ( y=0; y<m_N; y++ )
			{
				if ( m_matCost[x][y] == m_vecLx[x] + m_vecLy[y] && m_vecT[y] == false )
				{
					if ( m_vecYX[y] == -1 )
					{
						break;
					}
					m_vecT[y] = true;
					vecQ[wr++] = m_vecYX[y];
					AddToTree(m_vecYX[y], x);
				}
			}
			if ( y < m_N )
			{
				break;
			}
		}
		if ( y <m_N )
		{
			break;
		}

		UpdateLabel();
		wr = 0;
		rd = 0;
		for ( y=0; y<m_N; y++ )
		{
			if ( m_vecT[y] == false && m_vecSlack[y] == 0 )
			{
				if ( m_vecYX[y] == -1 )
				{
					x = m_vecSlackx[y];
					break;
				}
				else
				{
					m_vecT[y] = true;
					if ( m_vecS[m_vecYX[y]] == false )
					{
						vecQ[wr++] = m_vecYX[y];
						AddToTree(m_vecYX[y], m_vecSlackx[y]);
					}
				}
			}
		}
		if ( y < m_N )
		{
			break;
		}
	}
	if ( y < m_N )
	{
		m_maxMatch ++;
		for ( int cx=x, cy=y, ty; cx!=-2; cx=m_vecPrev[cx], cy=ty )
		{
			ty = m_vecXY[cx];
			m_vecYX[cy] = cx;
			m_vecXY[cx] = cy;
		}
		Augment();
	}
}
