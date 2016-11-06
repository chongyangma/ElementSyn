
#ifndef HUNGARIANALGORITHM_H
#define HUNGARIANALGORITHM_H

#define CHECK_FOR_INF
#define ONE_INDEXING
#define INF 10000000

#include "main.h"

class CHungarianAlgorithm
{
public:
	void FindOptimalAssignment(vector<int>& assignment, double* cost, vector<double>& distMatrixIn, const int nOfRows, const int nOfColumns);

	void FindOptimalAssignmentOld(vector<int>& assignment, double* cost, vector<double>& distMatrixIn, const int nOfRows, const int nOfColumns);

private:
	inline double haGetInf() { return 1e10; }
	inline bool haIsFinite(double d) { return d < haGetInf(); }

	void SetMatCost(vector<vector<double>>& matCost);

	void InitLabel(int n);

	void UpdateLabel();

	void AddToTree(int i, int previ);

	void Augment();

	inline vector<int>& GetXY() { return m_vecXY; }

	int m_N; // number of vertices in one part
	int m_maxMatch;
	vector<double> m_vecLx, m_vecLy; // labels of X and Y parts
	vector<int> m_vecXY, m_vecYX;
	vector<bool> m_vecS, m_vecT;
	vector<double> m_vecSlack;
	vector<double> m_vecSlackx;
	vector<int> m_vecPrev; // to record alternating paths
	vector<vector<double>> m_matCost; // cost matrix
};

#endif HUNGARIANALGORITHM_H
