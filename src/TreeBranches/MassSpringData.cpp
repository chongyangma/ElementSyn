
#include "MassSpringData.h"

CMass::CMass()
{
	m_srcIdx = -1;
	m_pos = Vec3f(0.0f, 0.0f, 0.0f);
	m_vel = Vec3f(0.0f, 0.0f, 0.0f);
	m_quat = Vec4f(1.0f, 0.0f, 0.0f, 0.0f);
	m_flagFixed = false;
}

void CMass::TimeIntegrationMass(const Flt dt)
{
	if ( m_flagFixed == true )
	{
		return;
	}
	m_vel += m_force * dt;
	m_pos += m_vel * dt;
}

GLUquadric* CMassSpringData::m_ptrQuadric = NULL;

CMassSpringData::CMassSpringData()
{
	m_srcIdx = -1;
	m_startIdx = 0;
}

void CMassSpringData::ComputeForce(const float springLength, const float stiffness, const float springFriction, const float gravity, const float airFriction)
{
	const int numOfMasses = int(m_vecMass.size());
	for ( int i=0; i<numOfMasses; i++ )
	{
		m_vecMass[i].ResetForce();
	}
	for ( int i=0; i<numOfMasses-1; i++ )
	{
		ComputeSpringForce(i, springLength, stiffness, springFriction);
	}
	for ( int i=0; i<numOfMasses; i++ )
	{
		m_vecMass[i].AddForce(Vec3f(0.0f, gravity, 0.0f));
		m_vecMass[i].AddForce(-m_vecMass[i].GetVel() * airFriction);
	}
}

void CMassSpringData::ComputeSpringForce(const int idx, const float springLength, const float stiffness, const float springFriction)
{
	CMass& mass1 = m_vecMass[idx];
	CMass& mass2 = m_vecMass[idx+1];
	Vec3f springVec = mass1.GetPos() - mass2.GetPos();
	Flt r = mag(springVec);
	Vec3f force = Vec3f(0.0f, 0.0f, 0.0f);
	if ( r > 0.0f )
	{
		force += springVec * (1.0f / r) * (r - springLength) * (-stiffness);
	}
	force += -(mass1.GetVel() - mass2.GetVel()) * springFriction;
	mass1.AddForce(force);
	mass2.AddForce(-force);
}

void CMassSpringData::TimeIntegrationMassSpring(const Flt dt)
{
	const int numOfMasses = int(m_vecMass.size());
	for ( int i=0; i<numOfMasses; i++ )
	{
		m_vecMass[i].TimeIntegrationMass(dt);
	}
}

void CMassSpringData::RenderMassSpringData(Vec3f trans)
{
	int numOfMasses = int(m_vecMass.size());
	if ( m_ptrQuadric == NULL ) m_ptrQuadric = gluNewQuadric();
	const Flt cylinderWd = 0.004f;
	const Flt sphereRad = 0.01f;
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	GLfloat objAmbient[]  = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat objDiffuse[]  = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat objSpecular[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	GLfloat objHLAmbient[]  = { 1.0f, 0.0f, 0.0f, 1.0f };
	GLfloat objHLDiffuse[]  = { 0.5f, 0.0f, 0.0f, 1.0f };
	GLfloat objHLSpecular[] = { 0.1f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT,  objAmbient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,  objDiffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, objSpecular);
	//glMateriali(GL_FRONT, GL_SHININESS, 0.5f);
	glPushMatrix();
	glTranslatef(trans[0], trans[1], trans[2]);
	for ( int i=0; i<numOfMasses-1; i++ )
	{
		Vec3f& p1 = m_vecMass[i].GetPos();
		Vec3f& p2 = m_vecMass[i+1].GetPos();
		Vec3f pMean = (p1 + p2) * 0.5f;
		Vec3f pd = p2 - p1;
		normalize(pd);
		Vec3f v1 = Vec3f(1.0f, 0.0f, 0.0f);
		Vec3f v2 = cross(v1, pd);
		normalize(v2);
		Vec3f v3 = cross(pd, v2);
		Flt mat[16];
		for ( int n=0; n<16; n++ ) mat[n] = 0.0f;
		mat[15] = 1.0f;
		mat[0] = v2[0]; mat[1] = v2[1]; mat[2] = v2[2];
		mat[4] = v3[0]; mat[5] = v3[1]; mat[6] = v3[2];
		mat[8] = pd[0]; mat[9] = pd[1]; mat[10]= pd[2];
		glPushMatrix();
		glTranslatef(p1[0], p1[1], p1[2]);
		glMultMatrixf(mat);
		//gluCylinder(m_ptrQuadric, cylinderWd, cylinderWd, dist(p1, p2), 16, 8);
		det::RenderCylinder(cylinderWd, cylinderWd, dist(p1, p2), 16);
		glPopMatrix();
	}
	glMaterialfv(GL_FRONT, GL_AMBIENT,  objHLAmbient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,  objHLDiffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, objHLSpecular);
	//glMateriali(GL_FRONT, GL_SHININESS, 0.5f);
	for ( int i=0; i<numOfMasses; i++ )
	{
		Vec3f& pi = m_vecMass[i].GetPos();
		glPushMatrix();
		glTranslatef(pi[0], pi[1], pi[2]);
		//gluSphere(m_ptrQuadric, sphereRad, 10, 10);
		glPopMatrix();
	}
	glPopMatrix();
	glDisable(GL_LIGHTING);
}

CMassSpringData CMassSpringData::NurbsInterpolateStrand(int interval /* = 5 */)
{
	int numOfMasses = (GetNumOfMasses()-1) * interval + 1;
	CMassSpringData strandNew;
	strandNew.ResizeMassSpringData(numOfMasses);
	CNURBSinterpolate nurbsInterp;
	vector<Flt> vecShapePts(GetNumOfMasses());
	vector<Vec3f> vecPosNew(numOfMasses);
	for ( int i=0; i<3; i++ )
	{
		for ( int j=0; j<GetNumOfMasses(); j++ )
		{
			vecShapePts[j] = GetMass(j).GetPos()[i];
		}
		vector<Flt> vecInterpPts = nurbsInterp.InterpolateUniformOpen(vecShapePts, interval);
		for ( int j=0; j<numOfMasses; j++ )
		{
			vecPosNew[j][i] = vecInterpPts[j];
		}
	}
	for ( int i=0; i<numOfMasses; i++ )
	{
		strandNew.GetMass(i).SetPos(vecPosNew[i]);
	}
	return strandNew;
}

CMassSpringData CMassSpringData::DownsampleStrand(int interval /* = 5 */)
{
	int numOfMasses = (GetNumOfMasses()-1) / interval + 1;
	CMassSpringData strandNew;
	strandNew.ResizeMassSpringData(numOfMasses);
	CNURBSinterpolate nurbsInterp;
	vector<Vec3f> vecPosNew(numOfMasses);
	for ( int i=0; i<numOfMasses; i++ )
	{
		Vec3f pi = GetMass(i*interval).GetPos();
		strandNew.GetMass(i).SetPos(pi);
	}
	return strandNew;
}

CMassSpringData CMassSpringData::TranslateStrand(Vec3f trans)
{
	int numOfMasses = GetNumOfMasses();
	CMassSpringData strandNew;
	strandNew.ResizeMassSpringData(numOfMasses);
	for ( int i=0; i<numOfMasses; i++ )
	{
		Vec3f pi = GetMass(i).GetPos();
		pi += trans;
		strandNew.GetMass(i).SetPos(pi);
	}
	return strandNew;
}

CMassSpringData CMassSpringData::DeformStrand(Flt d)
{
	int numOfMasses = GetNumOfMasses();
	CMassSpringData strandNew;
	strandNew.ResizeMassSpringData(numOfMasses);
	for ( int i=0; i<numOfMasses; i++ )
	{
		int idx0 = max(i - 1, 0);
		int idx1 = min(i + 1, numOfMasses - 1);
		Vec3f pr = GetMass(idx1).GetPos() - GetMass(idx0).GetPos();
		pr = cross(pr, Vec3f(0.0f, 0.0f, 1.0f));
		normalize(pr);
		pr *= d;
		Vec3f pi = GetMass(i).GetPos();
		pi += pr;
		strandNew.GetMass(i).SetPos(pi);
	}
	return strandNew;
}

CMassSpringData CMassSpringData::FixLengthAsNonStretchable()
{
	int numOfMasses = GetNumOfMasses();
	Flt yMax = -1e10;
	int midIdx = -1;
	for ( int i=0; i<numOfMasses; i++ )
	{
		Flt yTmp = GetMass(i).GetPos()[1];
		if ( yTmp > yMax )
		{
			yMax = yTmp;
			midIdx = i;
		}
	}
	CMassSpringData strandNew;
	strandNew.ResizeMassSpringData(numOfMasses);
	Vec3f p0 = GetMass(midIdx).GetPos();
	strandNew.GetMass(midIdx).SetPos(p0);
	Vec3f p1 = GetMass(midIdx-1).GetPos();
	strandNew.GetMass(midIdx-1).SetPos(p1);
	Flt length1 = dist(p0, p1);
	for ( int i=midIdx-2; i>=0; i-- )
	{
		Vec3f pr = GetMass(i).GetPos() - GetMass(i+1).GetPos();
		Flt l1 = (mag(pr) + length1) * 0.5f;
		normalize(pr);
		pr *= l1; //length1;
		p1 = p1 + pr;
		strandNew.GetMass(i).SetPos(p1);
	}
	Vec3f p2 = GetMass(midIdx+1).GetPos();
	strandNew.GetMass(midIdx+1).SetPos(p2);
	Flt length2 = dist(p0, p2);
	for ( int i=midIdx+2; i<numOfMasses; i++ )
	{
		Vec3f pr = GetMass(i).GetPos() - GetMass(i-1).GetPos();
		Flt l2 = (mag(pr) + length2) * 0.5f;
		normalize(pr);
		pr *= l2; //length2;
		p2 = p2 + pr;
		strandNew.GetMass(i).SetPos(p2);
	}
	return strandNew;
}

vector<CMassSpringData> CMassSpringData::LoadStrandsFromTXT(string fileName, int maxNumOfStrands /* = -1 */)
{
	vector<CMassSpringData> strands;
	ifstream fin(fileName.c_str());
	if ( fin.fail() == true )
	{
		cout << "Failed to load strands from TXT file " << fileName << "!\n";
		return strands;
	}
	int numOfStrands;
	int numOfMasses;
	fin >> numOfStrands;
	numOfStrands = (maxNumOfStrands > 0) ? maxNumOfStrands : numOfStrands; // For debug purpose
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData strand;
		fin >> numOfMasses;
		strand.ResizeMassSpringData(numOfMasses);
		for ( int i=0; i<numOfMasses; i++ )
		{
			CMass& mass = strand.GetMass(i);
			Flt px, py, pz;
			fin >> px >> py >> pz;
			mass.SetPos(Vec3f(px, py, pz));
			Flt vx, vy, vz;
			fin >> vx >> vy >> vz;
			mass.SetVel(Vec3f(vx, vy, vz));
			Flt fx, fy, fz;
			fin >> fx >> fy >> fz;
			mass.SetForce(Vec3f(fx, fy, fz));
		}
		strands.push_back(strand);
	}
	return strands;
}

vector<CMassSpringData> CMassSpringData::NurbsInterpolateStrands(vector<CMassSpringData>& strands, int interval /* = 5 */)
{
	vector<CMassSpringData> strandsNew(strands.size());
	for ( int i=0; i<int(strands.size()); i++ )
	{
		strandsNew[i] = (interval == 0) ? strands[i] : strands[i].NurbsInterpolateStrand(interval);
	}
	return strandsNew;
}

vector<CMassSpringData> CMassSpringData::DownsampleStrands(vector<CMassSpringData>& strands, int interval /* = 5 */)
{
	vector<CMassSpringData> strandsNew(strands.size());
	for ( int i=0; i<int(strands.size()); i++ )
	{
		strandsNew[i] = (interval == 0) ? strands[i] : strands[i].DownsampleStrand(interval);
	}
	return strandsNew;
}

bool CMassSpringData::DumpStrandsToTXT(string fileName, vector<CMassSpringData>& strands, const Vec3f trans)
{
	ofstream fout(fileName.c_str());
	if ( fout.fail() == true )
	{
		cout << "Failed to dump strands to file " << fileName << "!\n";
		return false;
	}
	int numOfStrands = int(strands.size());
	fout << numOfStrands << endl;
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& msData = strands[n];
		int numOfMasses = msData.GetNumOfMasses();
		fout << numOfMasses << endl;
		for ( int i=0; i<numOfMasses; i++ )
		{
			CMass& mass = msData.GetMass(i);
			Vec3f pos = mass.GetPos() + trans;
			Vec3f vel = mass.GetVel();
			Vec3f force = mass.GetForce();
			fout << pos[0] << "\t" << pos[1] << "\t" << pos[2] << endl;
			fout << vel[0] << "\t" << vel[1] << "\t" << vel[2] << endl;
			fout << force[0] << "\t" << force[1] << "\t" << force[2] << endl;
		}
	}
	return true;
}
bool CMassSpringData::DumpStrandsToOBJ(vector<CMassSpringData>& strands, string fileName /* =  */, Flt rad /* = 0.01f */, int sampleNum /* = 8 */, int interval /* = 1 */, int interp /* = 0 */)
{
	ofstream fout(fileName.c_str());
	if ( fout.fail() )
	{
		cout << "Failed to dump spaghetti to " << fileName << "!\n";
		return false;
	}

	int vAccum = 0;
	for ( int i=0; i<int(strands.size()); i++ )
	{
		CMassSpringData strand = (interp == 0) ? strands[i] : strands[i].NurbsInterpolateStrand(5);
		DumpStrandToOBJ(strand, fout, vAccum, rad, sampleNum, interval);
	}

	return true;
}

void CMassSpringData::DumpStrandToOBJ(CMassSpringData& strand, ofstream& fout, int& vAccum, Flt rad, int sampleNum, int interval)
{
	vector<Vec3f> vecVertices;
	int vCnt = 0;
	int ptNum = strand.GetNumOfMasses();
	Flt angleDelta = TWO_PI / Flt(sampleNum);
	Vec3f pr, pi;
	int iMax = -1;
	for ( int i=0; i<ptNum; i+=interval )
	{
		pi = strand.GetMass(i).GetPos();
		iMax = max(i, iMax);
		if ( i < ptNum - interval )
		{
			pr = strand.GetMass(i+interval).GetPos() - pi;
		}
		Vec3f d1 = cross(pr, Vec3f(0.f, 0.f, 1.f));
		Vec3f d2 = cross(d1, pr);
		normalize(d1);
		normalize(d2);
		Flt angleVal = 0.f;
		for ( int j=0; j<sampleNum; j++ )
		{
			Flt cv = cos(angleVal);
			Flt sv = sin(angleVal);
			Vec3f pj = pi + rad * (cv * d1 + sv * d2);
			fout << "v " << pj << endl;
			vecVertices.push_back(pj);
			angleVal += angleDelta;
			vCnt ++;
		}
	}
	Flt s = 0.5f;
	Vec3f v0 = (1.f + s) * strand.GetMass(0).GetPos() - s * strand.GetMass(interval).GetPos();
	Vec3f v1 = (1.f + s) * pi - s * strand.GetMass(iMax - interval).GetPos();
	fout << "v " << v0 << endl;
	fout << "v " << v1 << endl;
	vecVertices.push_back(v0);
	vecVertices.push_back(v1);

	int sampleCnt = 0;
	for ( int i=0; i<ptNum-interval; i+=interval )
	{
		vector<int> vecIndices0;
		vector<int> vecIndices1;
		vector<int> vecIndices2;
		vector<int> vecIndices3;
		for ( int j=0; j<sampleNum; j++ )
		{
			int idx0 = vAccum + j + 1 + sampleCnt;
			int idx1 = vAccum + (j + 1) % sampleNum + 1 + sampleCnt;
			int idx2 = idx0 + sampleNum;
			int idx3 = idx1 + sampleNum;
			vecIndices0.push_back(idx0);
			vecIndices1.push_back(idx1);
			vecIndices2.push_back(idx2);
			vecIndices3.push_back(idx3);
		}
		Flt distSumMin = 1e10;
		int bestShift = 0;
		for ( int j=0; j<sampleNum; j++ )
		{
			Flt distSumTmp = 0.f;
			for ( int k=0; k<sampleNum; k++ )
			{
				Vec3f& p0 = vecVertices[vecIndices0[k] - vAccum - 1];
				Vec3f& p2 = vecVertices[vecIndices2[(k + j) % sampleNum] - vAccum - 1];
				distSumTmp += dist2(p0, p2);
			}
			if ( distSumTmp < distSumMin )
			{
				distSumMin = distSumTmp;
				bestShift = j;
			}
		}
		for ( int j=0; j<sampleNum; j++ )
		{
			int idx0 = vecIndices0[j];
			int idx1 = vecIndices1[j];
			int idx2 = vecIndices2[(j + bestShift) % sampleNum];
			int idx3 = vecIndices3[(j + bestShift) % sampleNum];
			fout << "f " << idx0 << " " << idx3 << " " << idx1 << endl;
			fout << "f " << idx0 << " " << idx2 << " " << idx3 << endl;
		}
		sampleCnt += sampleNum;
	}

	// Close the two ends...
	for ( int j=0; j<sampleNum; j++ )
	{
		int idx0 = vAccum + vCnt + 1;
		int idx1 = vAccum + j + 1;
		int idx2 = vAccum + (j + 1) % sampleNum + 1;
		fout << "f " << idx0 << " " << idx1 << " " << idx2 << endl;
	}
	for ( int j=0; j<sampleNum; j++ )
	{
		int idx0 = vAccum + vCnt + 2;
		int idx1 = vAccum + vCnt - sampleNum + (j + 1) % sampleNum + 1;
		int idx2 = vAccum + vCnt - sampleNum + j + 1;
		fout << "f " << idx0 << " " << idx1 << " " << idx2 << endl;
	}

	vCnt += 2;
	vAccum += vCnt;
}

bool CMassSpringData::DumpStrandsToOBJz(vector<CMassSpringData>& strands, string fileName /* =  */, Flt rad /* = 0.01f */, int sampleNum /* = 8 */, int interval /* = 1 */, int interp /* = 0 */)
{
	ofstream fout(fileName.c_str());
	if ( fout.fail() )
	{
		cout << "Failed to dump spaghetti to " << fileName << "!\n";
		return false;
	}

	int vAccum = 0;
	for ( int i=0; i<int(strands.size()); i++ )
	{
		CMassSpringData strand = (interp == 0) ? strands[i] : strands[i].NurbsInterpolateStrand(5);
		DumpStrandToOBJz(strand, fout, vAccum, rad, sampleNum, interval);
	}

	return true;
}

void CMassSpringData::DumpStrandToOBJz(CMassSpringData& strand, ofstream& fout, int& vAccum, Flt rad, int sampleNum, int interval)
{
	vector<Vec3f> vecVertices;
	int vCnt = 0;
	int ptNum = strand.GetNumOfMasses();
	Flt angleDelta = TWO_PI / Flt(sampleNum);
	Vec3f pr, pi;
	int iMax = -1;
	for ( int i=0; i<ptNum; i+=interval )
	{
		pi = strand.GetMass(i).GetPos();
		iMax = max(i, iMax);
		if ( i < ptNum - interval )
		{
			pr = strand.GetMass(i+interval).GetPos() - pi;
		}
		Vec3f d1 = cross(pr, Vec3f(0.f, 1.f, 0.f));//cross(pr, Vec3f(0.f, 0.f, 1.f));
		Vec3f d2 = cross(d1, pr);
		normalize(d1);
		normalize(d2);
		Flt angleVal = 0.f;
		for ( int j=0; j<sampleNum; j++ )
		{
			Flt cv = cos(angleVal);
			Flt sv = sin(angleVal);
			Vec3f pj = pi + rad * (cv * d1 + sv * d2);
			fout << "v " << pj << endl;
			vecVertices.push_back(pj);
			angleVal += angleDelta;
			vCnt ++;
		}
	}
	Flt s = 0.0f; //0.5f;
	Vec3f v0 = (1.f + s) * strand.GetMass(0).GetPos() - s * strand.GetMass(interval).GetPos();
	Vec3f v1 = (1.f + s) * pi - s * strand.GetMass(iMax - interval).GetPos();
	fout << "v " << v0 << endl;
	fout << "v " << v1 << endl;
	vecVertices.push_back(v0);
	vecVertices.push_back(v1);

	int sampleCnt = 0;
	for ( int i=0; i<ptNum-interval; i+=interval )
	{
		vector<int> vecIndices0;
		vector<int> vecIndices1;
		vector<int> vecIndices2;
		vector<int> vecIndices3;
		for ( int j=0; j<sampleNum; j++ )
		{
			int idx0 = vAccum + j + 1 + sampleCnt;
			int idx1 = vAccum + (j + 1) % sampleNum + 1 + sampleCnt;
			int idx2 = idx0 + sampleNum;
			int idx3 = idx1 + sampleNum;
			vecIndices0.push_back(idx0);
			vecIndices1.push_back(idx1);
			vecIndices2.push_back(idx2);
			vecIndices3.push_back(idx3);
		}
		Flt distSumMin = 1e10;
		int bestShift = 0;
		for ( int j=0; j<sampleNum; j++ )
		{
			Flt distSumTmp = 0.f;
			for ( int k=0; k<sampleNum; k++ )
			{
				Vec3f& p0 = vecVertices[vecIndices0[k] - vAccum - 1];
				Vec3f& p2 = vecVertices[vecIndices2[(k + j) % sampleNum] - vAccum - 1];
				distSumTmp += dist2(p0, p2);
			}
			if ( distSumTmp < distSumMin )
			{
				distSumMin = distSumTmp;
				bestShift = j;
			}
		}
		for ( int j=0; j<sampleNum; j++ )
		{
			int idx0 = vecIndices0[j];
			int idx1 = vecIndices1[j];
			int idx2 = vecIndices2[(j + bestShift) % sampleNum];
			int idx3 = vecIndices3[(j + bestShift) % sampleNum];
			fout << "f " << idx0 << " " << idx3 << " " << idx1 << endl;
			fout << "f " << idx0 << " " << idx2 << " " << idx3 << endl;
		}
		sampleCnt += sampleNum;
	}

	// Close the two ends...
	for ( int j=0; j<sampleNum; j++ )
	{
		int idx0 = vAccum + vCnt + 1;
		int idx1 = vAccum + j + 1;
		int idx2 = vAccum + (j + 1) % sampleNum + 1;
		fout << "f " << idx0 << " " << idx1 << " " << idx2 << endl;
	}
	for ( int j=0; j<sampleNum; j++ )
	{
		int idx0 = vAccum + vCnt + 2;
		int idx1 = vAccum + vCnt - sampleNum + (j + 1) % sampleNum + 1;
		int idx2 = vAccum + vCnt - sampleNum + j + 1;
		fout << "f " << idx0 << " " << idx1 << " " << idx2 << endl;
	}

	vCnt += 2;
	vAccum += vCnt;
}
