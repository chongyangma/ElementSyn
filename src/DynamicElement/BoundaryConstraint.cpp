
#include "BoundaryConstraint.h"
#include "lodepng.h"

float CBoundaryConstraint::m_colorThresh = 0.5f;

CBoundaryConstraint::CBoundaryConstraint()
{
	m_bcPosMin = Vec3f(-1.0f,-1.0f,-1.0f);
	m_bcPosMax = Vec3f( 1.0f, 1.0f, 1.0f);
	m_bcFrameIdx = -1;
	m_bcInputIdx = 0;
}

bool CBoundaryConstraint::ApplyBoundaryConstraint(Vec3f& posOld, Vec3f& posNew)
{
	bool flag = InsideBoundary(posOld);
	if ( flag == true )
	{
		int idxMin = -1;
		Flt distSqrMin = 1e10;
		for ( int i=0; i<int(m_vecBoundarySample.size()); i++ )
		{
			Vec3f& pi = m_vecBoundarySample[i].m_samplePos;
			Flt distSqrTmp = dist2(posOld, pi);
			if ( distSqrTmp < distSqrMin )
			{
				idxMin = i;
				distSqrMin = distSqrTmp;
			}
		}
		Vec3f pos = m_vecBoundarySample[idxMin].m_samplePos;
		Vec3f norm = m_vecBoundarySample[idxMin].m_sampleNorm;
		Vec3f pr = posOld - pos;
		Flt d = dot(pr, norm);
		pr = pr - d * norm;
		posNew = pos + pr;
	}
	return flag;
}

bool CBoundaryConstraint::ApplyBoundaryConstraintNew(Vec3f& posOld, Vec3f& posNew)
{
	Vec3f cubicCen = (m_bcPosMin + m_bcPosMax) * 0.5f;
	Vec3f cubicSize = (m_bcPosMax - m_bcPosMin);
	Flt scaling = 1.0f / max(cubicSize[0], cubicSize[1]);
	Vec3f pos = posOld;
	pos = pos - cubicCen;
	pos *= scaling;
	pos += Vec3f(0.5f, 0.5f, 0.0f);
	Vec3f posTmp = pos;
	bool flag = ApplyBoundaryConstraint(pos, posTmp);
	if ( flag == true )
	{
		posNew = posTmp - Vec3f(0.5f, 0.5f, 0.0f);
		posNew *= 1.0f / scaling;
		posNew += cubicCen;
	}
	return flag;
}

bool CBoundaryConstraint::InsideBoundary(Vec3f& pos)
{
	return (m_grid.Sample(pos[0], pos[1])[0] < m_colorThresh);
}

bool CBoundaryConstraint::InsideBoundaryNew(Vec3f& pos)
{
	Vec3f cubicCen = (m_bcPosMin + m_bcPosMax) * 0.5f;
	Vec3f cubicSize = (m_bcPosMax - m_bcPosMin);
	Flt scaling = 1.0f / max(cubicSize[0], cubicSize[1]);
	Vec3f posNew = pos;
	posNew = posNew - cubicCen;
	posNew *= scaling;
	posNew += Vec3f(0.5f, 0.5f, 0.0f);
	return InsideBoundary(posNew);
}

bool CBoundaryConstraint::InsideBoundary(int px, int py)
{
	return (m_grid(px, py)[0] < m_colorThresh);
}

void CBoundaryConstraint::LoadBoundaryConstraintFromImage(string fileName)
{
	m_grid = LoadGrid2D3fFromCImage(fileName.c_str());
	m_gridContour = IdentifyGridContour(m_grid);
	SaveGrid2D3fAsCImage(m_gridContour, "contour.png");
	Triangulate();
}

CGrid2D3f CBoundaryConstraint::LoadGrid2D3fFromCImage(const char* fileName)
{
	CGrid2D3f grid;

	std::vector<unsigned char> image;
	unsigned width, height;
	unsigned error = lodepng::decode(image, width, height, fileName);
	if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;

	grid.Allocate(width, height);

	for ( unsigned j=0; j<height; j++  )
	{
		unsigned jn = height-1-j;
		for ( unsigned i=0; i<width; i++ )
		{
			for ( int d=0; d<3; d++ )
			{
				grid(i, jn)[d] = image[4 * (j * width + i) + d] / 255.0;
			}
		}
	}
	return grid;
}

void CBoundaryConstraint::SaveGrid2D3fAsCImage(const CGrid2D3f& grid, const char* fileName)
{
	Vec2ui size = grid.GetSize();
	std::vector<unsigned char> vecByte(4 * size[0] * size[1]);
	for ( uint j=0; j<size[1]; j++ )
	{
		for ( uint i=0; i<size[0]; i++ )
		{
			for ( int d=0; d<3; d++ )
			{
				vecByte[4 * (j * size[0] + i) + d] = (unsigned char)(grid(i, j)[d] * 255.0);
			}
		}
	}
	unsigned error = lodepng::encode(fileName, vecByte, size[0], size[1]);
	if(error) std::cout << "PNG encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

CGrid2D3f CBoundaryConstraint::IdentifyGridContour(const CGrid2D3f& gridInput)
{
	Vec2ui sizeInput = gridInput.GetSize();
	CGrid2D3f gridContour;
	gridContour.Allocate(sizeInput);
	vector<Vec2i> vecContourPoint;
	for (unsigned int i=0; i<sizeInput[0]; i++ )
	{
		for (unsigned int j=0; j<sizeInput[1]; j++ )
		{
			Vec3f clContour = Vec3f(0.0f, 0.0f, 0.0f);
			bool contourFlag = false;
			if ( InsideBoundary(i, j) == false )
			{
				int windowSize = 3;
				for ( int di=-windowSize/2; di<=windowSize/2; di++ )
				{
					for ( int dj=-windowSize/2; dj<=windowSize/2; dj++ )
					{
						int ni = (i+di+sizeInput[0]) % sizeInput[0];
						int nj = (j+dj+sizeInput[1]) % sizeInput[1];
						Vec3f clNew = gridInput(ni, nj);
						if ( InsideBoundary(ni, nj) == true )
						{
							clContour = Vec3f(1.0f, 1.0f, 1.0f);
							contourFlag = true;
						}
					}
				}
			}
			if ( contourFlag == true )
			{
				vecContourPoint.push_back(Vec2i(i, j));
			}
			gridContour(i, j) = clContour;
		}
	}
	cout << "Number of contour points: " << vecContourPoint.size() << endl;
	// Sample the contour...
	const float sampleDistSqr = 100;
	vector<BoundarySample> vecBoundarySample;
	// Detect corner samples first...
	for ( int i=0; i<int(vecContourPoint.size()); i++ )
	{
		Vec2i coord = vecContourPoint[i];
		int d = 2;
		int cnt = 0;
		if ( InsideBoundary(coord[0]-d, coord[1]-d) == false ) cnt ++;
		if ( InsideBoundary(coord[0]-d, coord[1]+d) == false ) cnt ++;
		if ( InsideBoundary(coord[0]+d, coord[1]-d) == false ) cnt ++;
		if ( InsideBoundary(coord[0]+d, coord[1]+d) == false ) cnt ++;
		if ( cnt == 1 || cnt == 3 ) // A corner point!
		{
			bool flag = false;
			for ( int n=0; n<int(vecBoundarySample.size()); n++ )
			{
				Vec3f pn = vecBoundarySample[n].m_samplePos;
				if ( (coord[0]-pn[0]) * (coord[0]-pn[0]) + (coord[1]-pn[1]) * (coord[1]-pn[1]) < sampleDistSqr )
				{
					flag = true;
					break;
				}
			}
			if ( flag == false )
			{
				BoundarySample sample;
				sample.m_samplePos = Vec3f(coord[0], coord[1], 0.0f);
				sample.m_sampleNorm = ComputeSampleNorm(coord[0], coord[1]);
				vecBoundarySample.push_back(sample);
			}
		}
	}
	cout << "Number of corner points: " << vecBoundarySample.size() << endl;
	// Add non-corner samples...
	for ( int i=0; i<int(vecContourPoint.size()); i++ )
	{
		Vec2i coord = vecContourPoint[i];
		bool flag = false;
		for ( int n=0; n<int(vecBoundarySample.size()); n++ )
		{
			Vec3f pn = vecBoundarySample[n].m_samplePos;
			if ( (coord[0]-pn[0]) * (coord[0]-pn[0]) + (coord[1]-pn[1]) * (coord[1]-pn[1]) < sampleDistSqr )
			{
				flag = true;
				break;
			}
		}
		if ( flag == false )
		{
			BoundarySample sample;
			sample.m_samplePos = Vec3f(coord[0], coord[1], 0.0f);
			sample.m_sampleNorm = ComputeSampleNorm(coord[0], coord[1]);
			vecBoundarySample.push_back(sample);
		}
	}
	for ( int n=0; n<int(vecBoundarySample.size()); n++ )
	{
		Vec3f pn = vecBoundarySample[n].m_samplePos;
		Vec2i coord = Vec2i(int(pn[0]), int(pn[1]));
		gridContour(coord[0], coord[1]) = Vec3f(0.0f, 0.0f, 1.0f);
	}
	cout << "Number of contour sample points: " << int(vecBoundarySample.size()) << endl;
	m_vecBoundarySample = vecBoundarySample;
	return gridContour;
}

void CBoundaryConstraint::Triangulate()
{
	vector<Vec2f> vecPoint;
	ofstream fout("tri.obj");
	for ( int i=0; i<int(m_vecBoundarySample.size()); i++ )
	{
		Vec3f pi = m_vecBoundarySample[i].m_samplePos;
		vecPoint.push_back(Vec2f(pi[0], pi[1]));
		fout << "v " << pi << endl;
	}
	CDelaunayTri triangulation;
	triangulation.InitMesh(vecPoint);
	triangulation.IncrementalDelaunay();
	MESH_PTR ptrMesh = triangulation.GetMeshPtr();
	TRIANGLE_PTR ptrTri = ptrMesh->pTriArr;
	while ( ptrTri != NULL )
	{
		int idx1 = ptrTri->i1 - 2;
		int idx2 = ptrTri->i2 - 2;
		int idx3 = ptrTri->i3 - 2;
		Vec3f p1 = m_vecBoundarySample[idx1-1].m_samplePos;
		Vec3f p2 = m_vecBoundarySample[idx2-1].m_samplePos;
		Vec3f p3 = m_vecBoundarySample[idx3-1].m_samplePos;
		Vec3f p = (p1 + p2 + p3) / 3.0f;
		bool insideFlag = InsideBoundary(p[0], p[1]);
		if ( insideFlag == true )
		{
			fout << "f " << idx1 << " " << idx2 << " " << idx3 << endl;
		}
		ptrTri = ptrTri->pNext;
	}
}

Vec3f CBoundaryConstraint::ComputeSampleNorm(int px, int py)
{
	Vec3f x1 = Vec3f(0.0f, 0.0f, 0.0f);
	const int nnTotal = 360;
	const float nnInv = 1.0f / float(nnTotal);
	float fx = px / float(m_grid.GetSize()[0]);
	float fy = py / float(m_grid.GetSize()[1]);
	for ( int nn=0; nn<nnTotal; nn++ )
	{
		float angleValue = atan(1.0f) * 8.0 * nn * nnInv;
		float cv = cos(angleValue);
		float sv = sin(angleValue);
		Vec3f pr = Vec3f(cv, sv, 0.0f);
		Vec3f p0 = Vec3f(fx, fy, 0.0f) + pr * 0.0001f;
		if ( InsideBoundary(p0) )
		{
			x1 = x1 + pr;
		}
	}
	if ( mag2(x1) > 0.0f )
	{
		normalize(x1);
	}
	return x1;
}
