
#ifndef ARCBALL_H
#define ARCBALL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if defined(WIN32)
  #include <windows.h>
#elif defined(__APPLE__)
struct POINT {
    int x;
    int y;
};
#endif
#include "main.h"

class CArcBall  
{
public:
	typedef enum AxisSet
	{
		NoAxes, CameraAxes, BodyAxes, OtherAxes, NSets
	} AxisSet;
private:
	typedef struct Quat_t
	{
		float x, y, z, w;
	} Quat_t;

	typedef Quat_t HVect;

	enum QuatPart
	{
		X, Y, Z, W, QuatLen
	};

	typedef float HMatrix[QuatLen][QuatLen];

	typedef float *ConstraintSet;

	typedef struct BallData
	{
		HVect		center;
		double		radius;
		Quat_t		qNow, qDown, qDrag;
		HVect		vNow, vDown, vFrom, vTo, vrFrom, vrTo;
		HMatrix		mNow, mDown, mDeltaNow;
		bool		showResult, dragging;
		ConstraintSet	sets[NSets];
		int		setSizes[NSets];
		AxisSet		axisSet;
		int		axisIndex;
	} BallData;

private: // variables
	HMatrix mId;
	float otherAxis[1][4];
	Quat_t qOne;

	BallData ball_data;
	HMatrix ball_matrix;

private: // functions
	void MakeIdentity(float *fpTheMatrix);
	// aux
	HMatrix* Qt_ToMatrix(Quat_t q, HMatrix& out);
	Quat_t Qt_Conj(Quat_t q);
	Quat_t Qt_Mul(Quat_t qL, Quat_t qR);
	HVect V3_(float x, float y, float z);
	float V3_Norm(HVect v);
	HVect V3_Unit(HVect v);
	HVect V3_Scale(HVect v, float s);
	HVect V3_Negate(HVect v);
	HVect V3_Add(HVect v1, HVect v2);
	HVect V3_Sub(HVect v1, HVect v2);
	float V3_Dot(HVect v1, HVect v2);
	HVect V3_Cross(HVect v1, HVect v2);
	HVect V3_Bisect(HVect v0, HVect v1);

	// math
	HVect MouseOnSphere(HVect mouse, HVect ballCenter, double ballRadius);
	HVect ConstrainToAxis(HVect loose, HVect axis);
	int NearestConstraintAxis(HVect loose, HVect *axes, int nAxes);
	Quat_t Qt_FromBallPoints(HVect from, HVect to);
	void Qt_ToBallPoints(Quat_t q, HVect *arcFrom, HVect *arcTo);

	// ball
	void Place(BallData* ball, HVect center, double radius);
	void UseSet(BallData* ball, AxisSet axisSet);
	void ShowResult(BallData* ball);
	void HideResult(BallData* ball);
	void Value();
	void Mouse(HVect vNow);
	void Update();
	void BeginDrag();
	void EndDrag();
	HVect convert(int x, int y, int nWinWidth, int nWinHeight);
	void Init();

public:
	float* GetInvertedBallMatrix();
	void InitBall();
	void MouseMove(POINT point, int nWinWidth, int nWinHeight);
	void MouseUp(POINT point, int nWinWidth, int nWinHeight);
	void MouseDown( POINT point, int nWinWidth, int nWinHeight);
	float* GetBallMatrix()
	{
		float* ptr = (float*)ball_matrix;
		return	ptr;
	}
	CArcBall();
	virtual ~CArcBall();
	void DumpArcBall(FILE* file);
	bool LoadArcBall(string fileName);
};

#endif // ARCBALL_H
