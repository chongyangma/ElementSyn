
#ifndef BULLETSPAGHETTI_H
#define BULLETSPAGHETTI_H

#include "SpaghettiSynConfig.h"

#include "GlutDemoApplication.h"
#include "LinearMath/btAlignedObjectArray.h"

#ifdef _DEBUG
#pragma comment(lib, "BulletCollision_vs2010_debug.lib")
#pragma comment(lib, "BulletDynamics_vs2010_debug.lib")
#pragma comment(lib, "BulletSoftBody_vs2010_debug.lib")
#pragma comment(lib, "ConvexDecomposition_vs2010_debug.lib")
#pragma comment(lib, "HACD_vs2010_debug.lib")
#pragma comment(lib, "LinearMath_vs2010_debug.lib")
#pragma comment(lib, "OpenGLSupport_vs2010_debug.lib")
#else
#pragma comment(lib, "BulletCollision_vs2010.lib")
#pragma comment(lib, "BulletDynamics_vs2010.lib")
#pragma comment(lib, "BulletSoftBody_vs2010.lib")
#pragma comment(lib, "ConvexDecomposition_vs2010.lib")
#pragma comment(lib, "HACD_vs2010.lib")
#pragma comment(lib, "LinearMath_vs2010.lib")
#pragma comment(lib, "OpenGLSupport_vs2010.lib")
#endif

class btBroadphaseInterface;
class btCollisionShape;
class btOverlappingPairCache;
class btCollisionDispatcher;
class btConstraintSolver;
struct btCollisionAlgorithmCreateFunc;
class btDefaultCollisionConfiguration;

class CBulletSpaghetti : public GlutDemoApplication
{
public:
	void initPhysics() {}

	void InitBulletSpaghetti(vector<CMassSpringData>& strands);

	void InitBulletSpaghettiNew(vector<CMassSpringData>& strands); // For final output

	void exitPhysics();

	virtual ~CBulletSpaghetti()
	{
		exitPhysics();
	}

	virtual void clientMoveAndDisplay();

	void EditSpaghettiFromBullet(vector<CMassSpringData>& strands);

	void EditSpaghettiFromBulletNew(vector<CMassSpringData>& strands); // For final output

private:
	btAlignedObjectArray<class CSpaghetti*> m_vecSpaghetti;

	//keep the collision shapes, for deletion/cleanup
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;

	btBroadphaseInterface*	m_broadphase;

	btCollisionDispatcher*	m_dispatcher;

	btConstraintSolver*	m_solver;

	btDefaultCollisionConfiguration* m_collisionConfiguration;
};

#endif BULLETSPAGHETTI_H
