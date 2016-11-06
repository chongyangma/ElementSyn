
#define CONSTRAINT_DEBUG_SIZE 0.2f

#include "btBulletDynamicsCommon.h"
#include "GlutStuff.h"
#include "GL_ShapeDrawer.h"

#include "LinearMath/btIDebugDraw.h"

#include "GLDebugDrawer.h"
#include "BulletSpaghetti.h"

class CSpaghetti
{
	btDynamicsWorld* m_ownerWorld;
	vector<btCollisionShape*> m_vecShapes;
	vector<btRigidBody*> m_vecBodies;
	vector<btTypedConstraint*> m_vecJoints;
	int m_numOfMasses;

	btRigidBody* localCreateRigidBody(btScalar mass, const btTransform& startTransform, btCollisionShape* shape)
	{
		bool isDynamic = (mass != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic)
		{
			shape->calculateLocalInertia(mass, localInertia);
		}

		btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

		btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, shape, localInertia);
		btRigidBody* body = new btRigidBody(rbInfo);
		body->setFriction(1.0f);

		m_ownerWorld->addRigidBody(body);

		return body;
	}

public:
	CSpaghetti(btDynamicsWorld* ownerWorld, vector<float>& vecPx, vector<float>& vecPy, vector<float>& vecPz)
		: m_ownerWorld (ownerWorld)
	{
		m_numOfMasses = int(vecPx.size());
		int numOfCapsuleShapes = m_numOfMasses - 1;
		btScalar massRadius = 0.01f;
		btScalar capsuleHeight = 0.05f;

		// Setup the geometry
		m_vecShapes.resize(numOfCapsuleShapes, NULL);
		for ( int i=0; i<numOfCapsuleShapes; i++ )
		{
			m_vecShapes[i] = new btCapsuleShape(massRadius, capsuleHeight);
		}

		// Setup all the rigid bodies
		m_vecBodies.resize(numOfCapsuleShapes, NULL);

		for ( int i=0; i<numOfCapsuleShapes; i++ )
		{
			float px = (vecPx[i] + vecPx[i+1]) * 0.5f;
			float py = (vecPy[i] + vecPy[i+1]) * 0.5f;
			float pz = (vecPz[i] + vecPz[i+1]) * 0.5f;
			btTransform transform;
			transform.setIdentity();
			transform.setOrigin(btVector3(btScalar(px), btScalar(py), btScalar(pz)));
			float mass = (i < 1) ? 0.0f : 1.0f;
			m_vecBodies[i] = localCreateRigidBody(mass, transform, m_vecShapes[i]);
		}

		// Setup some damping on the bodies
		for ( int i=0; i<numOfCapsuleShapes; i++ )
		{
			m_vecBodies[i]->setDamping(0.05, 0.85);
			m_vecBodies[i]->setDeactivationTime(0.8);
			//m_vecBodies[i]->setSleepingThresholds(1.6, 2.5);
		}

		// Now setup the constraints
		int numOfJoints = numOfCapsuleShapes - 1;
		m_vecJoints.resize(numOfJoints, NULL);
		for ( int i=0; i<numOfJoints; i++ )
		{
			Vec3f p0 = Vec3f(vecPx[i], vecPy[i], vecPz[i]);
			Vec3f p2 = Vec3f(vecPx[i+2], vecPy[i+2], vecPz[i+2]);
			Vec3f pr = (p2 - p0) * 0.5f;
			Vec3f prHalf = pr * 0.5f;
			btTransform localA, localB;
			localA.setIdentity();
			localB.setIdentity();
			//localA.setOrigin(btVector3(btScalar(0.),-btScalar(capsuleHeight*0.5f), btScalar(0.)));
			//localB.setOrigin(btVector3(btScalar(0.), btScalar(capsuleHeight*0.5f), btScalar(0.)));
			localA.setOrigin(btVector3( prHalf[0], prHalf[1], prHalf[2]));
			localB.setOrigin(btVector3(-prHalf[0],-prHalf[1],-prHalf[2]));
			btHingeConstraint* hingeC = new btHingeConstraint(*m_vecBodies[i], *m_vecBodies[i+1], localA, localB);
			hingeC->setLimit(btScalar(-HALF_PI), btScalar(HALF_PI));
			m_vecJoints[i] = hingeC;
			hingeC->setDbgDrawSize(CONSTRAINT_DEBUG_SIZE);
			m_ownerWorld->addConstraint(m_vecJoints[i], true);
		}
	}

	virtual	~CSpaghetti()
	{
		// Remove all constraints
		for ( int i=0; i<int(m_vecJoints.size()); i++ )
		{
			m_ownerWorld->removeConstraint(m_vecJoints[i]);
			delete m_vecJoints[i];
			m_vecJoints[i] = NULL;
		}

		// Remove all bodies and shapes
		for ( int i=0; i<int(m_vecBodies.size()); ++i)
		{
			m_ownerWorld->removeRigidBody(m_vecBodies[i]);

			delete m_vecBodies[i]->getMotionState();
			delete m_vecBodies[i];
			m_vecBodies[i] = 0;
			delete m_vecShapes[i];
			m_vecShapes[i] = 0;
		}
	}
};

void CBulletSpaghetti::InitBulletSpaghetti(vector<CMassSpringData>& strands)
{
	m_collisionConfiguration = new btDefaultCollisionConfiguration();

	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);

	btVector3 worldAabbMin(-10000,-10000,-10000);
	btVector3 worldAabbMax(10000,10000,10000);
	m_broadphase = new btAxisSweep3 (worldAabbMin, worldAabbMax);

	m_solver = new btSequentialImpulseConstraintSolver;

	m_dynamicsWorld = new btDiscreteDynamicsWorld(m_dispatcher,m_broadphase,m_solver,m_collisionConfiguration);
	//m_dynamicsWorld->getDispatchInfo().m_useConvexConservativeDistanceUtil = true;
	//m_dynamicsWorld->getDispatchInfo().m_convexConservativeDistanceThreshold = 0.01f;
	m_dynamicsWorld->setGravity(btVector3(0,CSpaghettiSynConfig::m_gravity,0));

	int numOfStrands = int(strands.size());
	for ( int n=0; n<numOfStrands; n++ )
	{
		int numOfMasses = strands[n].GetNumOfMasses();
		vector<float> vecPx(numOfMasses);
		vector<float> vecPy(numOfMasses);
		vector<float> vecPz(numOfMasses);
		for ( int i=0; i<numOfMasses; i++ )
		{
			Vec3f pi = strands[n].GetMass(i).GetPos();
			vecPx[i] = pi[0];
			vecPy[i] = pi[1];
			vecPz[i] = pi[2];
		}
		CSpaghetti* ptrSpaghetti = new CSpaghetti(m_dynamicsWorld, vecPx, vecPy, vecPz);
		m_vecSpaghetti.push_back(ptrSpaghetti);
	}

	clientResetScene();
}

void CBulletSpaghetti::InitBulletSpaghettiNew(vector<CMassSpringData>& strands)
{
	m_collisionConfiguration = new btDefaultCollisionConfiguration();

	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);

	btVector3 worldAabbMin(-10000,-10000,-10000);
	btVector3 worldAabbMax(10000,10000,10000);
	m_broadphase = new btAxisSweep3 (worldAabbMin, worldAabbMax);

	m_solver = new btSequentialImpulseConstraintSolver;

	m_dynamicsWorld = new btDiscreteDynamicsWorld(m_dispatcher,m_broadphase,m_solver,m_collisionConfiguration);
	//m_dynamicsWorld->getDispatchInfo().m_useConvexConservativeDistanceUtil = true;
	//m_dynamicsWorld->getDispatchInfo().m_convexConservativeDistanceThreshold = 0.01f;
	m_dynamicsWorld->setGravity(btVector3(0,CSpaghettiSynConfig::m_gravity,0));

	int numOfStrands = int(strands.size());
	for ( int n=0; n<numOfStrands; n++ )
	{
		int numOfMasses = strands[n].GetNumOfMasses();
		int midIdx = numOfMasses / 2;
		// The first part...
		int segmentLength1 = midIdx;
		vector<float> vecPx1(segmentLength1);
		vector<float> vecPy1(segmentLength1);
		vector<float> vecPz1(segmentLength1);
		for ( int i=0; i<segmentLength1; i++ )
		{
			Vec3f pi = strands[n].GetMass(midIdx-1-i).GetPos();
			vecPx1[i] = pi[0];
			vecPy1[i] = pi[1];
			vecPz1[i] = pi[2];
		}
		CSpaghetti* ptrSpaghetti1 = new CSpaghetti(m_dynamicsWorld, vecPx1, vecPy1, vecPz1);
		m_vecSpaghetti.push_back(ptrSpaghetti1);
		// The second part...
		int segmentLength2 = numOfMasses - midIdx - 1;
		vector<float> vecPx2(segmentLength2);
		vector<float> vecPy2(segmentLength2);
		vector<float> vecPz2(segmentLength2);
		for ( int i=0; i<segmentLength2; i++ )
		{
			Vec3f pi = strands[n].GetMass(midIdx+1+i).GetPos();
			vecPx2[i] = pi[0];
			vecPy2[i] = pi[1];
			vecPz2[i] = pi[2];
		}
		CSpaghetti* ptrSpaghetti2 = new CSpaghetti(m_dynamicsWorld, vecPx2, vecPy2, vecPz2);
		m_vecSpaghetti.push_back(ptrSpaghetti2);
	}

	clientResetScene();
}

void CBulletSpaghetti::exitPhysics()
{
	for ( int i=0; i<m_vecSpaghetti.size(); i++ )
	{
		CSpaghetti* ptrSpaghetti = m_vecSpaghetti[i];
		delete ptrSpaghetti;
	}

	//cleanup in the reverse order of creation/initialization

	//remove the rigid bodies from the dynamics world and delete them

	for ( int i=m_dynamicsWorld->getNumCollisionObjects()-1; i>=0; i-- )
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject( obj );
		delete obj;
	}

	//delete collision shapes
	for ( int j=0; j<m_collisionShapes.size(); j++ )
	{
		btCollisionShape* shape = m_collisionShapes[j];
		delete shape;
	}

	//delete dynamics world
	delete m_dynamicsWorld;

	//delete solver
	delete m_solver;

	//delete broadphase
	delete m_broadphase;

	//delete dispatcher
	delete m_dispatcher;

	delete m_collisionConfiguration;
}

void CBulletSpaghetti::clientMoveAndDisplay()
{
	if (m_dynamicsWorld)
	{
		m_dynamicsWorld->stepSimulation(1.0f / 60.0f);
	}
}

void CBulletSpaghetti::EditSpaghettiFromBullet(vector<CMassSpringData>& strands)
{
	btCollisionObjectArray objects = m_dynamicsWorld->getCollisionObjectArray();
	int numOfStrands = int(strands.size());
	int objCnt = 0;
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& strand = strands[n];
		int numOfMasses = strand.GetNumOfMasses();
		for ( int i=1; i<numOfMasses; i++ )
		{
			Vec3f posLast = strand.GetMass(i-1).GetPos();
			btCollisionObject* obj = objects[objCnt];
			btRigidBody* body = btRigidBody::upcast(obj);
			btTransform trans;
			body->getMotionState()->getWorldTransform(trans);
			btVector3 pos = trans.getOrigin();
			Vec3f posMid = Vec3f(pos[0], pos[1], pos[2]);
			Vec3f pi = posMid * 2.0f - posLast;
			strand.GetMass(i).SetPos(pi);
			objCnt ++;
		}
	}
}

void CBulletSpaghetti::EditSpaghettiFromBulletNew(vector<CMassSpringData>& strands)
{
	btCollisionObjectArray objects = m_dynamicsWorld->getCollisionObjectArray();
	int numOfStrands = int(strands.size());
	int objCnt = 0;
	for ( int n=0; n<numOfStrands; n++ )
	{
		CMassSpringData& strand = strands[n];
		int numOfMasses = strand.GetNumOfMasses();
		int midIdx = numOfMasses / 2;
		// First part...
		for ( int i=1; i<midIdx; i++ )
		{
			Vec3f posLast = strand.GetMass(midIdx-i).GetPos();
			btCollisionObject* obj = objects[objCnt];
			btRigidBody* body = btRigidBody::upcast(obj);
			btTransform trans;
			body->getMotionState()->getWorldTransform(trans);
			btVector3 pos = trans.getOrigin();
			Vec3f posMid = Vec3f(pos[0], pos[1], pos[2]);
			Vec3f pi = posMid * 2.0f - posLast;
			strand.GetMass(midIdx-i-1).SetPos(pi);
			objCnt ++;
		}
		// Second part...
		for ( int i=1; i<midIdx-1; i++ )
		{
			Vec3f posLast = strand.GetMass(midIdx+i).GetPos();
			btCollisionObject* obj = objects[objCnt];
			btRigidBody* body = btRigidBody::upcast(obj);
			btTransform trans;
			body->getMotionState()->getWorldTransform(trans);
			btVector3 pos = trans.getOrigin();
			Vec3f posMid = Vec3f(pos[0], pos[1], pos[2]);
			Vec3f pi = posMid * 2.0f - posLast;
			strand.GetMass(midIdx+i+1).SetPos(pi);
			objCnt ++;
		}
	}
}
