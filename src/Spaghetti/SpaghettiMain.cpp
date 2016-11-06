
#include "../DynamicElement/ArcBall.h"
#include "../DynamicElement/Camera.h"
#include "SpaghettiSyn.h"

bool pause = false; //true; //
int g_width = 720;
int g_height = 480;
int g_bmpCount = 0;
Flt g_mouseOldPx;
Flt g_mouseOldPy;
CArcBall* ptrArcBall = NULL;
CCamera* ptrCamera = NULL;
CSpaghettiSyn* ptrSynthesizer = NULL;

void SaveSnapshot(bool flag = true)
{
	if ( ptrSynthesizer->GetStepCount() < 1 && pause == true && flag == true ) return;
	if ( ptrSynthesizer->GetStepCount() > CSynConfigBase::m_stepCountMax && pause == false ) return;
	int wd_mod = g_width % 4;
	int wd = (wd_mod == 0) ? g_width : (g_width + 4 - wd_mod);
	GLvoid* ptrVoid = new BYTE[3 * wd * g_height];
	glReadPixels(0, 0, g_width, g_height, GL_BGR_EXT, GL_UNSIGNED_BYTE, ptrVoid);
	char fileName[MAX_PATH];
	sprintf_s(fileName, "%sSnapshot\\snapshot_%04d.png", CSpaghettiSynConfig::m_outputPrefix.c_str(), ptrSynthesizer->GetStepCount());
	CImage dstImg;
	dstImg.Create(g_width, g_height, 24);
	BYTE* ptrDst = (BYTE*)dstImg.GetBits();
	BYTE* ptrSrc = (BYTE*)ptrVoid;
	int dstPitch = dstImg.GetPitch();
	for ( int j=0; j<g_height; j++ )
	{
		for ( int i=0; i<g_width; i++ )
		{
			*(ptrDst + dstPitch * (g_height - 1 - j) + 3 * i + 0) = *(ptrSrc + 0);
			*(ptrDst + dstPitch * (g_height - 1 - j) + 3 * i + 1) = *(ptrSrc + 1);
			*(ptrDst + dstPitch * (g_height - 1 - j) + 3 * i + 2) = *(ptrSrc + 2);
			ptrSrc += 3;
		}
	}
	dstImg.Save(fileName);
	delete [] ptrVoid;
}

void DumpCameraAndArcBall()
{
	string fileName = "CameraAndArcBall.txt";
	FILE* file;
	file = fopen(fileName.c_str(), "w");
	if ( !file )
	{
		cout << "Failed to dump camera and arc-ball into config file " << fileName << "!\n";
		return;
	}
	ptrCamera->DumpCamera(file);
	ptrArcBall->DumpArcBall(file);
	fclose(file);
	cout << "Have dumped camera and arc-ball into config file " << fileName << "!\n";
}

void LoadCameraAndArcBall()
{
	string fileName = "CameraAndArcBall.txt";
	ptrCamera->LoadCamera(fileName);
	ptrArcBall->LoadArcBall(fileName);
	cout << "Have loaded camera and arc-ball from config file " << fileName << "!\n";
}

void DisplayFunc()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);  // White background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, float(g_width)/float(g_height), 0.01, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	// Set look-at matrix
	ptrCamera->Draw();
	// Rotate according to the arc-ball
	glPushMatrix();
	glMultMatrixf(ptrArcBall->GetBallMatrix());	
	// Draw the cubic framework...
	glColor3f(1.f, 1.f, 1.f);
	//glutWireCube(2.0);
	ptrSynthesizer->RenderOutput();
	glPopMatrix();
	SaveSnapshot();
	glutSwapBuffers();
	if ( ptrSynthesizer->GetStepCount() == CSpaghettiSynConfig::m_stepCountMax ) ptrSynthesizer->RestartSynthesis();
}

void IdleFunc()
{
	if ( pause == false )
	{
		ptrSynthesizer->UpdateOutput();
		//cout << "Step: " << ptrSynthesizer->GetStepCount() << endl;
	}
	glutPostRedisplay();
}

void ReshapeFunc(int width, int height)
{
	if( height == 0 )
	{
		height = 1;
	}
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, float(width)/float(height), 0.01, 1000);
	glMatrixMode(GL_MODELVIEW);
}

void KeyboardFunc(unsigned char key, int x, int y) 
{
	switch(key)
	{
	case 'a':
		ptrCamera->MoveLeft(0.1f);
		break;
	case 'd':
		ptrCamera->MoveRight(0.1f);
		break;
	case 'w':
		ptrCamera->MoveForward(0.1f);
		break;
	case 's':
		ptrCamera->MoveBackward(0.1f);
		break;
	case 'r':
		ptrCamera->MoveUp(0.1f);
		break;
	case 'f':
		ptrCamera->MoveDown(0.1f);
		break;
	case 'q':
		ptrCamera->RotateY(-3);
		break;
	case 'e':
		ptrCamera->RotateY(3);
		break;
	case 't':
		ptrCamera->RotateX(3);
		break;
	case 'g':
		ptrCamera->RotateX(-3);
		break;
	case 'n':
		ptrSynthesizer->UpdateOutput();
		break;
	case 'u':
		DumpCameraAndArcBall();
		break;
	case 'l':
		LoadCameraAndArcBall();
		break;
	case 'p':
		SaveSnapshot(false);
		break;
	case ' ':
		pause = !pause;
		break;
	case 27: // 'Esc'
		exit(0);
		break;
	default:
		break;
	}
}

void SpecialFunc(int key, int x, int y)
{
}

void MouseFunc(int button, int state, int x, int y)
{
	x = g_width  - x;
	y = g_height - y;
	POINT tPoint;
	tPoint.x = x;
	tPoint.y = y;
	if ( state == GLUT_DOWN )
	{
		ptrArcBall->MouseDown(tPoint, g_width, g_height);
	}
	else if ( state == GLUT_UP )
	{
		ptrArcBall->MouseUp(tPoint, g_width, g_height);
	}
	DisplayFunc();
}

void MouseMoveFunc(int x, int y)
{
	x = g_width  - x;
	y = g_height - y;
	POINT tPoint;
	tPoint.x = x;
	tPoint.y = y;
	ptrArcBall->MouseMove(tPoint, g_width, g_height);
	DisplayFunc();
}

void Initialize()
{
	ptrArcBall = new CArcBall;
	ptrArcBall->InitBall();
	ptrCamera = new CCamera(0, 0, 3);
	ptrCamera->RotateY(180);
	ptrSynthesizer = new CSpaghettiSyn;

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA); // set display mode
	glutInitWindowSize(g_width, g_height); // set window size
	glutInitWindowPosition(0, 0); // set window position on screen
	glutCreateWindow("Spaghetti"); // set window title

	glutMouseFunc(MouseFunc);
	glutMotionFunc(MouseMoveFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutSpecialFunc(SpecialFunc);
	glutDisplayFunc(DisplayFunc);
	glutReshapeFunc(ReshapeFunc);
	glutIdleFunc(IdleFunc);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glPointSize(5.f);
	glLineWidth(1.f);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	LoadCameraAndArcBall();
}

void ReleasePtr()
{
	DELETE_OBJECT(ptrArcBall);
	DELETE_OBJECT(ptrCamera);
	DELETE_OBJECT(ptrSynthesizer);
}

int main(int argc,char **argv)
{
	glutInit(&argc,argv);
	Initialize();
	glutMainLoop();
	ReleasePtr();
	return 0;
}
