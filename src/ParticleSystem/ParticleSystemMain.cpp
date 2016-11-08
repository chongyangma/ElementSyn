
#include "../DynamicElement/ArcBall.h"
#include "../DynamicElement/Camera.h"
#include "ParticleSystemSyn.h"

bool pause = false; //true; //
int g_width = 720; //640;
int g_height = 480;
int g_bmpCount = 0;
CArcBall* ptrArcBall = NULL;
CCamera* ptrCamera = NULL;
CParticleSystemSyn* ptrSynthesizer = NULL;

void SaveSnapshot()
{
	if ( ptrSynthesizer->GetStepCount() < 1 && pause == true ) return;
	if ( ptrSynthesizer->GetStepCount() > CSynConfigBase::m_stepCountMax && pause == false ) return;
	int wd_mod = g_width % 4;
	int wd = (wd_mod == 0) ? g_width : (g_width + 4 - wd_mod);
	GLvoid* ptrVoid = new BYTE[3 * wd * g_height];
	glReadPixels(0, 0, g_width, g_height, GL_BGR_EXT, GL_UNSIGNED_BYTE, ptrVoid);
	char fileName[MAX_PATH];
	sprintf_s(fileName, "%sSnapshot\\snapshot_%04d.png", CSynConfigBase::m_outputPrefix.c_str(), ptrSynthesizer->GetStepCount());
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
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	// Set look-at matrix
	ptrCamera->Draw();
	// Rotate according to the arc-ball
	glPushMatrix();
	glMultMatrixf(ptrArcBall->GetBallMatrix());
	ptrSynthesizer->RenderOutput();
	glPopMatrix();
	SaveSnapshot();
	glutSwapBuffers();
	if ( ptrSynthesizer->GetStepCount() == CSynConfigBase::m_stepCountMax ) ptrSynthesizer->RestartSynthesis();
}

void IdleFunc()
{
	if ( pause == false )
	{
		ptrSynthesizer->UpdateOutput();
		cout << ptrSynthesizer->GetStepCount() << endl;
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
	//gluPerspective(45, float(width)/float(height), 0.01, 1000);
	//gluOrtho2D(-3, 3, -2, 2);
	glOrtho(-3, 3, -2, 2, 0.01, 1000);
	glMatrixMode(GL_MODELVIEW);
}

void KeyboardFunc(unsigned char key, int x, int y) 
{
	switch(key) 
	{
	case 'a':
		ptrCamera->MoveLeft(0.2f);
		break;
	case 'd':
		ptrCamera->MoveRight(0.2f);
		break;
	case 'w':
		ptrCamera->MoveForward(0.2f);
		break;
	case 's':
		ptrCamera->MoveBackward(0.2f);
		break;
	case 'r':
		ptrCamera->MoveUp(0.2f);
		break;
	case 'f':
		ptrCamera->MoveDown(0.2f);
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
	case 'u':
		DumpCameraAndArcBall();
		break;
	case 'l':
		LoadCameraAndArcBall();
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

void MouseFunc(int button, int state, int x, int y)
{
	x = g_width  - x;
	y = g_height - y;
	POINT tPoint;
	tPoint.x = x;
	tPoint.y = y;
	if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
	{
		ptrArcBall->MouseDown(tPoint, g_width, g_height);
	}
	else if ( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
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

void Initialize(const char* config_file_path)
{
	ptrArcBall = new CArcBall;
	ptrArcBall->InitBall();
	ptrCamera = new CCamera(0, 0, 4);
	ptrCamera->RotateY(180);
	ptrSynthesizer = new CParticleSystemSyn(config_file_path);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA); // set display mode
	glutInitWindowSize(g_width, g_height); // set window size
	glutInitWindowPosition(0, 0); // set window position on screen
	glutCreateWindow("Group motion of soft bodies"); // set window title

	glutMouseFunc(MouseFunc);
	glutMotionFunc(MouseMoveFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutDisplayFunc(DisplayFunc);
	glutReshapeFunc(ReshapeFunc);
	glutIdleFunc(IdleFunc);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH_HINT);
	glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	LoadCameraAndArcBall();
}

void ReleasePtr()
{
	DELETE_OBJECT(ptrArcBall);
	DELETE_OBJECT(ptrCamera);
}

void usage(char** argv)
{
	std::cout << "Usage: " << argv[0] << " config_file_path\n";
}

int main(int argc, char **argv)
{
	if ( argc < 2 )
	{
		usage(argv);
		return -1;
	}
	glutInit(&argc, argv);
	Initialize(argv[1]);
	glutMainLoop();
	ReleasePtr();
	return 0;
}
