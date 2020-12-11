#include "Camera.h"

CCamera::CCamera(float x, float y, float z) :
    m_theta(90), m_phi(0), m_dirty(false), m_dir(0, -1, -1), m_pos(x, y, z)
{
}

void CCamera::MoveForward(float dist)
{
    CalcDir();
    m_pos += m_dir * dist;
}

void CCamera::MoveBackward(float dist)
{
    CalcDir();
    m_pos -= m_dir * dist;
}

void CCamera::MoveLeft(float dist)
{
    CalcDir();
    m_pos -= m_right * dist;
}

void CCamera::MoveRight(float dist)
{
    CalcDir();
    m_pos += m_right * dist;
}

void CCamera::MoveUp(float dist)
{
    CalcDir();
    m_pos += m_up * dist;
}

void CCamera::MoveDown(float dist)
{
    CalcDir();
    m_pos -= m_up * dist;
}

void CCamera::RotateY(float degrees)
{
    m_dirty = true;
    m_theta += degrees;
}

void CCamera::RotateX(float degrees)
{
    m_dirty = true;
    m_phi += degrees;
    m_phi = max(-85.0f, m_phi);
    m_phi = min(85.0f, m_phi);
}

Vec3f CCamera::GetDir()
{
    CalcDir();
    return m_dir;
}

void CCamera::Draw()
{
    CalcDir();
    Vec3f at = m_pos + m_dir * 100.f;
    gluLookAt(m_pos[0], m_pos[1], m_pos[2], at[0], at[1], at[2], m_up[0], m_up[1], m_up[2]);
}

void CCamera::DumpCamera(FILE* file)
{
    fprintf(file, "%s\t%f\n", "CAMERA_THETA", m_theta);
    fprintf(file, "%s\t%f\n", "CAMERA_PHI", m_phi);
    fprintf(file, "%s\t%f %f %f\n", "CAMERA_DIR", m_dir[0], m_dir[1], m_dir[2]);
    fprintf(file, "%s\t%f %f %f\n", "CAMERA_POS", m_pos[0], m_pos[1], m_pos[2]);
}

bool CCamera::LoadCamera(const string& fileName)
{
    ifstream fin(fileName.c_str());
    if (fin.fail() == true)
    {
        cout << "Failed to load camera parameters from config file " << fileName << "!\n";
        return false;
    }
    string param;
    while (fin >> param)
    {
        if (param == string("CAMERA_THETA"))
        {
            fin >> m_theta;
        }
        else if (param == string("CAMERA_PHI"))
        {
            fin >> m_phi;
        }
        else if (param == string("CAMERA_DIR"))
        {
            fin >> m_dir[0] >> m_dir[1] >> m_dir[2];
        }
        else if (param == string("CAMERA_POS"))
        {
            fin >> m_pos[0] >> m_pos[1] >> m_pos[2];
        }
    }
    m_dirty = true;
    return true;
}

Vec3f CCamera::GetPos()
{
    CalcDir();
    return m_pos;
}

void CCamera::CalcDir()
{
    if (!m_dirty) return;

    float rphi = Radians(m_phi);
    float rtheta = Radians(m_theta);
    float r = cos(rphi);
    m_dir[0] = r * cos(rtheta);
    m_dir[1] = sin(rphi);
    m_dir[2] = r * sin(rtheta);
    Vec3f vup(0.f, 1.f, 0.f);
    m_right = cross(m_dir, vup);
    normalize(m_dir);
    normalize(m_right);
    m_up = cross(m_right, m_dir);
    normalize(m_up);
    m_dirty = false;
}
