#ifndef CAMERA_H
#define CAMERA_H

#include "main.h"
#include "vec.h"

class CCamera
{
public:
    CCamera(float x, float y, float z);
    void MoveForward(float dist);
    void MoveBackward(float dist);
    void MoveLeft(float dist);
    void MoveRight(float dist);
    void MoveUp(float dist);
    void MoveDown(float dist);

    void RotateY(float degrees);
    void RotateX(float degrees);

    Vec3f GetDir();
    Vec3f GetPos();

    void Draw();

    void DumpCamera(FILE* file);
    bool LoadCamera(const string& fileName);

private:
    void CalcDir();
    float m_theta, m_phi;
    Vec3f m_pos;
    Vec3f m_dir, m_up, m_right;
    bool m_dirty;
};

#endif // CAMERA_H
