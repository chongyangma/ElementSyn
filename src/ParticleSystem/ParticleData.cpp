#include "ParticleData.h"

CParticleData::CParticleData()
{
    m_pos = Vec3f(0.0f, 0.0f, 0.0f);
    m_vel = Vec3f(0.0f, 0.0f, 0.0f);
    m_quat = Vec4f(1.0f, 0.0f, 0.0f, 0.0f);
    m_flagBoundary = false;
    m_flagFixed = false;
}

void CParticleData::RenderSoftBody(Vec3f trans /* = Vec3f */)
{
    Vec3f ambient = CParticleSystemConfig::m_particleAmbient;
    Vec3f diffuse = CParticleSystemConfig::m_particleDiffuse;
    Vec3f specular = CParticleSystemConfig::m_particleSpecular;
    if (GetFlagFixed() == true)
    {
        diffuse = Vec3f(1.f, 0.f, 0.f);
    }
    GLfloat objAmbient[] = { ambient[0], ambient[1], ambient[2], 1.0f };
    GLfloat objDiffuse[] = { diffuse[0], diffuse[1], diffuse[2], 1.0f };
    GLfloat objSpecular[] = { specular[0], specular[1], specular[2], 1.0f };
    glMaterialfv(GL_FRONT, GL_AMBIENT, objAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, objDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, objSpecular);
    glPushMatrix();
    glTranslatef(trans[0], trans[1], trans[2]);
    glutSolidSphere(CParticleSystemConfig::m_repulsionDist, 20, 20);
    glPopMatrix();
}

void CParticleData::TranslateSoftBody(Vec3f trans)
{
    for (int i = 0; i < int(m_vecSamplePos.size()); i++)
    {
        m_vecSamplePos[i] += trans;
    }
}
