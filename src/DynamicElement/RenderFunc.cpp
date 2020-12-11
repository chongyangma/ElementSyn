#include "RenderFunc.h"

namespace det
{
void RenderCylinder(float radius_bottom, float radius_top, float height, int slices)
{
    const float delta_theta = TWO_PI / float(slices);
    const float half_height = height * 0.5f;
    for (int i = 0; i < slices; i++)
    {
        float theta = float(i) * delta_theta;
        float next_theta = float(i + 1) * delta_theta;
        glBegin(GL_TRIANGLE_STRIP);
        glVertex3f(0.0, 0.0, half_height);
        glVertex3f(radius_top * cos(theta), radius_top * sin(theta), half_height);
        glVertex3f(radius_top * cos(next_theta), radius_top * sin(next_theta), half_height);
        glVertex3f(radius_bottom * cos(next_theta), radius_bottom * sin(next_theta), -half_height);
        glVertex3f(radius_bottom * cos(theta), radius_bottom * sin(theta), -half_height);
        glVertex3f(0.0, 0.0, -half_height);
        glEnd();
    }
}
} // namespace det
