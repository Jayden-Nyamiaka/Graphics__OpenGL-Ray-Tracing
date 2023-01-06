#include "light.h"

#include "glinclude.h"

void Light::OpenGLInit(int light_id) {
    glEnable(light_id);

    glLightfv(light_id, GL_AMBIENT, (float *) &color);
    glLightfv(light_id, GL_DIFFUSE, (float *) &color);
    glLightfv(light_id, GL_SPECULAR, (float *) &color);

    glLightf(light_id, GL_QUADRATIC_ATTENUATION, attenuation);
}

void Light::OpenGLRender(int light_id) {
    Eigen::Vector4f pos = position.cast<float>();
    glLightfv(light_id, GL_POSITION, pos.data());

    float max_color = std::max(std::max(color.r, color.g), color.b);
    float lcolor[3] = {
        color.r / max_color,
        color.g / max_color,
        color.b / max_color
    };

    glPointSize(8.0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, lcolor);
    glEnable(GL_POINT_SMOOTH);

    glBegin(GL_POINTS);
    glVertex4dv(position.data());
    glEnd();
}