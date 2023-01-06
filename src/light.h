#ifndef LIGHT_H
#define LIGHT_H

#include <Eigen/Dense>

#include "util.h"

class Light {
private:
public:
    Eigen::Vector4d position;

    Color color;
    float attenuation;
    
    Light(): position(Eigen::Vector4d::Zero()), color(Color()), attenuation(0.0) {};
    Light(Eigen::Vector4d position, Color color, float attenuation): position(position), color(color), attenuation(attenuation) {};

    void OpenGLInit(int light_id);
    void OpenGLRender(int light_id);
};

#endif // LIGHT_H