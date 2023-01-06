#ifndef CAMERA_H
#define CAMERA_H

#include "transform.h"

class Frustum {
public:
    double aspect_ratio;
    double fov;
    double near;
    double far;

    void OpenGLSetMatrix() const;
};

class Camera {
public:
    Translate translate;
    Rotate rotate;
    Frustum frustum;

    void OpenGLSetPosition() const;
};

#endif // CAMERA_H