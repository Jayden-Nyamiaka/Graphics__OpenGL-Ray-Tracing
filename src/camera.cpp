#include "camera.h"

#include "glinclude.h"

void Frustum::OpenGLSetMatrix() const {
    gluPerspective(fov, aspect_ratio, near, far);
}

void Camera::OpenGLSetPosition() const {
    translate.OpenGLTransform();
    rotate.OpenGLTransform();
}