#include "transform.h"

#include <iostream>

#include "glinclude.h"

using namespace Eigen;

/**
 * Rotate Transform Implementation.
 */

Rotate::Rotate(double axis_x, double axis_y, double axis_z, double angle) {
    double norm = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
    this->axis_x = axis_x / norm;
    this->axis_y = axis_y / norm;
    this->axis_z = axis_z / norm;
}

void Rotate::OpenGLTransform() const {
    glRotatef(angle, axis_x, axis_y, axis_z);
}

Matrix4d Rotate::GetMatrix() {
    double theta = angle * M_PI / 180;
    double cos1m = 1 - cos(theta);

    Matrix4d rotation;

    rotation <<
        // Row 1
        (axis_x * axis_x) + (1 - (axis_x * axis_x)) * cos(theta),
        (axis_x * axis_y) * cos1m - axis_z * sin(theta),
        (axis_x * axis_z) * cos1m + axis_y * sin(theta),
        0,
        // Row 2
        (axis_y * axis_x) * cos1m + axis_z * sin(theta),
        (axis_y * axis_y) + (1 - (axis_y * axis_y)) * cos(theta),
        (axis_y * axis_z) * cos1m - axis_x * sin(theta),
        0,
        // Row 3
        (axis_z * axis_x) * cos1m - axis_y * sin(theta),
        (axis_z * axis_y) * cos1m + axis_x * sin(theta),
        (axis_z * axis_z) + (1 - (axis_z * axis_z)) * cos(theta), 0,
        // Row 4
        0, 0, 0, 1;
    
    return rotation;
}

double Rotate::GetAngle() const {
    return angle;
}

Eigen::Vector3d Rotate::GetAxis() const {
    return Vector3d { axis_x, axis_y, axis_z };
}

void Rotate::SetAxisAngle(Eigen::Vector3d axis, double theta) {
    axis = axis.normalized();
    axis_x = axis[0];
    axis_y = axis[1];
    axis_z = axis[2];
    angle = theta;
}

/**
 * Translate Transformation Implementation.
 */


void Translate::OpenGLTransform() const {
    glTranslatef(delta_x, delta_y, delta_z);
}

Matrix4d Translate::GetMatrix() {
    Matrix4d translation;
    translation << 1, 0, 0, delta_x,
                   0, 1, 0, delta_y,
                   0, 0, 1, delta_z,
                   0, 0, 0, 1;

    return translation;
}

Eigen::Vector3d Translate::GetDelta() const {
    return Vector3d { delta_x, delta_y, delta_z };
}

void Translate::SetDelta(Eigen::Vector3d delta) {
    delta_x = delta[0];
    delta_y = delta[1];
    delta_z = delta[2];
}


/**
 * Scale Transformation Implementation.
 */


void Scale::OpenGLTransform() const {
    glScalef(scale_x, scale_y, scale_z);
}

Matrix4d Scale::GetMatrix() {
    Matrix4d scaling;
    scaling << scale_x, 0, 0, 0,
                0, scale_y, 0, 0,
                0, 0, scale_z, 0,
                0, 0, 0,  1;

    return scaling;
}

Eigen::Vector3d Scale::GetScale() const {
    return Vector3d { scale_x, scale_y, scale_z };
}

void Scale::SetScale(Eigen::Vector3d scale) {
    scale_x = scale[0];
    scale_y = scale[1];
    scale_z = scale[2];
}
