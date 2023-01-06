#include "util.h"

#include <iostream>

#include "glinclude.h"

using namespace Eigen;

/**
 * Color implementation
 */

Eigen::Vector3f Color::ToVector() const {
    return Vector3f { r, g, b };
}

/**
 * Arcball implementation
 */

Arcball::Arcball(int x, int y) {
    x_max = x;
    y_max = y;
    enabled = false;
    start = Vector3f::Zero();
    base = Quaternionf::Identity();
    curr = Quaternionf::Identity();
}

void Arcball::SetRes(int x, int y) {
    x_max = x;
    y_max = y;
}

void Arcball::Start(int x, int y) {
    enabled = true;

    float x_pos = ((float)x / x_max) * 2 - 1;
    float y_pos = -((float)y / y_max) * 2 + 1;
    float z_pos = sqrt(1 - std::min(1.0f, x_pos * x_pos + y_pos * y_pos));

    start = { x_pos, y_pos, z_pos };
}

void Arcball::Move(int x, int y) {
    if (!enabled)
        return;

    float x_pos = ((float)x / x_max) * 2 - 1;
    float y_pos = -((float)y / y_max) * 2 + 1;
    float z_pos = sqrt(1 - std::min(1.0f, x_pos * x_pos + y_pos * y_pos));

    Vector3f dir = { x_pos, y_pos, z_pos };

    curr = Quaternionf::FromTwoVectors(start, dir);
}

void Arcball::End() {
    enabled = false;

    base = (curr * base).normalized();
    curr = Quaternionf::Identity();
}

void Arcball::Apply() {
    Quaternionf rotation = enabled ? (curr * base).normalized() : base;

    float s = rotation.w();
    float theta = 2 * acos(s) * 180 / M_PI;
    Vector3f axis = rotation.coeffs().head<3>() / sqrt(1 - s * s);

    glRotatef(theta, axis.x(), axis.y(), axis.z());
}
