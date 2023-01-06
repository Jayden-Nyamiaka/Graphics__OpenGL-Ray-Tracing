#include "object.h"

#include <cmath>
#include <iostream>

#include "glinclude.h"

using namespace std;
using namespace Eigen;

inline int sign(float x) {
    return (x > 0) ? 1 : -1;
}

/* Calculates the parametric superquadric sine analog */
inline float pSin(float u, float p) {
    float sin_u = sinf(u);
    return sign(sin_u) * powf(fabs(sin_u), p);
}

/* Calculates the parametric superquadric cosine analog */
inline float pCos(float u, float p) {
    float cos_u = cosf(u);
    return sign(cos_u) * powf(fabs(cos_u), p);
}

/**
 * Material Implementation
 */


Material::Material() {
    ambient = Color(0.2, 0.2, 0.2);
    diffuse = Color(0.6, 0.6, 0.6);
    specular = Color(0.4, 0.4, 0.4);

    shininess = 10;
    reflected = 0;
    refracted = 0;
}

void Material::SetOpenGLMaterial() {
    glMaterialfv(GL_FRONT, GL_AMBIENT, (float *) &ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, (float *) &diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, (float *) &specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
}

/**
 * Base Object Implementation
 */


void Object::OpenGLRender() {
    for (auto it = transforms.rbegin(); it != transforms.rend(); it++) {
        (*it)->OpenGLTransform();
    }
}

/**
 * Superquadric Implementation
 */


Superquadric::Superquadric() {
    Superquadric(1.0, 1.0);
}

Superquadric::Superquadric(double e0, double e1) {
    exp0 = e0;
    exp1 = e1;
    patch_u = 15;
    patch_v = 15;

    mat = Material();
}

void Superquadric::SetMaterial(const Material &material) {
    mat = material;
}

const Material& Superquadric::GetMaterial() const {
    return mat;
}

Vector3f Superquadric::GetVertex(float u, float v) {
    float cos_v = pCos(v, exp1);
    float sin_v = pSin(v, exp1);

    Vector3f vert;
    vert << cos_v * pCos(u, exp0), cos_v * pSin(u, exp0), sin_v;

    return vert;
}

Vector3d Superquadric::GetNormal(const Vector3d &vertex) {
    double x = vertex(0);
    double y = vertex(1);
    double z = vertex(2);

    Vector3d normal;
    normal(0) = (x == 0.0) ? 0.0 :
        2.0 * x * pow(x * x, 1.0 / this->exp0 - 1.0) *
        pow(pow(x * x, 1.0 / this->exp0) + pow(y * y, 1.0 / this->exp0),
            this->exp0 / this->exp1 - 1.0) / this->exp1;
    normal(1) = (y == 0.0) ? 0.0 :
        2.0 * y * pow(y * y, 1.0 / this->exp0 - 1.0) *
        pow(pow(x * x, 1.0 / this->exp0) + pow(y * y, 1.0 / this->exp0),
            this->exp0 / this->exp1 - 1.0) / this->exp1;
    normal(2) = (z == 0.0) ? 0.0 :
        2.0 * z * pow(z * z, 1.0 / this->exp1 - 1.0) / this->exp1;

    normal.normalize();

    return normal;
}

void Superquadric::Tesselate(std::vector<Eigen::Vector3f> &vertices, std::vector<Eigen::Vector3f> &normals) {
    buffer_start = vertices.size();
    
    float half_pi = M_PI / 2;
    float du = 2 * M_PI / patch_u, dv = M_PI / patch_v;

    auto add_vert = [&](float u, float v) mutable {
        Vector3f vertex = GetVertex(u, v);
        Vector3f normal = GetNormal(vertex.cast<double>()).cast<float>();
        vertices.push_back(vertex);
        normals.push_back(normal);
    };

    float u, v;

    // Tesselate the latitude ranges except top and bottom.
    // Uses GL_TRIANGLE_STRIPs
    for (v = dv - half_pi; v < half_pi - dv; v += dv) {
        for (u = -M_PI; u <= M_PI; u += du) {
            add_vert(u, v + dv);
            add_vert(u, v);
        }
        add_vert(-M_PI, v + dv);
        add_vert(-M_PI, v);
    }

    // Tesselate the bottom ring.
    // Uses GL_TRIANGLE_FAN
    add_vert(M_PI, -half_pi);
    v = dv - half_pi;
    for (u = M_PI; u > -M_PI; u -= du) {
        add_vert(u, v);
    }
    add_vert(-M_PI, v);

    // Tesselate the top ring.
    // Uses GL_TRIANGLE_FAN
    add_vert(-M_PI, half_pi);
    v *= -1;
    for (u = -M_PI; u < M_PI; u += du) {
        add_vert(u, v);
    }
    add_vert(M_PI, v);

    buffer_end = vertices.size();
}

void Superquadric::OpenGLRender() {
    glPushMatrix();

    // Apply transforms.
    Object::OpenGLRender();

    // Render the superquadric.
    size_t offset = 2 * (patch_u + 1);
    size_t start = buffer_start;

    mat.SetOpenGLMaterial();

    for (int j = 1; j < patch_v - 1; j++) {
        glDrawArrays(GL_TRIANGLE_STRIP, start, offset);
        start += offset;
    }

    offset = patch_u + 2;
    glDrawArrays(GL_TRIANGLE_FAN, start, offset);

    start += offset;
    glDrawArrays(GL_TRIANGLE_FAN, start, offset);

    glPopMatrix();
}


/**
 * Assembly Implementation
 */

Assembly::Assembly() {}

void Assembly::OpenGLRender() {
    glPushMatrix();
    Object::OpenGLRender();

    for (auto &child : children) {
        child->OpenGLRender();
    }

    glPopMatrix();
}

void Assembly::Tesselate(std::vector<Eigen::Vector3f> &vertices, std::vector<Eigen::Vector3f> &normals) {
    for (auto &child : children) {
        child->Tesselate(vertices, normals);
    }
}

/**
 * Ray implementation
 */

Vector3d Ray::At(float time) {
    return origin + time * direction;
}

Ray Ray::Transformed(const Matrix4d &transform) const {
    Vector4d new_origin = Vector4d::Ones();
    new_origin.head(3) = origin;
    Vector3d new_direction = Vector3d::Ones();

    new_origin = transform * new_origin;
    new_direction = transform.block<3, 3>(0, 0) * direction;

    return { new_origin.head(3), new_direction };
}

void Ray::Transform(const Matrix4d &transform) {
    Vector4d new_origin = Vector4d::Ones();
    new_origin.head(3) = origin;

    new_origin = transform * new_origin;
    direction = transform.block<3, 3>(0, 0) * direction;

    origin = new_origin.head(3);
}

float Ray::Normalize() {
    float norm = direction.norm();
    direction /= norm;
    return norm;
}

Intersection::Intersection() {
    location = Ray();
    obj = nullptr;
}

Intersection::Intersection(Ray loc, Superquadric *sq) {
    location = loc;
    obj = sq;
}