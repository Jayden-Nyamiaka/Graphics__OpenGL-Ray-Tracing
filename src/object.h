#ifndef OBJECT_H
#define OBJECT_H

#include <memory>
#include <experimental/optional>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "transform.h"
#include "util.h"

// Some forward declarations...
class Ray;
class Intersection;

class Material {
public:
    Color ambient;
    Color diffuse;
    Color specular;
    
    float shininess;
    float reflected;
    float refracted;

    Material();

    void SetOpenGLMaterial();
};

class Object {
protected:
    std::vector<std::unique_ptr<Transformation>> transforms;
public:
    // Disallow copying.
    Object(const Object&) = delete;

    Object(): transforms(std::vector<std::unique_ptr<Transformation>>()) {};
    
    // Extension of transforms
    Eigen::Matrix4d getForwardTransformMatrix();
    Eigen::Matrix4d getInverseTransformMatrix();

    virtual void Tesselate(std::vector<Eigen::Vector3f> &vertices, std::vector<Eigen::Vector3f> &normals) = 0;
    virtual void OpenGLRender();
    
    virtual bool IOTest(const Eigen::Vector3d &point) = 0;
    virtual std::pair<double, Intersection> ClosestIntersection(const Ray &ray) = 0;

    template <class T>
    void AddTransform(const T &t) {
        transforms.push_back(std::make_unique<T>(t));
    }
};

class Superquadric: public Object {
private:
    double exp0;
    double exp1;
    float patch_u;
    float patch_v;
    Material mat;

    // Stores the location in the OpenGL vertex and normal buffers.
    size_t buffer_start;
    size_t buffer_end;
public:
    // Disallow copying.
    Superquadric(const Superquadric&) = delete;

    Superquadric();
    Superquadric(double e0, double e1);
    void Tesselate(std::vector<Eigen::Vector3f> &vertices, std::vector<Eigen::Vector3f> &normals);
    void OpenGLRender();

    void SetMaterial(const Material &material);
    const Material& GetMaterial() const;

    bool IOTest(const Eigen::Vector3d &point);
    std::pair<double, Intersection> ClosestIntersection(const Ray &ray);

    Eigen::Vector3f GetVertex(float u, float v);
    Eigen::Vector3d GetNormal(const Eigen::Vector3d &vertex);
};

class Assembly: public Object {
private:
    std::vector<std::shared_ptr<Object>> children;
public:
    // Disallow copying.
    Assembly(const Assembly&) = delete;

    Assembly();
    void OpenGLRender();
    void Tesselate(std::vector<Eigen::Vector3f> &vertices, std::vector<Eigen::Vector3f> &normals);

    bool IOTest(const Eigen::Vector3d &point);
    std::pair<double, Intersection> ClosestIntersection(const Ray &ray);

    void AddChild(std::shared_ptr<Object> obj) {
        children.push_back(obj);
    }
};

// Define ray and intersection now that we have the needed types.

class Ray {
public:
    Eigen::Vector3d origin;
    Eigen::Vector3d direction;

    Ray() = default;

    Eigen::Vector3d At(float time);
    Ray Transformed(const Eigen::Matrix4d &transform) const;
    void Transform(const Eigen::Matrix4d &transform);
    float Normalize();
};

class Intersection {
public:
    Ray location;
    Superquadric *obj;

    Intersection();
    Intersection(Ray location, Superquadric *obj);
};


#endif // OBJECT_H