#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Eigen/Dense>

class Transformation {
public:
    virtual void OpenGLTransform() const = 0;
    virtual Eigen::Matrix4d GetMatrix() = 0;
};

class Rotate: public Transformation {
private:
    double axis_x;
    double axis_y;
    double axis_z;
    double angle;
public:
    Rotate() : axis_x(1), axis_y(0), axis_z(0), angle(0) {};
    Rotate(double axis_x, double axis_y, double axis_z, double angle);
    Eigen::Matrix4d GetMatrix();
    void OpenGLTransform() const;
    double GetAngle() const;
    Eigen::Vector3d GetAxis() const;
    void SetAxisAngle(Eigen::Vector3d axis, double theta);
};

class Translate: public Transformation {
private:
    double delta_x;
    double delta_y;
    double delta_z;
public:
    Translate() : delta_x(0), delta_y(0), delta_z(0) {};
    Translate(double delta_x, double delta_y, double delta_z) : delta_x(delta_x), delta_y(delta_y), delta_z(delta_z) {};
    Eigen::Matrix4d GetMatrix();
    void OpenGLTransform() const;
    Eigen::Vector3d GetDelta() const;
    void SetDelta(Eigen::Vector3d delta);
};

class Scale: public Transformation {
private:
    double scale_x;
    double scale_y;
    double scale_z;
public:
    Scale() : scale_x(1), scale_y(1), scale_z(1) {};
    Scale(double scale_x, double scale_y, double scale_z) : scale_x(scale_x), scale_y(scale_y), scale_z(scale_z) {};
    Eigen::Matrix4d GetMatrix();
    void OpenGLTransform() const;
    Eigen::Vector3d GetScale() const;
    void SetScale(Eigen::Vector3d scale);
};

#endif // TRANSFORM_H