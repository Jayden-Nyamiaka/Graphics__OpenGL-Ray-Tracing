#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Dense>

class Color {
public:
    float r;
    float g;
    float b;

    Color(): r(0.0f), g(0.0f), b(0.0f) {};
    Color(float r, float g, float b): r(r), g(g), b(b) {};

    Eigen::Vector3f ToVector() const;
};

class Arcball {
private:
    bool enabled;
    int x_max, y_max;
    Eigen::Quaternionf base;
    Eigen::Quaternionf curr;
    Eigen::Vector3f start;
public:
    Arcball() = default;
    Arcball(int x, int y);
    void SetRes(int x, int y);
    void Start(int x, int y);
    void Move(int x, int y);
    void End();
    void Apply();
};

#endif