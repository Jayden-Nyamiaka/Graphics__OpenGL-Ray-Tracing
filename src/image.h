#ifndef IMAGE_H
#define IMAGE_H

#include <vector>

#include <Eigen/Dense>

class Image {
public:
    int xres, yres;
    std::vector<Eigen::Vector3f> pixels;

    Image();
    Image(int width, int height);

    void SetPixel(int x, int y, const Eigen::Vector3f &color);
    bool SaveImage(const std::string &path);
};

#endif