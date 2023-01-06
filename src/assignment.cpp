#include "object.h"
#include "scene.h"

#include <iostream>

#include "image.h"

using namespace Eigen;
using namespace std;

const int MAX_ITERS = 10000;
const int XRES = 500;
const int YRES = 500;

/**
 * IOTest Code
 */

// PART 1
bool Superquadric::IOTest(const Vector3d &point) {
    // Transforms the point from world space to body coordinates
    /* Explanation: point is given in world space, thus we have big X.
     * The inside outside function takes little x, so we need to get little x.
     * big X = Transform Matrix O on (x) = T(R(S(x))) where TRS are transforms applied in order S, R, T.
     * So Inverse Transform Matrix O on (X) = x = S^-1(R^- 1(T^-1(big X))).
     * Notice the 1st object transform applied is the most right in the matrix multiplication.
     * That's what's happening here. */
    Matrix4d worldToBodySpace = Matrix4d::Identity();
    for (int i = 0; i < transforms.size(); i++) {
        Matrix4d transform = transforms[i]->GetMatrix();
        worldToBodySpace = worldToBodySpace * transform.inverse();
    }
    Vector4d body_point = worldToBodySpace * Vector4d(point, 1.0);

    // Uses transformed body coordinate to compute Superquadric inside-outside function
    double x = body_point[0];
    double y = body_point[1];
    double z = body_point[2];
    double in_out = -1.0 + pow(z*z, 1.0/exp1) + pow( pow(x*x, 1.0/exp0) + pow(y*y, 1.0/exp0) , exp0/exp1):
    return (in_out < 0);
}

bool Assembly::IOTest(const Vector3d &point) {
    /**
     * PART 1
     * TODO: Implement the IO Test function for an assembly (recursively call
     *       IOTest on the children). Make sure to apply any transformations
     *       to the assembly before calling IOTest on the children.
     */
    return false;
}

/**
 * Closest Intersection Code
 */

pair<double, Intersection> Superquadric::ClosestIntersection(const Ray &ray) {
    /**
     * PART 1
     * TODO: Implement a ray-superquadric intersection using Newton's method.
     *       Make sure to apply any transformations to the superquadric before
     *       performing Newton's method.
     */
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
    return closest;
}

pair<double, Intersection> Assembly::ClosestIntersection(const Ray &ray) {
    /**
     * PART 1
     * TODO: Implement a ray-assembly intersection by recursively finding
     *       intersection with the assembly's children. Make sure to apply any
     *       transformations to the assembly before calling ClosestIntersection
     *       on the children.
     */
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
    return closest;
}

/**
 * Raytracing Code
 */

void Scene::Raytrace() {
    Image img = Image(XRES, YRES);

    for (int i = 0; i < XRES; i++) {
        for (int j = 0; j < YRES; j++) {
            /**
             * PART 2
             * TODO: Implement raytracing using the code from the first part
             *       of the assignment. Set the correct color for each pixel
             *       here.
             */
            img.SetPixel(i, j, Vector3f::Ones());
        }
    }

    // Outputs the image.
    if (!img.SaveImage("rt.png")) {
        cerr << "Error: couldn't save PNG image" << std::endl;
    } else {
        cout << "Done!\n";
    }
}
