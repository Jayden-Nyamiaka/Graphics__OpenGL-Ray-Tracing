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

bool Superquadric::IOTest(const Vector3d &point) {
    /**
     * PART 1
     * TODO: Implement the IO Test function for a superquadric. Remember to
     *       apply any transformations to the superquadric before testing
     *       against the IO function.
     */
    return false;
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
