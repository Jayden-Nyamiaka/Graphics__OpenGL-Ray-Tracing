#include "object.h"
#include "scene.h"

#include <iostream>

#include "image.h"

using namespace Eigen;
using namespace std;

const int MAX_ITERS = 10000;
const int XRES = 500;
const int YRES = 500;


// Helper functions
double inside_outside_func(double x, double y, double z, double exp, double n) {
    return -1.0 + pow(z*z, 1.0/n) + pow( pow(x*x, 1.0/exp) + pow(y*y, 1.0/exp) , exp/n);
}

Vector3d gradient_inside_outside_func(double x, double y, double z, double exp, double n) {
    double derivative_x = 2.0 * x * pow(x*x, 1.0/exp-1.0) * pow( pow(x*x, 1.0/exp) + pow(y*y, 1.0/exp), exp/n-1.0);
    double derivative_y = 2.0 * y * pow(y*y, 1.0/exp-1.0) * pow( pow(x*x, 1.0/exp) + pow(y*y, 1.0/exp), exp/n-1.0);
    double derivative_z = 2.0 * z * pow(z*z, 1.0/n-1.0);
    return (1.0 / n) * Vector3d(derivative_x, derivative_y, derivative_z);
}

static const double CLOSE_ENOUGH_BOUND = 1.0 / 20.0;
bool close_enough(double x) {
    return (abs(x) < CLOSE_ENOUGH_BOUND);
}

int sign(double x) {
    return (x < 0) ? -1 : 1;
}

/**
 * IOTest Code
 */


// PART 1.1
bool Superquadric::IOTest(const Vector3d &point) {
    // Transforms the point from world space to body coordinates
    /* Explanation: point is given in world space, thus we have big X.
     * The inside outside function takes little x, so we need to get little x.
     * big X = Transform Matrix O on (x) = T(R(S(x))) where TRS are transforms applied in order S, R, T.
     * So Inverse Transform Matrix O on (X) = x = S^-1(R^- 1(T^-1(big X))).
     * Notice the 1st object transform applied is the most right in the matrix multiplication.
     * That's what's happening here. */
    Vector4d vec_point(point[0], point[1], point[2], 1.0);
    Vector4d body_point = getInverseTransformMatrix() * vec_point;

    // Compute Superquadric inside-outside function using transformed body coordinate
    return (inside_outside_func(body_point[0], body_point[1], body_point[2], exp0, exp1) < 0);
}


// PART 1.2
bool Assembly::IOTest(const Vector3d &point) {
    // Transforms the point from world space to assembly coordinates
    Vector4d vec_point(point[0], point[1], point[2], 1.0);
    Vector3d assembly_point = (getInverseTransformMatrix() * vec_point).head<3>();

    /* Recursively calls IOTest on all children assembly and primitive objects
     * using tranformed assembly coordinate */
    for (size_t i = 0; i < children.size(); i++) {

        bool inside = children[i]->IOTest(assembly_point);

        // Returns true if the point is inside any of its children objects
        if (inside) {
            return true;
        }
    }

    // Returns false if none of the children objects have returned true
    return false;
}



/**
 * Closest Intersection Code
 */


// Part 1.3
pair<double, Intersection> Superquadric::ClosestIntersection(const Ray &ray) {
    // Initialize the return pair to no intersection
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());

    // Transform ray = a * t + b parameters to body coordinates
    Ray transformed_ray = ray;
    transformed_ray.Transform(getInverseTransformMatrix());
    transformed_ray.Normalize();
    Vector3d vec_a = transformed_ray.origin;
    Vector3d vec_b = transformed_ray.direction;
    /*Matrix4d worldToBodySpace = getInverseTransformMatrix();
    Vector4d vec_origin(ray.origin[0], ray.origin[1], ray.origin[2], 1.0);
    Vector4d vec_direction(ray.direction[0], ray.direction[1], ray.direction[2], 1.0);
    Vector3d vec_a = (worldToBodySpace * vec_origin).head<3>();
    Vector3d vec_b = (worldToBodySpace * vec_direction).head<3>();
    */ // maybe normalize ray too

    // Calculates initial value of t
    double a = vec_a.dot(vec_a);
    double b = 2.0 * vec_a.dot(vec_b);
    double c = vec_b.dot(vec_b) - 3.0;
    double discriminant = b*b - 4.0*a*c;
    // No collision: if discriminant is negative, solutions for t are imaginary
    if (discriminant < 0) {
        return closest;
    }
    double t_1 = (-b - sign(b) * sqrt(discriminant)) / (2.0* a);
    double t_2 = (2.0 * c) / (-b - sign(b) * sqrt(discriminant));
    // No collision: if both t are negative, the obj is behing the camera
    if (t_1 < 0 && t_2 < 0) {
        return closest;
    }
    double t_pos, t_neg;
    if (b < 0) {
        t_pos = t_1;
        t_neg = t_2;
    } else {
        t_pos = t_2;
        t_neg = t_1;
    }
    double t;
    if (t_neg >= 0) {
        t = t_neg;
    } else {
        t = t_pos;
    }

    // Iteratively traces the ray, finding point of intersection if it exists
    do {
        // Note: the ray equation gives us our location at t
        Vector3d loc = vec_a * t + vec_b;

        // Returns an intersection if io function is approx 0 (means we're on surface of obj)
        double io_value = inside_outside_func(loc[0], loc[1], loc[2], exp0, exp1);
        if (close_enough(io_value)) {
            closest.first = t;
            closest.second.location = ray;
            closest.second.obj = this;
        }

        // Stopping Condition: if the derivative isn't negative, we've missed the obj
        double io_derivative = vec_a.dot(gradient_inside_outside_func(loc[0], loc[1], loc[2], exp0, exp1));
        if (io_derivative >= 0) {
            break;
        }

        // Updates t using Newton's Method (first-order Taylor Series approximation)
        t -= io_value / io_derivative;
    } while(true);

    return closest;
}


// Part 1.4
pair<double, Intersection> Assembly::ClosestIntersection(const Ray &ray) {
    // Initialize the return pair to no intersection
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());

    // Transform ray to assembly coordinates
    Ray transformed_ray = ray;
    transformed_ray.Transform(getInverseTransformMatrix());
    transformed_ray.Normalize();

    /* Recursively calls closest intersection on all children objs, finding 
     * the intersection with the smallest non-negative value of t */
    for (size_t i = 0; i < children.size(); i++) {
        pair<double, Intersection> collision = children[i]->ClosestIntersection(transformed_ray);
        if (collision.first < closest.first) {
            closest = collision;
        }
    }
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
