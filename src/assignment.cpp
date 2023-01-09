#include "object.h"
#include "scene.h"

#include <iostream>

#include "image.h"

using namespace Eigen;
using namespace std;

const size_t MAX_ITERS = 10000;
const size_t XRES = 500;
const size_t YRES = 500;


// HELPER FUNCTIONS
Eigen::Matrix4d Object::getForwardTransformMatrix() {
    Eigen::Matrix4d forwardTransform = Matrix4d::Identity();
    for (size_t i = 0; i < transforms.size(); i++) {
        forwardTransform = transforms[i]->GetMatrix() * forwardTransform;
    }
    return forwardTransform;
}

// Every Object has transforms so this returns the inverse transform matrix for that
Eigen::Matrix4d Object::getInverseTransformMatrix() {
    Eigen::Matrix4d inverseTransform = Matrix4d::Identity();
    for (size_t i = 0; i < transforms.size(); i++) {
        inverseTransform *= transforms[i]->GetMatrix().inverse();
    }
    return inverseTransform;
}

double Superquadric::IOFunction(Vector3d pos) {
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    return -1.0 + pow(z*z, 1.0/exp1) + pow( pow(x*x, 1.0/exp0) + pow(y*y, 1.0/exp0) , exp0/exp1);
}

Vector3d Superquadric::IOGradient(Vector3d pos) {
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double derivative_x = 2.0 * x * pow(x*x, 1.0/exp0-1.0) * pow( pow(x*x, 1.0/exp0) + pow(y*y, 1.0/exp0), exp0/exp1-1.0);
    double derivative_y = 2.0 * y * pow(y*y, 1.0/exp0-1.0) * pow( pow(x*x, 1.0/exp0) + pow(y*y, 1.0/exp0), exp0/exp1-1.0);
    double derivative_z = 2.0 * z * pow(z*z, 1.0/exp1-1.0);
    return (1.0 / n) * Vector3d(derivative_x, derivative_y, derivative_z);
}

static double close_enough_bound = 1.0 / 20.0;
bool close_enough(double x) {
    return (abs(x) < close_enough_bound);
}

int sign(double x) {
    return (x < 0) ? -1 : 1;
}

float deg2rad(float angle)
{
    return angle * M_PI / 180.0;
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
    return (IOFunction(body_point) < 0);
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
// Note: Returned ray is not normalized
pair<double, Intersection> Superquadric::ClosestIntersection(const Ray &ray) {
    // Initialize the return pair to no intersection
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());

    // Transform ray = a * t + b parameters to body coordinates
    Ray transformed_ray = ray.Transformed(getInverseTransformMatrix());

    // Calculates initial value of t
    Vector3d vec_a = transformed_ray.origin;
    Vector3d vec_b = transformed_ray.direction;
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
        Vector3d loc = transformed_ray.At(t);

        // Returns an intersection if io function is approx 0 (means we're on surface of obj)
        double io_value = IOFunction(loc);
        if (close_enough(io_value)) {
            // Updates closest (lowest positive) value of t so far and sets returning obj
            closest.first = t;
            closest.second.obj = this;

            // Sets the returning ray's location to the point of intersection and 
            // its direction to the surface normal on the Superquadric at the point of intersection
            closest.second.location.origin = loc;
            closest.second.location.direction = IOGradient(loc);

            // Transforms the returning ray from Body Coordinates to Assembly / World Space
            closest.second.location.Transform(getForwardTransformMatrix());

            return closest;
        }

        // Stopping Condition: if the derivative isn't negative, we've missed the obj
        double io_derivative = vec_a.dot(IOGradient(loc));
        if (io_derivative >= 0) {
            break;
        }

        // Updates t using Newton's Method (first-order Taylor Series approximation)
        t -= io_value / io_derivative;
    } while(true);

    return closest;
}


// Part 1.4
// Note: Returned ray is normalized
pair<double, Intersection> Assembly::ClosestIntersection(const Ray &ray) {
    // Initialize the return pair to no intersection
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());

    // Transform ray to assembly coordinates
    Ray transformed_ray = ray.Transformed(getInverseTransformMatrix());

    /* Recursively calls closest intersection on all children objs, finding 
     * the intersection with the smallest non-negative value of t */
    for (size_t i = 0; i < children.size(); i++) {
        pair<double, Intersection> collision = children[i]->ClosestIntersection(transformed_ray);
        if (collision.first < closest.first) {
            closest = collision;
        }
    }
    
    // Transforms the closest returning ray from Assembly Space to World Space
    closest.second.location.Transform(getForwardTransformMatrix());

    // Normalized final ray
    closest.second.location.Normalize();

    return closest;
}



/**
 * Raytracing Code
 */


// Computes point lighting calculation for a point in World Space if there's a collision
Vector3f pointLighting(Ray &ray, 
                       Vector3d camera_pos, 
                       vector<Light> &lights, 
                       Superquadric *obj, 
                       vector<shared_ptr<Object>> root_objects) {
    // Gets camera direction in World Space
    Vector3d camera_dir = camera_pos - ray.origin;
    camera_dir.normalize();

    // Gets the material properties for our colliding Superquadric
    const Material &mat = obj->GetMaterial();

    // Defines Color Component Sums for Diffuse & Specular Light Reflection
    Vector3f diffuse_total = Vector3f::Zero();
    Vector3f specular_total = Vector3f::Zero();

    // Does the Shadowing and Point Lighting Calculation for every point light
    for (size_t light_idx = 0; light_idx < lights.size(); light_idx++) {
        Vector3d light_pos = lights[light_idx].position.head<3>();
        

        /* Shadowing: Tests if light ray is obstructed by another Object and 
         * only does lighting computations for this light if there's no obstruction */
        Ray light_ray = Ray();
        light_ray.origin = light_pos;
        light_ray.direction = ray_origin - light_pos;

        double t = 0;
        bool obstructed = false;
        do {
            // Note: the ray equation gives us our location at t
            Vector3d current_location = light_ray.At(t);

            // Marks light obstructed and breaks out if we're inside any other object
            for (size_t obj_idx; obj_idx < root_objects.size(); obj_idx++)c{
                if (root_objects[obj_idx]->IOTest(current_location)) {
                    obstructed = true;
                    break;
                }

            }

            if (obstructed) {
                break;
            }

            /* Updates t using Newton's Method (first-order Taylor Series approximation)
             * based on the inside outside function of the object we're trying to color */
            double io_value = obj->IOFunction(current_location);
            double io_derivative = light_ray.direction.dot(obj->IOFunction(current_location));
            t -= io_value / io_derivative;
            
        } while (t < 1.0 - close_enough_bound);

        // If the light is obstructed, skip to the next point light
        if (obstructed) {
            continue;
        }


        // Point Lighting Computation
        Vector3f light_color = lights[light_idx].color.ToVector();
        Vector3d light_dir = light_pos - ray.origin;
        light_dir.normalize();
        

        // Computes Attenuation Factor
        double x_dif = (ray.origin[0] - light_pos[0]);
        double y_dif = (ray.origin[1] - light_pos[1]);
        double z_dif = (ray.origin[2] - light_pos[2]);
        double dist_to_light_squared = x_dif*x_dif + y_dif*y_dif + z_dif*z_dif;
        double attenuation_factor = 
            1.0 / (1.0 + lights[light_idx].attenuation * dist_to_light_squared);

        // Updates Diffuse Total with Diffuse Reflection for this point light
        diffuse_total += light_color * max(0.0, attenuation_factor * ray.direction.dot(light_dir));

        // Updates Specular Total with Specular Reflection for this point light
        Vector3d sum_dir = camera_dir + light_dir;
        sum_dir.normalize();
        double specular_factor = pow(max(0.0, ray.direction.dot(sum_dir)), 1.0 * mat.shininess);
        specular_total += light_color * (attenuation_factor * specular_factor);
    }

    // Sums Ambient, Diffuse, & Specular keeping rgb values within [0, 1]
    Vector3f color = (mat.ambient.ToVector() + 
                      diffuse_total.cwiseProduct(mat.diffuse.ToVector()) + 
                      specular_total.cwiseProduct(mat.specular.ToVector())
                     ).cwiseMin(1.0f);
    return color;
}


// Part 2
void Scene::Raytrace() {
    Image img = Image(XRES, YRES);
    double height = 2 * camera.frustum.near * tan(0.5 * deg2rad(camera.frustum.fov));
    double width = camera.frustum.aspect_ratio * height;

    double pixel_height = height / YRES;
    double pixel_width = width / XRES;

    // Sets close_enough_bound based on the shortest length of a pixel
    close_enough_bound = (1.0 / 20.0) * (pixel_height < pixel_width ? pixel_height : pixel_width);

    // Gets basis vectors applying camera translations 
    // NOTEEE: (may have to use transformations not translations - so including rotation)
    Matrix4d camera_transform = camera.rotate.GetMatrix().inverse();
    Vector3d basis_e1 = (camera_transform * Vector4d(0, 0, -1, 1)).head<3>().normalized();
    Vector3d basis_e2 = (camera_transform * Vector4d(1, 0, 0, 1)).head<3>().normalized();
    Vector3d basis_e3 = (camera_transform * Vector4d(0, 1, 0, 1)).head<3>().normalized();

    for (size_t i = 0; i < XRES; i++) {
        for (size_t j = 0; j < YRES; j++) {
            Vector3f pixel_color = Vector3f::Zero();

            // Computes x and y positions of the current pixel
            double x = pixel_width * i - 0.5 * width;
            double y = pixel_height * j - 0.5 * height;

            // Computes the ray to send out at the current pixel
            Ray ray = Ray();
            // Note: This is the position of the camera in World Space
            ray.origin = -1.0 * camera.translate.GetDelta();
            ray.direction = camera.frustum.near * basis_e1 + x * basis_e2 + y * basis_e3;
    
            // Finds the closest intersection of our ray
            pair<double, Intersection> closest = ClosestIntersection(ray);
            
            if (closest.first != INFINITY) {
                // Uses helper method to compute the lighting at this pixel
                pixel_color = pointLighting(closest.second.location, 
                                            ray.origin, 
                                            lights, 
                                            closest.second.obj,
                                            root_objects);
            }

            img.SetPixel(i, j, pixel_color);
        }
    }

    // Outputs the image.
    if (!img.SaveImage("rt.png")) {
        cerr << "Error: couldn't save PNG image" << std::endl;
    } else {
        cout << "Done!\n";
    }
}
