#ifndef SCENE_H
#define SCENE_H

#include <fstream>
#include <map>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "camera.h"
#include "light.h"
#include "object.h"

class Scene {
private:
    std::vector<Light> lights;
    std::vector<Eigen::Vector3f> vertex_buffer;
    std::vector<Eigen::Vector3f> normal_buffer;
    std::vector<std::shared_ptr<Object>> root_objects;
    std::map<std::string, std::shared_ptr<Object>> objects;

    Camera camera;

    unsigned int buffer_array;
    unsigned int buffer_objects[2];
public:
    Scene() = default;
    Scene(std::ifstream &scene_file);

    void ReloadObjects();
    void ReloadLighting();
    void OpenGLSetup();
    void OpenGLRender();
    void IOTest();
    void DrawIntersectTest();
    void Raytrace();

    std::pair<float, Intersection> ClosestIntersection(const Ray &incoming) const;

    void SetCamera(const Camera &cam) {
        camera = cam;
    }

    const Camera& GetCamera() const {
        return camera;
    }

    const std::vector<Light> GetLights() const {
        return lights;
    }

    void AddLight(const Light &light) {
        lights.push_back(light);
    }

    void AddObject(std::string name, std::shared_ptr<Object> obj, bool is_root = false) {
        objects.emplace(name, obj);
        if (is_root) {
            root_objects.push_back(obj);
        }
    }

    std::shared_ptr<Object> GetObject(std::string name) {
        if (objects.find(name) == objects.end()) {
            return nullptr;
        }
        
        return objects[name];
    }
};

#endif // SCENE_H