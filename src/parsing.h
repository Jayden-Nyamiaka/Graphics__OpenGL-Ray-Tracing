#ifndef PARSING_H
#define PARSING_H

#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>

#include "light.h"
#include "object.h"
#include "scene.h"
#include "util.h"

namespace YAML {

template <>
struct convert<Eigen::Vector3d> {
    static Node encode(const Eigen::Vector3d &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);

        return node;
    }

    static bool decode(const Node &node, Eigen::Vector3d &rhs) {
        if (!node.IsSequence() || node.size() != 3) {
            return false;
        }

        rhs[0] = node[0].as<double>();
        rhs[1] = node[1].as<double>();
        rhs[2] = node[2].as<double>();

        return true;
    }
};

template <>
struct convert<Eigen::Vector4d> {
    static Node encode(const Eigen::Vector4d &rhs) {
        Node node;
        node.push_back(rhs[0]);
        node.push_back(rhs[1]);
        node.push_back(rhs[2]);
        node.push_back(rhs[3]);

        return node;
    }

    static bool decode(const Node &node, Eigen::Vector4d &rhs) {
        if (!node.IsSequence() || node.size() != 4) {
            return false;
        }

        rhs[0] = node[0].as<double>();
        rhs[1] = node[1].as<double>();
        rhs[2] = node[2].as<double>();
        rhs[3] = node[3].as<double>();

        return true;
    }
};

template <>
struct convert<Rotate> {
    static Node encode(const Rotate &rhs) {
        Node node;
        node["axis"] = rhs.GetAxis();
        node["angle"] = rhs.GetAngle();

        return node;
    }

    static bool decode(const Node &node, Rotate &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.SetAxisAngle(node["axis"].as<Eigen::Vector3d>(), node["angle"].as<float>());

        return true;
    }
};

template <>
struct convert<Translate> {
    static Node encode(const Translate &rhs) {
        Node node;
        node["delta"] = rhs.GetDelta();

        return node;
    }

    static bool decode(const Node &node, Translate &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.SetDelta(node["delta"].as<Eigen::Vector3d>());

        return true;
    }
};

template <>
struct convert<Scale> {
    static Node encode(const Scale &rhs) {
        Node node;
        node["scale"] = rhs.GetScale();

        return node;
    }

    static bool decode(const Node &node, Scale &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.SetScale(node["scale"].as<Eigen::Vector3d>());

        return true;
    }
};

template <>
struct convert<Color> {
    static Node encode(const Color &rhs) {
        Node node;
        node.push_back(rhs.r);
        node.push_back(rhs.g);
        node.push_back(rhs.b);

        return node;
    }

    static bool decode(const Node &node, Color &rhs) {
        if (!node.IsSequence() || node.size() != 3) {
            return false;
        }

        rhs.r = node[0].as<float>();
        rhs.g = node[1].as<float>();
        rhs.b = node[2].as<float>();

        return true;
    }
};

template <>
struct convert<Material> {
    static Node encode(const Material &rhs) {
        Node node;
        node["ambient"] = rhs.ambient;
        node["diffuse"] = rhs.diffuse;
        node["specular"] = rhs.specular;
        node["shininess"] = rhs.shininess;

        return node;
    } 

    static bool decode(const Node &node, Material &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.ambient = node["ambient"].as<Color>();
        rhs.diffuse = node["diffuse"].as<Color>();
        rhs.specular = node["specular"].as<Color>();
        rhs.shininess = node["shininess"].as<float>();

        return true;
    }
};

template <>
struct convert<Light> {
    static Node encode(const Light &rhs) {
        Node node;
        node["position"] = rhs.position;
        node["color"] = rhs.color;
        node["attenuation"] = rhs.attenuation;

        return node;
    }

    static bool decode(const Node &node, Light &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.position = node["position"].as<Eigen::Vector4d>();
        rhs.color = node["color"].as<Color>();
        rhs.attenuation = node["attenuation"].as<float>();

        return true;
    }
};

template <>
struct convert<Frustum> {
    static Node encode(const Frustum &rhs) {
        Node node;
        node["aspect_ratio"] = rhs.aspect_ratio;
        node["fov"] = rhs.fov;
        node["near"] = rhs.near;
        node["far"] = rhs.far;

        return node;
    }

    static bool decode(const Node &node, Frustum &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.aspect_ratio = node["aspect_ratio"].as<double>();
        rhs.fov = node["fov"].as<double>();
        rhs.near = node["near"].as<double>();
        rhs.far = node["far"].as<double>();

        return true;
    }
};

template <>
struct convert<Camera> {
    static Node encode(const Camera &rhs) {
        Node node;
        node["translate"] = rhs.translate;
        node["rotate"] = rhs.rotate;
        node["frustum"] = rhs.frustum;

        return node;
    }

    static bool decode(const Node &node, Camera &rhs) {
        if (!node.IsMap()) {
            return false;
        }

        rhs.translate = node["translate"].as<Translate>();
        rhs.rotate = node["rotate"].as<Rotate>();
        rhs.frustum = node["frustum"].as<Frustum>();

        return true;
    }
};

template<>
struct convert<Scene> {
    static Node encode(const Scene &scene) {
        return Node();
    }

    static bool decode(const Node &node, Scene &scene) {
        if (!node.IsMap()) {
            return false;
        }

        // Parse Camera
        if (!node["camera"].IsMap()) {
            return false;
        }
        
        std::cout << "Parsing camera...\n";
        scene.SetCamera(node["camera"].as<Camera>());

        // Parse Lights
        if (!node["lights"].IsSequence()) {
            return false;
        }

        std::cout << "Parsing lights...\n";
        for (auto &child : node["lights"]) {
            scene.AddLight(child.as<Light>());
        }

        // Parse Objects
        if (!node["objects"].IsSequence()) {
            return false;
        }

        std::cout << "Parsing objects...\n";

        std::map<std::string, std::vector<std::string>> assemblies;

        // Create pointers for all of the objects.
        for (auto &obj_data : node["objects"]) {
            if (!obj_data.IsMap()) {
                return false;
            }

            std::string name = obj_data["name"].as<std::string>();
            std::string type = obj_data["type"].as<std::string>();
            std::shared_ptr<Object> object;

            std::cout << "Parsing " << type << " " << name << "...\n";

            if (type == "superquadric") {
                float exp0 = obj_data["exp0"].as<float>();
                float exp1 = obj_data["exp1"].as<float>();
                auto obj = std::make_shared<Superquadric>(exp0, exp1);

                Material mat = obj_data["material"].as<Material>();
                obj->SetMaterial(mat);

                object = obj;
            } else if (type == "assembly") {
                auto obj = std::make_shared<Assembly>();

                assemblies.emplace(name, obj_data["children"].as<std::vector<std::string>>());

                object = obj;
            }

            // Parse transforms.
            if (!obj_data["transforms"].IsSequence()) {
                return false;
            }

            for (auto &transform : obj_data["transforms"]) {
                std::string ttype = transform["type"].as<std::string>();
                if (ttype == "rotate") {
                    object->AddTransform(transform.as<Rotate>());
                } else if (ttype == "scale") {
                    object->AddTransform(transform.as<Scale>());
                } else if (ttype == "translate") {
                    object->AddTransform(transform.as<Translate>());
                }
            }

            bool root = obj_data["root"].as<bool>();
            scene.AddObject(name, object, root);
        }

        // Go back and build the assemblies tree.
        for (auto &entry : assemblies) {
            auto assembly = std::static_pointer_cast<Assembly>(scene.GetObject(entry.first));
            for (auto child : entry.second) {
                assembly->AddChild(scene.GetObject(child));
            }
        }

        return true;
    }
};

}  // namespace YAML

#endif  // PARSING_H