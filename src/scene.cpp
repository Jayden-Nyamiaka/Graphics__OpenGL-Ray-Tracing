#include "scene.h"

#include <iostream>

#include "glinclude.h"

const float iotest_min = -4.0;
const float iotest_max = 4.0;
const float iotest_inc = 0.8;

void Scene::ReloadObjects() {
    // Clear buffers then tesselate.
    vertex_buffer.clear();
    normal_buffer.clear();

    for (auto &obj : root_objects) {
        obj->Tesselate(vertex_buffer, normal_buffer);
    }
}

void Scene::ReloadLighting() {
    for (unsigned int i = 0; i < lights.size(); i++) {
        int light = GL_LIGHT0 + i;

        lights[i].OpenGLInit(light);
        lights[i].OpenGLRender(light);
    }
}

void Scene::OpenGLSetup() {
    // Setup lights.
    glEnable(GL_LIGHTING);
    ReloadLighting();

    glGenVertexArrays(1, &buffer_array);
    glBindVertexArray(buffer_array);
    glGenBuffers(2, buffer_objects);

    ReloadObjects();
}

void Scene::OpenGLRender() {
    // Setup lights.
    for (unsigned int i = 0; i < lights.size(); i++) {
        int light = GL_LIGHT0 + i;

        lights[i].OpenGLRender(light);
    }

    glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[0]);
    glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * vertex_buffer.size(), vertex_buffer.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[1]);
    glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * normal_buffer.size(), normal_buffer.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);

    for (auto &obj : root_objects) {
        obj->OpenGLRender();
    }
}

void Scene::IOTest() {
    glPointSize(6.0);
    glEnable(GL_POINT_SMOOTH);

    float inside_color[3] = { 0.0, 0.6, 1.0 };

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, inside_color);
    glBegin(GL_POINTS);

    for (float x = iotest_min; x <= iotest_max; x += iotest_inc) {
        for (float y = iotest_min; y <= iotest_max; y += iotest_inc) {
            for (float z = iotest_min; z <= iotest_max; z += iotest_inc) {
                Eigen::Vector3d point = {x, y, z};
                for (auto &obj : root_objects) {
                    bool inside = obj->IOTest(point);
                    if (inside) {
                        glVertex3dv(point.data());
                    }
                }
            }
        }
    }

    glEnd();
}

void Scene::DrawIntersectTest() {
    Ray incoming = Ray { Eigen::Vector3d(0, 0, -5), Eigen::Vector3d(0, 0, 1) };
    auto closest = ClosestIntersection(incoming);
    Ray intersection = closest.second.location;

    float outside_color[3] = { 1.0, 1.0, 1.0 };

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, outside_color);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glVertex3dv(incoming.origin.data());
    glVertex3dv(incoming.At(1.0).data());
    glEnd();

    glBegin(GL_LINES);
    glVertex3dv(intersection.origin.data());
    glVertex3dv(intersection.At(1.0).data());
    glEnd();
}

std::pair<float, Intersection> Scene::ClosestIntersection(const Ray &incoming) const {
    std::pair<float, Intersection> closest = std::make_pair(INFINITY, Intersection());

    for (auto &obj : root_objects) {
        auto temp = obj->ClosestIntersection(incoming);

        if (temp.first < closest.first) {
            closest = temp;
        }
    }

    return closest;
}