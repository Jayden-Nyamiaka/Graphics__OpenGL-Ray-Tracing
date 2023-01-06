#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "glinclude.h"
#include "light.h"
#include "object.h"
#include "scene.h"
#include "transform.h"
#include "util.h"

#include "parsing.h"

// Renderer Usage String
const std::string usage = "Usage: renderer <xres> <yres> [scene_file.yaml]";

int xres;
int yres;

Scene scene;

GLenum shader;

bool wireframe = false;
bool io_test = false;
bool intersect_test = false;
bool scene_cam = false;

Arcball arcball;

/**
 * Handles GLUT window reshape event.
 */
void reshape(int width, int height) {
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;

    glViewport(0, 0, width, height);
    arcball.SetRes(width, height);
    glutPostRedisplay();
}

void init() {
    // Set up backface culling and depth buffering.
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Set up the projection view.
    glMatrixMode(GL_PROJECTION);
    gluPerspective(60.0, 1.0, 1.0, 80.0);

    glMatrixMode(GL_MODELVIEW);

    scene.OpenGLSetup();
}

void shader_setup() {
    std::string vert_source, frag_source;

    std::ifstream vert_file("shaders/vertex.glsl");
    if (!vert_file)
        std::cerr << "Error opening vertex shader program\n";

    std::ifstream frag_file("shaders/fragment.glsl");
    if (!frag_file)
        std::cerr << "Error opening fragment shader program\n";

    getline(vert_file, vert_source, '\0');
    const char *vert_shader_source = vert_source.c_str();

    getline(frag_file, frag_source, '\0');
    const char *frag_shader_source = frag_source.c_str();

    char buf[1024];
    GLsizei blah;

    // Initialize shaders
    GLenum vertex, fragment;

    shader = glCreateProgram();

    vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, &vert_shader_source, NULL);
    glCompileShader(vertex);

    GLint isCompiled = 0;
    glGetShaderiv(vertex, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        GLint maxLength = 0;
        glGetShaderiv(vertex, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> errorLog(maxLength);
        glGetShaderInfoLog(vertex, maxLength, &maxLength, &errorLog[0]);

        // Provide the infolog in whatever manor you deem best.
        // Exit with failure.
        for (auto &entry : errorLog)
            std::cout << entry;
        glDeleteShader(vertex); // Don't leak the shader.
        return;
    }

    fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, &frag_shader_source, NULL);
    glCompileShader(fragment);

    glGetShaderiv(fragment, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        GLint maxLength = 0;
        glGetShaderiv(fragment, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> errorLog(maxLength);
        glGetShaderInfoLog(fragment, maxLength, &maxLength, &errorLog[0]);

        // Provide the infolog in whatever manor you deem best.
        // Exit with failure.
        for (auto &entry : errorLog)
            std::cout << entry;
        glDeleteShader(fragment); // Don't leak the shader.
        return;
    }

    glAttachShader(shader, vertex);
    glAttachShader(shader, fragment);

    glLinkProgram(shader);
    std::cout << "Enabling fragment program: " << gluErrorString(glGetError())
         << "\n";
    glGetProgramInfoLog(shader, 1024, &blah, buf);
    std::cout << buf;

    std::cout << "Enabling program object\n";
    glUseProgram(shader);
}

/**
 * Handles GLUT window draw event.
 */
void display() {
    // Reset GL stuff.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glLoadIdentity();

    if (scene_cam) {
        scene.GetCamera().OpenGLSetPosition();
    } else {
        glTranslatef(0.0, 0.0, -8.0);
        
        arcball.Apply();
    }

    if (wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    scene.OpenGLRender();

    if (io_test) {
        scene.IOTest();
    }

    if (intersect_test) {
        scene.DrawIntersectTest();
    }
    
    // Swap in the new buffer.
    glutSwapBuffers();
}

/**
 * Handles GLUT mouse press event.
 */
void mouse_pressed(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        // Handles left click.
        arcball.Start(x, y);
    } else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        // Handles left click release.
        arcball.End();
    }
}

/**
 * Handles GLUT mouse move event.
 */
void mouse_moved(int x, int y) {
    if (!scene_cam) {
        arcball.Move(x, y);
        glutPostRedisplay();
    }
}

/**
 * Handles GLUT key press event.
 */
void key_pressed(unsigned char key, int x, int y) {
    switch (key) {
        case 't': {
            wireframe = !wireframe;
            if (wireframe) {
                glDisable(GL_CULL_FACE);
            } else {
                glEnable(GL_CULL_FACE);
            }
            break;
        };
        case 'o': {
            io_test = !io_test;
            break;
        };
        case 'i': {
            intersect_test = !intersect_test;
            break;
        };
        case 'c': {
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();

            scene_cam = !scene_cam;
            if (scene_cam) {
                scene.GetCamera().frustum.OpenGLSetMatrix();
            } else {
                gluPerspective(60.0, 1.0, 1.0, 80.0);
            }

            glMatrixMode(GL_MODELVIEW);
            break;
        };
        case 'r': {
            std::thread draw_thread(&Scene::Raytrace, &scene);
            draw_thread.detach();
            break;
        };
        case 'q': {
            exit(0);
        };
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {
    // Grab the xres and yres from arguments.
    if (argc >= 3) {
        try {
            xres = std::stoi(argv[1]);
            yres = std::stoi(argv[2]);
        } catch (int e) {
            std::cerr << usage << "\n";
            exit(e);
        }
    } else {
        std::cerr << usage << "\n";
        exit(1);
    }

    std::cout << yres << "\n";

    arcball = Arcball(xres, yres);

    // Initialize GLUT.
    glutInit(&argc, argv);

    // Logic for centering window on screen.
    int screen_width  = glutGet(GLUT_SCREEN_WIDTH);
    int screen_height = glutGet(GLUT_SCREEN_HEIGHT);
    int pos_x = (screen_width - xres) / 2;
    int pos_y = (screen_height - yres) / 2;

    // Set up GLUT window.
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(xres, yres);
    glutInitWindowPosition(pos_x, pos_y);
    glutCreateWindow("CS 171 Renderer");

    YAML::Node ln;
    if (argc > 3) {
        ln = YAML::LoadFile(argv[3]);
    } else {
        ln = YAML::LoadFile("scenes/sphere.yaml");
    }
    scene = ln.as<Scene>();

    shader_setup();
    init();

    // Setup GLUT functions.
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse_pressed);
    glutMotionFunc(mouse_moved);
    glutKeyboardFunc(key_pressed);

    // Start drawing.
    glutMainLoop();
}