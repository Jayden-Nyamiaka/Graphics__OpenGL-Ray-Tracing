attribute vec3 vertex_in;
attribute vec3 normal_in;

varying vec4 vertex;
varying vec3 normal;

void main() {
    // Transform vertex and normal to camera space.
    vertex = gl_ModelViewMatrix * vec4(vertex_in, 1.0);
    normal = normalize(gl_NormalMatrix * normal_in);

    // Set NDC position.
    gl_Position = gl_ProjectionMatrix * vertex;
}