varying vec4 vertex;
varying vec3 normal;

struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    float shininess;
};

uniform Material material;

void main() {
    // Color Sums
    vec4 ambient = vec4(0.0);
    vec4 diffuse = vec4(0.0);
    vec4 specular = vec4(0.0);

    vec3 l, e, r;

    for (int i = 0; i < gl_MaxLights; i++) {
        // Get the l, e, and r vectors.
        l = normalize(vec3(gl_LightSource[i].position - vertex));
        e = normalize(vec3(-vertex));
        r = normalize(reflect(-l, normal));

        ambient += gl_FrontLightProduct[i].ambient;
        diffuse += gl_FrontLightProduct[i].diffuse * max(dot(normal, l), 0.0);
        specular += gl_FrontLightProduct[i].specular *
                    pow(max(dot(r, e), 0.0), gl_FrontMaterial.shininess);
    }

    gl_FragColor = clamp(ambient + specular + diffuse, 0.0, 1.0);
}