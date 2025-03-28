in vec3 attr_pos;
uniform mat4 modelViewProjection;

void main() {
    gl_Position = modelViewProjection * vec4(attr_pos, 1.0);
}
