#version 460 core
#include "strand.glsl"

#include "../scene_graph/params.glsl"
#include "../scene_graph/camera.glsl"
#include "../scene_graph/lights.glsl"

#include "../utils/cubic_solver.glsl"
#include "../utils/fresnel.glsl"
#include "../utils/math.glsl"

#include "../shading/kajiya-kay-diffuse.glsl"
#include "../shading/marschener.glsl"


layout(location = 0) in PipelineIn {
    vec4 position;
    vec3 tangent;
    float thickness;
} fs_in;

layout(location = 0) out vec4 color;

void main() {
    vec3 tangent = normalize(fs_in.tangent);
    vec3 viewDirection = normalize(fs_in.position.xyz - camera.position);
    vec3 lightDirection = normalize(lights[0].origin - fs_in.position.xyz);

    vec3 diffuse_colors = kajiya_kay_diffuse(hair_color, tangent, lightDirection);

    vec3 abs_coef = vec3(abs_coef_R, abs_coef_G, abs_coef_B);
    vec3 light_bulb_color = lights[0].intensity;
    vec3 specular_colors = marschener(tangent, viewDirection, lightDirection, light_bulb_color, abs_coef);

    color = vec4(specular_colors + diffuse_colors, 1.0f);
}

