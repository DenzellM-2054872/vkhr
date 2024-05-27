#version 460 core
#include "../scene_graph/camera.glsl"
#include "../self-shadowing/approximate_deep_shadows.glsl"
#include "../volumes/local_ambient_occlusion.glsl"

#include "../transparency/ppll.glsl"
#include "../level_of_detail/scheme.glsl"
#include "../anti-aliasing/gpaa.glsl"

#include "../scene_graph/lights.glsl"
#include "../scene_graph/shadow_maps.glsl"
#include "../scene_graph/params.glsl"

#include "strand.glsl"

#include "../shading/kajiya-kay-diffuse.glsl"
#include "../shading/marschener.glsl"

layout(early_fragment_tests) in;

layout(location = 0) in PipelineIn {
    vec4 position;
    vec3 tangent;
    float thickness;
} fs_in;

layout(push_constant) uniform Object {
    mat4 model;
} object;



layout(location = 0) out vec4 color;

void main() {
    ivec2 pixel = ivec2(gl_FragCoord.xy);
    uint node = ppll_next_node();
    if (node == PPLL_NULL_NODE) discard;

    float coverage = gpaa(gl_FragCoord.xy, fs_in.position,
                          camera.projection * camera.view,
                          camera.resolution, strand_width);

    coverage *= hair_alpha; // Alpha used for transparency.
    if (coverage < 0.001) discard; // Shading not worth it!

    coverage *= fs_in.thickness * STRAND_SCALING; // Slowly fades the strand at the tip.

    vec3 tangent = normalize(fs_in.tangent);
    //in the paper both of these vectors are defined as going from hair to point
    vec3 viewDirection = normalize(camera.position - fs_in.position.xyz);
    vec3 lightDirection = normalize(lights[0].origin  - fs_in.position.xyz);

    vec3 diffuse_colors = kajiya_kay_diffuse(hair_color, tangent, lightDirection);

    vec3 abs_coef = vec3(abs_coef_R, abs_coef_G, abs_coef_B);
    vec3 light_bulb_color = lights[0].intensity;
    
    vec3 specular_colors = marschener(tangent, viewDirection, lightDirection, light_bulb_color, abs_coef, hair_color);
    vec3 shading = specular_colors + diffuse_colors;

    vec4 shadow_space_fragment = lights[0].matrix * fs_in.position;


    float occlusion = approximate_deep_shadows(shadow_maps[0],
                                            shadow_space_fragment,
                                            deep_shadows_kernel_size,
                                            deep_shadows_stride_size,
                                            15000.0f, hair_alpha);

    color = vec4(shading * occlusion, coverage);


    ppll_node_data(node, color, gl_FragCoord.z);
    ppll_link_node(pixel, node);

    discard; // Fragments resolved in next pass.

}

