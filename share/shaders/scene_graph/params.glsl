#ifndef VKHR_PARAMS_GLSL
#define VKHR_PARAMS_GLSL

#define YES  1
#define NO   0

layout(binding = 4) uniform Params {
    int shading_model;

    int deep_shadows_kernel_size;
    int deep_shadows_sampling_type;
    int deep_shadows_stride_size;
    int deep_shadows_on;

    int pcf_shadows_kernel_size;
    int pcf_shadows_sampling_type;
    float pcf_shadows_bias;
    int pcf_shadows_on;

    int shadow_technique;

    float isosurface;
    float raycast_steps;
    float occlusion_radius;
    float ao_exponent;
    float ao_max;

    float magnified_distance;
    int renderer;
    float minified_distance;

    int benchmarking;


    float longitudinal_width;       // BetaR which sould be between 5 & 10 degrees
    float longitudinal_shift;       // AlpaR which sould be between -10 & -5 degrees
    float refraction_index;         // Nu should be around 1.55

    float abs_coef_R;               //(r, g, b) 0.2 - inf
    float abs_coef_G;               //(r, g, b) 0.2 - inf
    float abs_coef_B;               //(r, g, b) 0.2 - inf

    float caustic_intensity_limit;  // 0.5
    float caustic_width;            // 10 to 25 degrees
    float glint_scale_fac;          // 0.5 to 5
    float caustic_merge_range;      // 0.2 to 0.4
};

#endif
