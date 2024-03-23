#ifndef VKHR_KAJIYA_KAY_DIFF_GLSL
#define VKHR_KAJIYA_KAY_DIFF_GLSL

#define KAJIYA_KAY 0

// Based on "Rendering Hair with Three Dimensional Textures" by J. T. Kajiya and T. L. Kay.
vec3 kajiya_kay_diffuse(vec3 diffuse, vec3 tangent, vec3 light) {
    float cosTL = dot(tangent, light);
    float cosTL_squared = cosTL*cosTL;
    float one_minus_cosTL_squared = 1.0f - cosTL_squared;
    float sinTL = sqrt(one_minus_cosTL_squared);

    return diffuse * sinTL;
}

#endif