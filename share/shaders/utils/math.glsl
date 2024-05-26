#ifndef VKHR_MATH_GLSL
#define VKHR_MATH_GLSL

#define M_PI 3.14159265358979323846
#define M_E  2.71828182845904523536
#define M_EPS 0.00001
#define M_INT_PRES 

vec3 vec_pow(vec3 vec, float expon){
	return vec3(pow(vec.r, expon), pow(vec.g, expon), pow(vec.b, expon));
}

float gauss(float x_m, float s)
{
    const float inv_sqrt_2pi = 0.3989422804014327;
    return (1 / abs(s)) * inv_sqrt_2pi * exp(- (x_m * x_m) / (2 * s * s));
}
#endif
