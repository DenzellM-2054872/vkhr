#ifndef VKHR_FRESNEL_GLSL
#define VKHR_FRESNEL_GLSL

float fresnelPar(float eta_2, float theta)
{
    float eta_1 = 1;

    float cos_gamma_i = cos(theta);
    float a = ((eta_1 / eta_2) * sin(theta));
    float b = a * a;

    if (b > 1) return 1;

    float cos_gamma_t = sqrt(1 - b);

    float R = (eta_2 * cos_gamma_i - eta_1 * cos_gamma_t) / 
			  (eta_2 * cos_gamma_i + eta_1 * cos_gamma_t);

    return min(1, R * R);
}


float fresnelPer(float eta_2, float theta)
{
    float eta_1 = 1;

    float cos_gamma_i = cos(theta);
    float a = ((eta_1 / eta_2) * sin(theta));
    float b = a * a;

    if (b > 1) return 1;

    float cos_gamma_t = sqrt(1 - b);

    float R = (eta_1 * cos_gamma_i - eta_2 * cos_gamma_t) / 
			  (eta_1 * cos_gamma_i + eta_2 * cos_gamma_t);

    return min(1, R * R);
}


float fresnel(float eta_parallel, float eta_perpendic, float theta)
{
    return 0.5f * (fresnelPer(eta_perpendic, theta) + 
				    fresnelPar(eta_parallel, theta));
}
#endif