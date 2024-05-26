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

float fresnelPar2(float eta_2, float theta)
{
    float eta_1 = 1;

    float cos_gamma_i = cos(theta);
    float cos_gamma_t = asin(sin(theta) / eta_2);

    float R =   (eta_2 * cos_gamma_i - eta_1 * cos_gamma_t) /
                (eta_2 * cos_gamma_i + eta_1 * cos_gamma_t);

    return min(1, R * R);
}

float fresnelPer2(float eta_2, float theta)
{
    float eta_1 = 1;

    float cos_gamma_i = cos(theta);
    float cos_gamma_t = asin(sin(theta) / eta_2);

    float R =   (eta_1 * cos_gamma_i - eta_2 * cos_gamma_t) /
                (eta_1 * cos_gamma_i + eta_2 * cos_gamma_t);

    return min(1, R * R);
}

float F_02(float eta_parallel, float eta_perpendic, float theta)
{

    return 0.5f * (fresnelPer2(eta_perpendic, theta) +
        fresnelPar2(eta_parallel, theta));
}
float F_0(float eta_parallel, float eta_perpendic, float theta)
{

    return 0.5f * (fresnelPer(eta_perpendic, theta) + 
				   fresnelPar(eta_parallel, theta));
}

float fresnel2(float eta, float x)
{
    float f_0 = pow(1 - eta, 2) / pow(1 + eta, 2);
    return f_0 + (1 - f_0) * pow(1 - x, 5);
}


float fresnel(float eta_parallel, float eta_perpendic, float theta)
{
    float cos_theta = cos(theta);
    
    float flipped = 1.0f - cos_theta;
    float flipped_squared = flipped * flipped;

    float fresnel_0 = fresnelPar(eta_parallel, theta);
    float fresnel_90 = fresnelPer(eta_perpendic, theta);

    return fresnel_0 + (fresnel_90 - fresnel_0) * (flipped_squared * flipped * flipped_squared);
}

#endif