#ifndef VKHR_MARSCHENER_GLSL
#define VKHR_MARSCHENER_GLSL

#include "../scene_graph/params.glsl"

float normal_pdf(float x, float m, float s)
{
    const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m);

    return (1 / abs(s)) * inv_sqrt_2pi * exp(- (a * a) / (2 * s * s));
}


vec3 T(vec3 absorb, float gamma_t)
{
    float l = 1 + cos(2 * gamma_t);
    return exp(-2 * absorb * l);
}


vec3 absorption(vec3 absorb, int p, float h, float refrac, float eta_parallel, float eta_perpendic)
{
    float gamma_i = asin(h);

    if (p == 0)
        return vec3(fresnel(eta_parallel, eta_perpendic, gamma_i));

    float gamma_t = asin(h / eta_perpendic);
    float fres = fresnel(eta_parallel, eta_perpendic, gamma_i);
    float inv_fres = fresnel(1 / eta_parallel, 1 / eta_perpendic, gamma_t);

    vec3 t = T(absorb, gamma_t);
    vec3 t_pow = vec3(pow(t.r, p), pow(t.g, p), pow(t.b, p));
    return (1 - fres) * (1 - fres) * pow(inv_fres, p - 1) * t_pow;
}


vec4 calc_roots(int p, float eta_perpendic, float phi)
{
    float c = asin(1 / eta_perpendic);
    return cubic_solver(-8 * (p * c / (M_PI * M_PI * M_PI)), 0, (6 * p * c / M_PI) - 2, p * M_PI - phi);
}


float B_index(float theta, float eta)
{
    float sin_theta = sin(theta);
    return sqrt(eta * eta - sin_theta * sin_theta) / cos(theta);
}


float InverseFirstDerivate(int p, float eta_perpendic, float h)
{
    float gamma_i = asin(h);
    float c = asin( 1 / eta_perpendic );
    float d_gamma = (6 * p * c / M_PI - 2) - 
		3 * 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i * gamma_i;

    return sqrt(1 - h * h) / d_gamma;
}

vec3 NP(int p, float phi, float theta_d, float eta, vec3 absorb)
{
    float eta_perpendic = B_index(theta_d, eta);
    float eta_parallel = (eta * eta) / eta_perpendic;

    vec4 roots = calc_roots(p , eta_perpendic, phi);
    vec3 res = vec3(0);

    for (int i = 0; i < roots[3]; i++ )
    {
        float gamma_i = roots[i];
        float h = sin(gamma_i);

        vec3 finalAbsorption = absorption(absorb, p, h, eta, eta_perpendic, eta_parallel);
        float inverseDerivateAngle = InverseFirstDerivate(p, eta_perpendic, h);
        res += finalAbsorption * 2 * abs(inverseDerivateAngle); //0.5 here
    }
    vec3 final_res = vec3(min(res.r, 1), min(res.g, 1), min(res.b, 1));
    return final_res;
}

float distribution(int p, float eta, float theta_d, float h)
{
    float cos_theta = cos(theta_d);
    float sin_theta = sin(theta_d);

    float eta_tick = pow(abs(eta * eta) - abs(sin_theta * sin_theta),0.5) / cos_theta;
    float c = asin(1/eta_tick);

    float gamma_i = asin(h);
    float d_phi = ((6 * p * c) / M_PI - 2) - 3 * (8 * p * c) / pow(M_PI, 3) * pow(gamma_i, 2);
    float d_h = sqrt(1 - pow(h, 2));
    return d_phi / d_h;
}

vec3 marschener(vec3 tangent, vec3 viewDirection, vec3 lightDirection, vec3 light_bulb_color, vec3 abs_coef){
    float long_width_R = longitudinal_width;     //Beta R
    float long_width_TT = long_width_R / 2;      //Beta TT
    float long_width_TRT = long_width_R * 2;     //Beta TRT

    float long_shift_R = longitudinal_shift;     //Alpha R
    float long_shift_TT = -long_shift_R / 2;     //Alpha TT
    float long_shift_TRT = -3 * long_shift_R / 2;//Alpha TRT

    //the projection of the light/view direction onto the normal of the normal plane (the tangent)
    float dotLightTangent = dot(lightDirection, tangent);
    float dotViewTangent = dot(viewDirection, tangent);
    vec3 LT_P = normalize(lightDirection - dotLightTangent);
    vec3 VT_P = normalize(viewDirection - dotViewTangent);

    vec3 Normal = normalize(LT_P + VT_P);
    vec3 Binormal = normalize(cross(Normal, tangent));

    float theta_r = acos(dot(VT_P, viewDirection)); //the angle of the view direction and the normal plane (0 = perpendicular to the hair, 90 = tangent) but in radians
    float theta_i = acos(dot(LT_P, lightDirection));//the angle of the light direction and the normal plane

    float theta_h = (theta_r + theta_i) / 2; //half angle
    float theta_d = (theta_r - theta_i) / 2; //difference angle

    float phi_r = acos(dot(Binormal, VT_P)); //the angle of the view direction and the binormal (0 = normal, 90 = binormal)
    float phi_i = acos(dot(Binormal, LT_P)); //the angle of the light direction and the binormal
    float phi = phi_r - phi_i; //the relative azimuth

    float M_R   = normal_pdf(theta_h, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(theta_h, long_shift_TT, long_width_TT);
    float M_TRT = normal_pdf(theta_h, long_shift_TRT, long_width_TRT);
    
    vec3 N_R = NP(0, phi, theta_d, refraction_index, abs_coef);
    vec3 N_TT = NP(1, phi, theta_d, refraction_index, abs_coef);
    vec3 N_TRT = NP(2, phi, theta_d, refraction_index, abs_coef);

    float cos_squared_theta_d = pow(cos(theta_d), 2);

    vec3 spec = M_R * N_R / cos_squared_theta_d + M_TT * N_TT / cos_squared_theta_d + M_TRT * N_TRT / cos_squared_theta_d;

    return light_bulb_color * spec;
}
#endif