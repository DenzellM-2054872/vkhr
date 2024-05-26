#ifndef VKHR_MARSCHENER_GLSL
#define VKHR_MARSCHENER_GLSL

#include "../scene_graph/params.glsl"
#include "../utils/fresnel.glsl"
#include "../utils/cubic_solver.glsl"

float normal_pdf(float x, float m, float s)
{
    const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m);

    return (1 / abs(s)) * inv_sqrt_2pi * exp(- (a * a) / (2 * s * s));
}

float normal_pdf(float x_m, float s)
{
    const float inv_sqrt_2pi = 0.3989422804014327;

    return (1 / abs(s)) * inv_sqrt_2pi * exp(-(x_m * x_m) / (2 * s * s));
}

vec3 T(vec3 absorb, float gamma_t)
{
    float l = 1 + cos(2 * gamma_t);
    return exp(-2 * absorb * l);
}


vec3 attenuation(vec3 absorb, int p, float h, float eta_parallel, float eta_perpendic, float theta_d)
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

    //karis says to do it like this (frostbite has no specification)
    //float gamma_t = asin(h / eta_perpendic);
    //float fres = fresnel(eta_parallel, eta_perpendic, cos(theta_d) * sqrt(1-(h * h)));

    //vec3 t = T(absorb, gamma_t);
    //vec3 t_pow = vec3(pow(t.r, p), pow(t.g, p), pow(t.b, p));

    //return (1 - fres) * (1 - fres) * pow(fres, p - 1) * t_pow;
}


vec4 calc_roots(int p, float eta_perpendic, float phi)
{
    float c = asin(1 / eta_perpendic);
    return cubic_solver(-8 * (p * c / (M_PI * M_PI * M_PI)), 0, (6 * p * c / M_PI) - 2, p * M_PI - phi);
}


float B_index(float theta, float eta)
{
    //if(eta == 1.55) return 1.19 / cos(theta) + 0.36 * cos(theta);

    float sin_theta = sin(theta);
    return sqrt(eta * eta - sin_theta * sin_theta) / cos(theta);
}

float C_index(float theta, float eta)
{
    float sin_theta = sin(theta);
    //return eta  * cos(theta) / (sqrt(1 - pow(eta, -2) * pow(sin(theta), 2)));
    return eta * eta * cos(theta) / (sqrt(eta * eta - sin_theta * sin_theta));
}


float inv_first_der(int p, float eta_perpendic, float h)
{
    float gamma_i = asin(h);
    float c = asin( 1 / eta_perpendic );
    float d_gamma = ((6 * p * c / M_PI) - 2) - 
		3 * 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i * gamma_i;

    return sqrt(1 - h * h) / d_gamma;
    //return pow((6 * p * c / M_PI * sqrt(1 - h * h)) - (2 / sqrt(1 - h * h)) - (3 * 8 * (p * c / (M_PI * M_PI * M_PI * sqrt(1 - pow(h, 6))))), -1);
    //return 1 / d_gamma;
}

float first_der(int p, float eta_perpendic, float h)
{
    float gamma_i = asin(h);
    float c = asin(1 / eta_perpendic);
    float d_gamma = ((6 * p * c / M_PI) - 2) -
        3 * 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i * gamma_i;

    return d_gamma / sqrt(1 - h * h);
}

float inv_second_der(int p, float eta_perpendic, float h)
{
    float gamma_i = asin(h);
    float c = asin(1 / eta_perpendic);
            
    float f = ((6 * p * c / M_PI) - 2) - 3 * 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i * gamma_i;
    float df = -2 * 3 * 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i;

    float g = sqrt(1 - h * h);
    float dg = h / max( g , M_EPS);

    return (g * g) / (max( df * g - f * dg , M_EPS ));
}

float angle_polynomial(int p, float eta_perpendic, float h)
{
    float gamma_i = asin(h);
    float c = asin(1 / eta_perpendic);

    return ((6 * p * c / M_PI) - 2) * gamma_i - 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i * gamma_i * gamma_i + p * M_PI;
}

vec3 NP(int p, float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 absorb)
{
    vec4 roots = calc_roots(p , eta_perpendic, phi);
    vec3 res = vec3(0);
    if (p == 2 && roots[3] == 4) return vec3(1, 0, 0);
    for (int i = 0; i < roots[3]; i++ )
    {
        float gamma_i = roots[i];
        float h = sin(gamma_i);
//        if(h < M_EPS) continue;

        vec3 finalAbsorption = attenuation(absorb, p, h, eta_perpendic, eta_parallel, theta_d);
        float inverseDerivateAngle = 1 / abs((2 * first_der(p, eta_perpendic, h)));
        res += finalAbsorption * inverseDerivateAngle; 
    }
    return vec3(min(res.r, 1), min(res.g, 1), min(res.b, 1));
}

vec3 NP_R(float phi, float theta_d, float eta_perpendic, float eta_parallel)
{
    float gamma_i = -phi / 2;
    float h = sin(gamma_i);

    //this one is from the implementation i found
    vec3 res = vec3(sqrt(1 - h * h));
    res *= vec3(fresnel(eta_perpendic, eta_parallel, gamma_i));


    return vec3(min(1, res.r),min(1, res.g), min(1, res.b));
}

vec3 NP_R_K(float phi, float theta_d, float eta_perpendic, float eta_parallel, float view_light_angle)
{
    vec3 res = vec3(0.25 * cos(phi / 2));
    res *= vec3(fresnel2(refraction_index, sqrt(0.5 * (1 + view_light_angle))));
//    res *= 3.1;   //this makes the graphs closer; expect more of these magical scalers going forward
    return vec3(min(1, res.r),min(1, res.g), min(1, res.b));
}

float calc_h_tt(float phi, float a){
    float top = sign(phi) * cos(phi / 2);
    float bottom = sqrt(1 + a * a - 2 * a * sign(phi) * sin(phi / 2));
    return abs(top/bottom);
}

float calc_h_2_tt(float phi, float a){
    float top = 0.5 + 0.5 * cos(phi);
    float bottom = 1 + a * a - 2 * a * sqrt(0.5 - 0.5 * cos(phi));
    return top/bottom;
}
vec3 NP_TT_K(float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 hair_color)
{

    float a = 1 / eta_perpendic;
    //float h = calc_h_tt(phi, a);
    float h = 0;

    float exponent = sqrt(1 - h * (a * a)) / (2 * cos(theta_d));
    vec3 t = vec3(pow(hair_color.r, exponent), pow(hair_color.g, exponent), pow(hair_color.b, exponent));

    float fres_angle = -0.65 * cos(theta_d) * sqrt(1 - h);
    float fres = fresnel2(refraction_index, fres_angle);

    //float fres = fresnel(eta_parallel, eta_perpendic, fres_angle);
    //given that p = 1: p-1 = 0 -> pow(x, 0) = 1 and thus can be ignored
    vec3 att = (1 - fres) * (1 - fres) * t;
    
    //from the karis paper (quite difrent then the original)
    float distrib = clamp(exp(-3.65 * cos(phi) - 3.98),0 , 0.01);
    vec3 res = att * 0.3 * distrib; 
    //vec3 res = att  * 1.2 * distrib; //another magic scaler :)

    vec3 final_res = vec3(min(res.r, 1), min(res.g, 1), min(res.b, 1));
    return final_res;
}

vec3 NP_TRT(float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 absorb)
{
    float delta_hM = caustic_intensity_limit;
    float w_c = caustic_width;
    float k_G = glint_scale_fac;
    float delta_eta_tick = caustic_merge_range;

    float delta_h, t, phi_c, h_c;

    if (eta_perpendic < 2)
    {
        float c = asin(1 / eta_perpendic);
        h_c = sqrt((4 - eta_perpendic * eta_perpendic) / 3);
        phi_c = (-8 * (2 * c / (M_PI * M_PI * M_PI))) * pow(asin(h_c), 3) + ((6 * 2 * c / M_PI) - 2) * asin(h_c) + 2 * M_PI;
        //phi_c = sqrt((6 * 2 * c / M_PI - 2) / (3 * 8 * (2 * c / (M_PI * M_PI * M_PI))));
        //h_c = abs(sin(phi_c));

        float inv_der_angle = inv_second_der(2, eta_perpendic, h_c);
        delta_h = min(delta_hM, 2 * sqrt(2 * w_c * inv_der_angle));
        t = 1;
    } else {
        phi_c = 0;
        h_c = 0;
        delta_h = delta_hM;
        t = 1 - smoothstep(2, 2 + delta_eta_tick, eta_perpendic);
    }

    //phi_c = angle_polynomial(2, eta_perpendic, h_c);
    vec3 res = NP(2, phi, theta_d, eta_perpendic, eta_parallel, absorb);
    vec3 final_abs = attenuation(absorb, 2, h_c, eta_perpendic, eta_parallel, theta_d);

    res = res * (1 - t * normal_pdf(phi - phi_c, w_c) / normal_pdf(0, w_c));
    res = res * (1 - t * normal_pdf(phi + phi_c, w_c) / normal_pdf(0, w_c));
    res = res + t * k_G * final_abs * delta_h * (normal_pdf(phi - phi_c, w_c) + normal_pdf(phi + phi_c, w_c));

    return vec3(min(1, res.r),min(1, res.g), min(1, res.b));
}

vec3 NP_TRT_K(float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 hair_color, float beta_n)
{
    // should be clamped but idk to what values
    float scale = 1;
    float distrib = scale * exp(scale * (17 * cos(phi) - 16.78));

    float h = sqrt(3) / 2;

    float exponent = 0.8 / cos(theta_d);
    vec3 t = vec3(pow(hair_color.r, exponent), pow(hair_color.g, exponent), pow(hair_color.b, exponent));

    vec3 t_pow = vec3(pow(t.r, 2), pow(t.g, 2), pow(t.b, 2));

    float fres_angle = cos(theta_d) * sqrt(1 - pow(h, 2));
    float fres = fresnel2(refraction_index, fres_angle);

    //given that p = 2: p-1 = 1 -> pow(x, 1) = x and thus can be ignored
    vec3 att = (1 - fres) * (1 - fres) * fres * t_pow;

    vec3 res = att * 1 * distrib; // another magic scaler :>

    return  vec3(min(1, res.r),min(1, res.g), min(1, res.b));
}

vec3 get_R(float phi, float theta_d, float eta_perpendic, float eta_parallel, float view_light_angle, vec3 abs_coef, bool using_karis) {
    if (using_karis) return NP_R_K( phi, theta_d, eta_perpendic, eta_parallel, view_light_angle);
    return NP(0, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
}

vec3 get_TT(float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 hair_color, vec3 abs_coef, bool using_karis) {
    if (using_karis) return NP_TT_K(phi, theta_d, eta_perpendic, eta_parallel, hair_color);
    return  NP(1, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
}

vec3 get_TRT(float phi, float theta_d, float eta_perpendic, float eta_parallel, float long_width_R, vec3 hair_color, vec3 abs_coef, bool using_karis) {
    if (using_karis) return NP_TRT_K(phi, theta_d, eta_perpendic, eta_parallel, hair_color, long_width_R);
     return NP_TRT(phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
}

vec3 marschener(vec3 tangent, vec3 viewDirection, vec3 lightDirection, vec3 light_bulb_color, vec3 abs_coef, vec3 hair_color){
    float long_width_R = radians(longitudinal_width);     //Beta R im pretty sure these are fine in degrees -- spoiler allert they werent
    float long_width_TT = long_width_R / 2;      //Beta TT
    float long_width_TRT = long_width_R * 2;     //Beta TRT

    float long_shift_R = radians(longitudinal_shift);     //Alpha R
    float long_shift_TT = -long_shift_R / 2;     //Alpha TT
    float long_shift_TRT = -3 * long_shift_R / 2;//Alpha TRT

    //the projection of the light/view direction onto the normal of the normal plane (the tangent)
    float dotLightTangent = dot(lightDirection, tangent);
    float dotViewTangent = dot(viewDirection, tangent);
    vec3 LT_P = normalize(lightDirection - tangent * dotLightTangent);
    vec3 VT_P = normalize(viewDirection - tangent * dotViewTangent);

    vec3 Normal = normalize(LT_P + VT_P);
    vec3 Binormal = normalize(cross(Normal, tangent));

    float theta_r = acos(dot(VT_P, viewDirection)); //the angle of the view direction and the normal plane (0 = perpendicular to the hair, 90 = tangent) but in radians
    float theta_i = acos(dot(LT_P, lightDirection));//the angle of the light direction and the normal plane

    float theta_h = (theta_r + theta_i) / 2; //half angle
    float theta_d = (theta_r - theta_i) / 2; //difference angle

    float phi_r = acos(dot(Binormal, VT_P)); //the angle of the view direction and the binormal (0 = normal, 90 = binormal)
    float phi_i = acos(dot(Binormal, LT_P)); //the angle of the light direction and the binormal
    float phi = phi_r - phi_i; //the relative azimuth
    

    float eta_perpendic = B_index(theta_d, refraction_index);
    float eta_parallel =  C_index(theta_d, refraction_index);

    float M_R   = normal_pdf(theta_h, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(theta_h, long_shift_TT, long_width_TT);
    float M_TRT = normal_pdf(theta_h, long_shift_TRT, long_width_TRT);

    vec3 N_R = vec3(0), N_TT = vec3(0), N_TRT = vec3(0);
    if(karis_mode == 1)
    {
        if(enable_r == 1)  N_R   = NP_R_K(phi, theta_d, refraction_index, eta_parallel, dot(lightDirection, viewDirection));
        if(enable_tt == 1) N_TT  = NP_TT_K(phi, theta_d, eta_perpendic, eta_parallel, hair_color);
        if(enable_trt == 1)N_TRT = NP_TRT_K(phi, theta_d, eta_perpendic, eta_parallel, hair_color, long_width_R);

    }
    else
    {
        if(enable_r == 1){
            N_R = NP(0, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
            //N_R   = NP_R(phi, theta_d, eta_perpendic, eta_parallel);
        }


        if(enable_tt == 1) N_TT  = NP(1, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);

        if(enable_trt == 1){
            //N_TRT = NP(2, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
            N_TRT = NP_TRT(phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
        }
    }


    float cos_squared_theta_d = pow(cos(theta_d), 2);

    vec3 spec = M_R * N_R / cos_squared_theta_d + M_TT * N_TT / cos_squared_theta_d + M_TRT * N_TRT / cos_squared_theta_d;
    //vec3 spec = M_TRT * N_TRT / cos_squared_theta_d;

    return light_bulb_color * spec;
    //return vec3(phi_i / (2 * M_PI));
}



vec3 marschener(float phi, float theta_i, float theta_r, vec3 light_bulb_color, vec3 abs_coef, vec3 hair_color){
    float theta_h = (theta_r + theta_i) / 2; //half angle
    float theta_d = (theta_r - theta_i) / 2; //difference angle

    float eta_perpendic = B_index(theta_d, refraction_index);
    float eta_parallel =  C_index(theta_d, refraction_index);

    float long_width_R = radians(longitudinal_width);     //Beta R im pretty sure these are fine in degrees -- spoiler allert they werent
    float long_shift_R = radians(longitudinal_shift);     //Alpha R

    float M_R   = normal_pdf(theta_h, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(theta_h, -long_shift_R / 2, long_width_R / 2);
    float M_TRT = normal_pdf(theta_h, -3 * long_shift_R / 2, long_width_R * 2);

    vec3 N_R = vec3(0), N_TT = vec3(0), N_TRT = vec3(0);

    if(enable_r == 1){
        N_R = NP(0, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
        //N_R   = NP_R(phi, theta_d, eta_perpendic, eta_parallel);
    }

    if(enable_tt == 1) N_TT  = NP(1, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);

    if(enable_trt == 1){
        //N_TRT = NP(2, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
        N_TRT = NP_TRT(phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
    }

    float cos_squared_theta_d = pow(cos(theta_d), 2);

    vec3 spec = M_R * N_R / cos_squared_theta_d + M_TT * N_TT / cos_squared_theta_d + M_TRT * N_TRT / cos_squared_theta_d;

    return light_bulb_color * spec;
}
#endif