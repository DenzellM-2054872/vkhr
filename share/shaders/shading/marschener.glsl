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

float fresnel2(float eta, float x)
{
    float F_0 = pow(1 - eta, 2) / pow(1 + eta, 2);
    return F_0 + (1 - F_0) * pow(1 - x, 5);
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
    if(eta == 1.55) return 1.19 / cos(theta) + 0.36 * cos(theta);

    float sin_theta = sin(theta);
    return sqrt(eta * eta - sin_theta * sin_theta) / cos(theta);
}


float inv_first_der(int p, float eta_perpendic, float h)
{
    float gamma_i = asin(h);
    float c = asin( 1 / eta_perpendic );
    float d_gamma = ((6 * p * c / M_PI) - 2) - 
		3 * 8 * (p * c / (M_PI * M_PI * M_PI)) * gamma_i * gamma_i;

    return sqrt(1 - h * h) / d_gamma;
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

    for (int i = 0; i < roots[3]; i++ )
    {
        float gamma_i = roots[i];
        float h = sin(gamma_i);
        //if(h > M_EPS) continue;

        vec3 finalAbsorption = attenuation(absorb, p, h, eta_perpendic, eta_parallel, theta_d);
        float inverseDerivateAngle = inv_first_der(p, eta_perpendic, h);
        res += finalAbsorption * 2 * abs(inverseDerivateAngle); //0.5 here
    }
    vec3 final_res = vec3(min(res.r, 1), min(res.g, 1), min(res.b, 1));
    return final_res;
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
    //this is how it is described in the frostbite paper
    //at certain angles this stops working mb figure out why later
    // why does this work with -abs but not normal abs
    vec3 res = vec3(0.25 * cos(phi / 2));
    res *= vec3(fresnel(eta_perpendic, eta_parallel, sqrt(-abs(0.5 * (radians(1) + view_light_angle)))));

    return vec3(min(1, res.r),min(1, res.g), min(1, res.b));
}

float calc_h(float phi, float inv_eta_tick)
{
    float top = sign(phi) * cos(phi / 2);
    float bottom = sqrt(1 + pow(inv_eta_tick, 2) - 2 * inv_eta_tick * sign(phi) * sin(phi / 2));

    return abs(top / bottom);
}

vec3 NP_TT_K(float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 absorb)
{
    float a = 1/eta_perpendic;
    float h = calc_h(phi, a);
    //float h = (1 + a * (0.6 - 0.8 * cos(phi))) * cos(phi / 2)

    vec3 finalAbsorption = attenuation(absorb, 1, h, eta_perpendic, eta_parallel, theta_d);
    
    //from the karis paper (quite difrent then the original)
    float distrib = exp(-3.65 * cos(phi) - 3.98);
    vec3 res = finalAbsorption * distrib;

    vec3 final_res = vec3(min(res.r, 1), min(res.g, 1), min(res.b, 1));
    return final_res;
}

vec3 NP_TRT(float phi, float theta_d, float eta_perpendic, float eta_parallel, vec3 absorb)
{
    float delta_hM = caustic_intensity_limit;
    float w_c = caustic_width * M_PI / 180;
    float k_G = glint_scale_fac;
    float delta_eta_tick = caustic_merge_range;

    float h_c = 0;
    float delta_h, t, phi_c;

    if (eta_perpendic < 2)
    {
        float c = asin(1 / eta_perpendic);
        phi_c = sqrt((6 * 2 * c / M_PI - 2) / 
			(3 * 8 * (2 * c / (M_PI * M_PI * M_PI))));
        h_c = abs(sin(phi_c));

        float inv_der_angle = inv_second_der(2, eta_perpendic, h_c);

        delta_h = min(delta_hM, 2 * sqrt(2 * w_c * inv_der_angle));
        t = 1;
    } else {

        phi_c = 0;
        delta_h = delta_hM;
        t = 1 - smoothstep(2, 2 + delta_eta_tick, eta_perpendic);
    }

    phi_c = angle_polynomial(2, eta_perpendic, h_c);
    vec3 res = NP(2, phi, theta_d, eta_perpendic, eta_parallel, absorb);
    vec3 final_abs = attenuation(absorb, 2, h_c, eta_perpendic, eta_parallel, theta_d);

    res = res * (1 - t * normal_pdf(phi, phi_c, w_c) / normal_pdf(0, 0, w_c));
    res = res * (1 - t * normal_pdf(phi, -phi_c, w_c) / normal_pdf(0, 0, w_c));
    res = res + t * k_G * final_abs * delta_h * (normal_pdf(phi, phi_c, w_c) + normal_pdf(phi, -phi_c, w_c));

    return vec3(min(1, res.r),min(1, res.g), min(1, res.b));
}


vec3 marschener(vec3 tangent, vec3 viewDirection, vec3 lightDirection, vec3 light_bulb_color, vec3 abs_coef){
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
    //the karis presentation implies this to be the case
    //float phi = dot(viewDirection, lightDirection); //the relative azimuth

    float eta_perpendic = B_index(theta_d, refraction_index);
    float eta_parallel = (refraction_index * refraction_index) / eta_perpendic;

    float M_R   = normal_pdf(theta_h, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(theta_h, -long_shift_TT, long_width_TT);
    float M_TRT = normal_pdf(theta_h, -long_shift_TRT, long_width_TRT);
    vec3 N_R = vec3(0), N_TT = vec3(0), N_TRT = vec3(0);
    if(karis_mode == 1)
    {
        if(enable_r == 1)  N_R   = NP_R_K(phi, theta_d, eta_perpendic, eta_parallel, dot(viewDirection, lightDirection));
        if(enable_tt == 1) N_TT  = NP_TT_K(phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
        if(enable_trt == 1)N_TRT = NP_TRT(phi, theta_d, eta_perpendic, eta_parallel, abs_coef);

    }
    else
    {
        if(enable_r == 1){
            // N_R = NP(0, phi, theta_d, refraction_index, abs_coef);
            N_R   = NP_R(phi, theta_d, eta_perpendic, eta_parallel);
        }


        if(enable_tt == 1) N_TT  = NP(1, phi, theta_d, eta_perpendic, eta_parallel, abs_coef);

        if(enable_trt == 1){
            //vec3 N_TRT = NP(2, phi, theta_d, refraction_index, abs_coef);
            N_TRT = NP_TRT(phi, theta_d, eta_perpendic, eta_parallel, abs_coef);
        }
    }


    float cos_squared_theta_d = pow(cos(theta_d), 2);

    vec3 spec = M_R * N_R / cos_squared_theta_d + M_TT * N_TT / cos_squared_theta_d + M_TRT * N_TRT / cos_squared_theta_d;
    //vec3 spec = M_TRT * N_TRT / cos_squared_theta_d;

    return light_bulb_color * spec;
    //return vec3(phi_i / (2 * M_PI));
}
#endif