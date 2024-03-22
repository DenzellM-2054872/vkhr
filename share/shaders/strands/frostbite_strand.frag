#version 460 core
#extension GL_EXT_debug_printf : enable
#include "strand.glsl"
#include "../scene_graph/params.glsl"
#include "../scene_graph/camera.glsl"
#include "../scene_graph/lights.glsl"

layout(location = 0) in PipelineIn {
    vec4 position;
    vec3 tangent;
    float thickness;
} fs_in;

layout(location = 0) out vec4 color;

float normal_pdf(float x, float m, float s)
{
    const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m);

    return inv_sqrt_2pi / s * exp((a * a) / -2 * (s * s) );
}

float fresnel(float nu, float x)
{
    float F_0 = pow(1 - nu, 2) / pow(1 + nu, 2);

    return F_0 + (1 - F_0) * pow(1 - x, 5);
}

vec3 absorption(int p, float nu, float theta_d, float h)
{
    vec3 abs_coef = vec3(abs_coef_R, abs_coef_G, abs_coef_B);
    float cos_theta = cos(theta_d);
    if(h != 0)
    {
        float sin_theta = sin(theta_d);
        float nu_tick = sqrt(abs(nu * nu) - abs(sin_theta * sin_theta)) / cos_theta;
        float gamma_t = asin(h / nu_tick);
        float gamma_i = asin(h);
        vec3 T = exp(-2 * abs_coef * (1 + cos(2 * gamma_t)));
        vec3 T_pow = vec3(pow(T.x, p), pow(T.y, p), pow(T.z, p));

        vec3 A = T_pow * pow(1 - fresnel(nu, cos_theta * pow(1 - abs(h * h), 0.5)), 2) * pow(fresnel(nu, cos_theta * pow(1 - abs(h * h), 0.5)), p - 1);
        return A;

    }
    vec3 T = exp(-4 * abs_coef);
    vec3 T_pow = vec3(pow(T.x, p), pow(T.y, p), pow(T.z, p));

    vec3 A = T_pow * pow(1 - fresnel(nu, cos_theta), 2) * pow(fresnel(nu, cos_theta), p - 1);
    return A;
}

float distribution(int p, float nu, float theta_d, float h)
{
    float pi = 3.14159265359;
    float cos_theta = cos(theta_d);
    float sin_theta = sin(theta_d);

    float nu_tick = pow(abs(nu * nu) - abs(sin_theta * sin_theta),0.5) / cos_theta;
    float c = asin(1/nu_tick);

    float gamma_i = asin(h);
    float d_phi = ((6 * p * c) / pi - 2) - 3 * (8 * p * c) / pow(pi, 3) * pow(gamma_i, 2);
    float d_h = sqrt(1 - pow(h, 2));
    return d_phi / d_h;
}
void main() {
    float pi = 3.14159265359;

    float long_width_R = longitudinal_width;     //Beta R
    float long_width_TT = long_width_R / 2;      //Beta TT
    float long_width_TRT = long_width_R * 2;     //Beta TRT

    float long_shift_R = longitudinal_shift;     //Alpha R
    float long_shift_TT = -long_shift_R / 2;     //Alpha TT
    float long_shift_TRT = -3 * long_shift_R / 2;//Alpha TRT

    vec3 tangent = normalize(fs_in.tangent);
    vec3 viewDirection = normalize(fs_in.position.xyz - camera.position);
    vec3 lightDirection = normalize(lights[0].origin - fs_in.position.xyz);

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

    /*
    //karis
    float sinIR = sin(theta_i) + sin(theta_r);
    float M_R   = normal_pdf(sinIR, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(sinIR, long_shift_TT, long_width_TT);
    float M_TRT = normal_pdf(sinIR, long_shift_TRT, long_width_TRT);  
    */
    
    float M_R   = normal_pdf(theta_h, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(theta_h, long_shift_TT, long_width_TT);
    float M_TRT = normal_pdf(theta_h, long_shift_TRT, long_width_TRT);
    
    float D_R = 0.25 * cos(0.5 * phi);
    float A_R = fresnel(refraction_index, sqrt(0.5 * (1 + dot(viewDirection, lightDirection))));
    float N_R = D_R * A_R;

    float D_TT = 1 / abs(2 * distribution(1, refraction_index, theta_d, 0));
    vec3 A_TT = absorption(1, refraction_index, theta_d, 0);
    vec3 N_TT = D_TT * A_TT;

    float D_TRT = 1 / abs(2 * distribution(2, refraction_index, theta_d,  0.8660254038)); //big number is root 3, over 2
    vec3 A_TRT = absorption(2, refraction_index, theta_d, 0.8660254038); //big number is root 3, over 2
    vec3 N_TRT = D_TRT * A_TRT;

    float cos_squared_theta_d = pow(cos(theta_d), 2);
    vec3 spec = M_R * N_R / cos_squared_theta_d + M_TT * N_TT / cos_squared_theta_d + M_TRT * N_TRT / cos_squared_theta_d;
    color = vec4(hair_color + spec, 1.0f);
}

