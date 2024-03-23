#version 460 core
#extension GL_EXT_debug_printf : enable
#include "strand.glsl"
#include "../scene_graph/params.glsl"
#include "../scene_graph/camera.glsl"
#include "../scene_graph/lights.glsl"

#define M_PI 3.1415926535897932384626433832795
#define M_EPS 0.00001

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

    return (1 / abs(s)) * inv_sqrt_2pi * exp(- (a * a) / (2 * s * s));
}

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

vec4 LinearSolver(float a, float b)
{
    vec4 roots = vec4(0);

    if (abs(a) > M_EPS)
    {
        roots[0] = -b / a;
        roots[3] = 1;
    }

    return roots;
}

vec4 QuadraticSolver(float a, float b, float c)
{
    vec4 roots;
    if (abs(a) < M_EPS)
        return LinearSolver(b, c);
    else
    {
        roots = vec4(0);

        float D = b * b - 4 * a * c;

        if (abs(D) < M_EPS)
        {
            roots[0] = -b / (2 * a);
            roots[1] = -b / (2 * a);
            roots[3] = 2;
        }
        else if (D > 0)
        {
            float delta = sqrt(D);
            roots[0] = (-b + delta) / (2 * a);
            roots[1] = (-b - delta) / (2 * a);
            roots[3] = 2;
        }
    }

    return roots;
}

vec4 NormalizedCubicSolver(float A, float B, float C)
{
    vec4 roots;

    if (abs(C) < M_EPS)	//	x = 0 solution
    {
        roots = QuadraticSolver(1, A, B);
		roots[int(roots[3])] = 0;
		roots[3]++;
    }
    else
    {
        roots = vec4(0);

        float Q = (3 * B - A * A) / 9;
        float R = (9 * A * B - 27 * C - 2 * A * A * A) / 54;
        float D = Q * Q * Q + R * R;

        if (D > 0)	// 1 root
        {
            float sqrtD = sqrt(D);
            float s = sign(R + sqrtD) * pow(abs(R + sqrtD), 1.0f / 3.0f);
            float t = sign(R - sqrtD) * pow(abs(R - sqrtD), 1.0f / 3.0f);

            roots[0] = (-A / 3 + (s + t));
            roots[3] = 1;
        }
        else	// 3 roots
        {
            float theta = acos(R / sqrt(-(Q * Q * Q)));
            float sqrtQ = sqrt(-Q);
            roots[0] = (2 * sqrtQ * cos(theta / 3) - A / 3);
            roots[1] = (2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - A / 3);
            roots[2] = (2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - A / 3);
            roots[3] = 3;
        }
    }

    return roots;
}
vec4 cubic_solver(float a, float b, float c, float d)
{
    if (abs(a) < M_EPS)
        return QuadraticSolver(b, c, d);
    else
        return NormalizedCubicSolver(b / a, c / a, d / a);
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
    vec3 abs_coef = vec3(abs_coef_R, abs_coef_G, abs_coef_B);
 
    float M_R   = normal_pdf(theta_h, long_shift_R, long_width_R);
    float M_TT  = normal_pdf(theta_h, long_shift_TT, long_width_TT);
    float M_TRT = normal_pdf(theta_h, long_shift_TRT, long_width_TRT);
    
    vec3 N_R = NP(0, phi, theta_d, refraction_index, abs_coef);
    vec3 N_TT = NP(1, phi, theta_d, refraction_index, abs_coef);
    vec3 N_TRT = NP(2, phi, theta_d, refraction_index, abs_coef);

    float cos_squared_theta_d = pow(cos(theta_d), 2);

    float cosTL_squared = dotLightTangent*dotLightTangent;
    float one_mietas_cosTL_squared = 1.0f - cosTL_squared;
    float sinTL = sqrt(one_mietas_cosTL_squared);
    vec3 diffuse_colors  = hair_color  * sinTL;

    vec3 light_bulb_color = lights[0].intensity;
    vec3 spec = M_R * N_R / cos_squared_theta_d + M_TT * N_TT / cos_squared_theta_d + M_TRT * N_TRT / cos_squared_theta_d;

    vec3 specular_colors = light_bulb_color * spec;

    color = vec4(specular_colors + diffuse_colors, 1.0f);
}

