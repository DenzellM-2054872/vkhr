#ifndef VKHR_DUAL_SCATTERING_GLSL
#define VKHR_DUAL_SCATTERING_GLSL

#include "../scene_graph/params.glsl"
#include "../shading/marschener.glsl"
#include "../utils/fresnel.glsl"
#include "../utils/cubic_solver.glsl"
#include "../utils/math.glsl"

vec3 a_f(float theta_d){

	float res = 0;
	float h = M_PI / M_INT_PRES;
	float cos_theta = cos(theta_d);
	for(float omicron = (theta_d - M_PI/2); omicron < theta_d + M_PI / 2; omicron += h){
		for(float phi = -M_PI/2; phi <=  M_PI/2; phi += h){
			//numerical integration mid point rule
			res += marschener(phi + h/2, theta_d, omicron+ h/2, lights[0].intensity ,vec3(abs_coef_R, abs_coef_G, abs_coef_B), hair_color) * h * cos_theta;
		}
	}
	return 1/M_PI * res;
}

vec3 a_b(float theta_d){
	float res = 0;
	float h = M_PI / M_INT_PRES;
	float cos_theta = cos(theta_d);
	for(float omicron = (theta_d + M_PI/2); omicron < theta_d + 3 * M_PI / 2; omicron += h){
		for(float phi = -M_PI/2; phi <=  M_PI/2; phi += h){
			res += marschener(phi + h/2, theta_d, omicron+ h/2, lights[0].intensity ,vec3(abs_coef_R, abs_coef_G, abs_coef_B), hair_color) * h * cos_theta;
		}
	}
	return 1/M_PI * res;
}

vec3 A_1(float theta_d){
	vec3 af2 = vec_pow(a_f(theta_d), 2);

	return (a_b(theta_d) * af2) / (1 - af2);
}

vec3 A_3(float theta_d){
	vec3 af2 = vec_pow(a_f(theta_d), 2);

	return (vec_pow(a_b(theta_d), 3) * af2) / (1 - af2);
}

vec3 A_b(float theta_d){
	return A_1(theta_d) + A_3(theta_d);
}

//this way might be faster depending on how expensive ab & af are
vec3 A_b2(float theta_d){
	vec3 af2 = vec_pow(a_f(theta_d), 2);
	vec3 ab = a_b(theta_d);

	return (ab * af2) / (1 - af2) + (vec_pow(ab), 3) * af2) / (1 - af2);
}


vec3 delta_b(float theta_d, float alpha_b, float alpha_f){
	vec3 af2 = vec_pow(a_f(theta_d), 2);
	vec3 ab2 = vec_pow(a_b(theta_d), 2);

	vec3 tB = alpha_b * (1 - (2 * ab2) / vec_pow(1 - af2, 2));
	vec3 tF = alpha_f * ((2 * vec_pow(1 - af2, 2) + 4 * af2 * ab2) / vec_pow(1 - af2, 3));

	return tB + tF;
}

vec3 sigma_b(float theta_d, float beta_b, float beta_f){
	vec3 af = a_f(theta_d);
	vec3 ab = a_b(theta_d);

	float beta_b2 = beta_b * beta_b;
	float beta_f2 = beta_f * beta_f;

	(1 + 0.7 * vec_pow(af, 2)) * (ab * sqrt(2 * beta_f2 + beta_b2) + vec_pow(ab, 3) * sqrt(2 * beta_f2 + 3 * beta_b2)) / (ab + vec_pow(ab, 3) * (2*beta_f + 3*beta_b));
}

vec3 S_b(float theta_o, float theta_i, float alpha, float beta){
	float theta = (theta_o - theta_i) / 2
	return (1/(M_PI * cos(theta)) * gauss(theta_o + theta_i - delta_b(theta, alpha, alpha), sigma_b(theta, beta, beta));
}

vec3 f_back(float theta_o, float theta_i, float alpha, float beta){
	float theta = (theta_o - theta_i) / 2;

	return 2/cos(theta) * A_b(theta) * S_b(theta_o, theta_i, alpha, beta);

}

vec3 Dual_scat(float phi, float theta_o, float theta_i, float eta_perpendic, float eta_parallel, float view_light_angle, vec3 abs_coef, vec3 hair_color){
	float alpha = radians(longitudinal_shift);
	float beta = radians(longitudinal_width);
	float theta = (theta_o - theta_i) / 2;

	//float fback = f_back(theta_o, theta_i, alpha, beta);
	float fback = 2 * A_b(theta) * gauss(theta + theta_o - delta_b(theta_d, alpha, alpha), )


	float direct_frac = 0; //either 0 or 1 depending on coverage 
	vec3 Fdirect = vec3(0);
	if (direct_frac != 0) {
		float M_R = gauss(theta - alpha, pow(beta, 2));
		float M_TT = gauss(theta + alpha / 2, pow(beta / 2, 2));
		float M_TRT = gauss(theta + 3 * alpha / 2, pow(beta * 2, 2));

		vec3 N_R = get_R(phi, theta, eta_perpendic, eta_parallel, view_light_angle, abs_coef, true);
		vec3 N_TT = get_TT(phi, theta, eta_perpendic, eta_parallel, hair_color, abs_coef, true);
		vec3 N_TRT = get_TRT(phi, theta, eta_perpendic, eta_parallel, hair_color, abs_coef, true);

		vec3 fdirects = M_R * N_R + M_TT * N_TT + M_TRT * N_TRT;

		Fdirect = direct_frac * (fdirects + 0.7 * fback);
	}


	

	
}