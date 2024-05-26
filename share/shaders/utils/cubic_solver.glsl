#ifndef VKHR_CUBIC_SOLVER_GLSL
#define VKHR_CUBIC_SOLVER_GLSL
#include "math.glsl"

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
    vec4 roots = vec4(0);
    if (abs(a) < M_EPS)
        return LinearSolver(b, c);
    else
    {
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
#endif