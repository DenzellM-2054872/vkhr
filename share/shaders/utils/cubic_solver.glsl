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

float cbrt(float x) { return sign(x) * pow(abs(x), 1.0 / 3.0); }

vec4 NormalizedCubicSolver(float A, float B, float C, float D) {
    float u = B / (3.0 * A);

    // Depress to x^3 + px + q by substituting x-b/3a
    // This can be found by substituting x+u and solving for u so that the x^2
    // term gets eliminated (then of course dividing by the leading coefficient)
    float p = (C - B * u) / A;
    float q = (D - (C - 2.0 * B * B / (9.0 * A)) * u) / A;

    // Everything blows up when p=0 so give this case special treatment
    if (abs(p) < M_EPS) { return vec4(cbrt(-q) - u, 0, 0, 1); }

    float h = 0.25 * q * q + p * p * p / 27.0;
    if (h > 0.0) { // Check depressed cubic discriminant
        h = sqrt(h);
        float o = -0.5 * q;
        float x = cbrt(o - h) + cbrt(o + h) - u; // Cardano's formula (see https://en.wikipedia.org/wiki/Cubic_equation)
        return vec4(x, 0, 0, 1);
    }

    // Solve by mapping an inverse smoothstep between the critical points
    // I found a whole lot simplified right out so now it probably looks rather obfuscated
    float m = sqrt(-p / 3.0);
    float x = -2.0 * m * sin(asin(1.5 * q / (p * m)) / 3.0);

    // Factor out the root to solve for the rest as a quadratic
    h = sqrt(-3.0 * x * x - 4.0 * p);
    float y = 0.5 * (h - x);
    float z = 0.5 * (-h - x);

    return vec4(x - u, y - u, z - u, 3);
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

vec4 testRoots(float a, float b, float c, float d, vec4 roots)
{

    for (int i = 0; i < roots[3]; i++) {
        if (abs((((a * roots[i] + b) * roots[i]) + c) * roots[i] + d) > 0) {
            return vec4(4);
        }
    }

    return roots;
}


vec4 cubic_solver(float a, float b, float c, float d)
{
    vec4 roots;
    if (abs(a) < M_EPS)
        roots = QuadraticSolver(b, c, d);
    else {
        roots = NormalizedCubicSolver(b / a, c / a, d / a);
       // roots =  NormalizedCubicSolver(a, b, c, d);

    }
    return roots;
    //return testRoots(a, b, c, d, roots);
}
#endif