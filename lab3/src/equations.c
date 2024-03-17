#include "../include/equations.h"

float newtons_method_1d(float x0, 
						float (*func)(float), 
						float (*derivative_func)(float), 
						float accuracy,
						uint8_t *iterations)
{
	bool terminate = false;
	float x; // Искомое решение

	while(!terminate){
		x = x0 - (func(x0) / derivative_func(x0));

		if(fabsf(x - x0) < accuracy){
			terminate = true;
		}

		x0 = x;

		(*iterations)++;
	}
	
	return x;
}

float *newtons_method_2d(float x0,
						float y0,
						float (*f)(float, float),
						float (*g)(float, float),
						float (*dfdx)(float, float),
						float (*dfdy)(float, float),
						float (*dgdx)(float, float),
						float (*dgdy)(float, float),
						float accuracy,
						uint8_t *iterations)
{
	float *solution = malloc(N * sizeof(float));
	
	bool terminate = false;

	float denominator; // Знаменатель дроби
	float b1k, b2k;

	while(!terminate){
		denominator = dfdx(x0, y0) * dgdy(x0, y0) - dgdx(x0, y0) * dfdy(x0, y0);
		
		b1k = -1.0 * f(x0, y0) + dfdx(x0, y0) * x0 + dfdy(x0, y0) * y0;
		b2k = -1.0 * g(x0, y0) + dgdx(x0, y0) * x0 + dgdy(x0, y0) * y0;

		solution[0] = (b1k * dgdy(x0, y0) - b2k * dfdy(x0, y0)) / denominator; 	// x
		solution[1] = (b2k * dfdx(x0, y0) - b1k * dgdx(x0, y0)) / denominator; 	// y
	
		float current_accuracy = sqrtf(powf(solution[0] - x0, 2) + powf(solution[1] - y0, 2));

		if(current_accuracy < accuracy){
			terminate = true;
		}

		x0 = solution[0];
		y0 = solution[1];

		(*iterations)++;
	}

	return solution;
} 
