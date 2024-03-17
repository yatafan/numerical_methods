#include <stdio.h>

#include "include/equations.h"
#include "include/vector.h"

#define ACC 0.0001

#define X0 1.0

#define SOLUTIONS 4 //Количество решений системы

float function(float x){
	return powf(5.0, x) - powf(4.0, x) - 9;
}

float derivative_function(float x){
	return powf(5.0, x) * logf(5.0) - powf(4.0, x) * logf(4.0);
}

float f(float x, float y), g(float x, float y);
float dfdx(float x, float y), dfdy(float x, float y), dgdx(float x, float y), dgdy(float x, float y);



int main(){
	uint8_t iterations = 0;
	float x = newtons_method_1d(X0, function, derivative_function, ACC, &iterations);
	
	// Вектора начальных приближений для x и y
	Vector *x0 = new_vector(SOLUTIONS);
	Vector *y0 = new_vector(SOLUTIONS);

	x0->data[0] = 4.0; y0->data[0] = -2.0;
	x0->data[1] = -2.0; y0->data[1] = 4.0;
	x0->data[2] = 2.5; y0->data[2] = 3.5;
	x0->data[3] = 3.5; y0->data[3] = 2.5;
	
	float *system_solution;
	
	printf("\n-=Одномерный метод Ньютона=-\n\n");
	printf("точность: %g\n", ACC);
	printf("начальное приближение:\t%g\n", X0);
	printf("вычисленное решение:\t%g\n", x);
	printf("количество итераций:\t%d\n", iterations);

	iterations = 0;

	printf("\n-=Двумерный метод Ньютона=-\n\n");
	printf("точность: %g\n", ACC);

	for(uint8_t i = 0; i < x0->size; i++){
		printf("начальное приближение: \t%g %g\n", x0->data[i], y0->data[i]);
	
		system_solution = newtons_method_2d(x0->data[i],
											y0->data[i],
											f, g, dfdx, dfdy, dgdx, dgdy,
											ACC,
											&iterations);
	
		printf("вычисленное решение: \t%g %g\n", system_solution[0], system_solution[1]);
		printf("количество итераций: \t%d\n\n", iterations);
	
		free(system_solution);
		iterations = 0;
	}

	free_vector(x0);
	free_vector(y0);
	return 0;
}

float f(float x, float y){
	return powf(x, 3.0) + powf(y, 3.0) - 35;
}

float g(float x, float y){
	return powf(x, 2.0) + powf(y, 2.0) - 13;
}

float dfdx(float x, float y){ return 3.0 * powf(x, 2.0); }

float dfdy(float x, float y){ return 3.0 * powf(y, 2.0); }

float dgdx(float x, float y){ return 2.0 * x; }

float dgdy(float x, float y){ return 2.0 * y; }
