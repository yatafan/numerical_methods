#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#define N 2 // Размерность системы решений

// Одномерный метод Ньютона
float newtons_method_1d(float x0, 							// Начальное приближение
						float (*func)(float), 				// Функиця f(x)
						float (*derivative_func)(float), 	// Производная функции
						float accuracy, 					// Точность
						uint8_t *iterations); 				// Количество итераций

// Двумерный метод Ньютона
float *newtons_method_2d(float x0, 						// Начальное приближение x
						float y0, 						// Начальное приближение y
						float (*f)(float, float), 		// f(x,y)
						float (*g)(float, float), 		// g(x,y)
						float (*dfdx)(float, float), 	// Частная производная f по x
						float (*dfdy)(float, float), 	// Частная производная f по y
						float (*dgdx)(float, float), 	// Частная производная g по x
						float (*dgdy)(float, float), 	// Частная производная g по y
						float accuracy, 				// Точность
						uint8_t *iterations); 			// Количество итераций
