#include "matrix.h"
#include "vector.h"

#include <stdbool.h>

//Умножает матрицу на вектор
Vector *multiply_matrix_by_vector(Matrix *m, Vector *v);

//Метод простых итераций
Vector *simple_iterations(Matrix *m, Vector *x0, Vector *b, const float t, const float accuracy, unsigned int *iterations);

//Метод Якоби
Vector *jacobi_method(Matrix *m, Vector *x0, Vector *b, const float accuracy, unsigned int *iterations);

//Метод Гаусса-Зейделя
Vector *gauss_seidel_method(Matrix *m, Vector *x0, Vector *b, const float accuracy, unsigned int *iterations);
