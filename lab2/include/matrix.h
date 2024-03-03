#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef struct{
	uint8_t size; //Размерность матрицы
	float **data; //Данные, хранящиеся в матрице
} Matrix;

//Выделяет память под новую матрицу
Matrix *new_matrix(uint8_t _size);

//Создает копию матрицы
Matrix *copy_matrix(Matrix *m);

//Возвращает единичную матрицу
Matrix *identity_matrix(uint8_t _size);

//Освобождает память из-под матрицы
void free_matrix(Matrix *m);

//Выводит данные матрицы на экран
void print_matrix(Matrix *m);

//Возвращает транспонированную матрицу
Matrix *transposed_matrix(Matrix *m);

//Диагональ матрицы
Matrix *matrix_diagonal(Matrix *m);

//Нижний треугольник матрицы
Matrix *lower_triangle(Matrix *m);

//Верхний треугольник матрицы
Matrix *upper_triangle(Matrix *m);

//Возвращает матрицу равную произведению двух матриц
Matrix *matrix_multiply(Matrix *m1, Matrix *m2);

//Умножает матрицу на константу
Matrix *multiply_on_number(const float num, Matrix *m);

//Вычитает из первой матрицы вторую, возвращает результат вычитания
Matrix *matrix_subtruct(Matrix *m1, Matrix *m2);

//Складывает матрицы и возвращает результат сложения
Matrix *matrix_add(Matrix *m1, Matrix *m2);

//Считает норму матрицы первого типа
float matrix_norm1(const float max_eigenvalue);

//Считает норму матрицы второго типа
float matrix_norm2(Matrix *m);

//Считает норму матрицы третьего типа
float matrix_norm3(Matrix *m);

//Число обусловленности с нормой первого типа
float condition_number1(const float eigen1, const float eigen2);

//Число обусловленности с нормой второго типа
float condition_number2(Matrix *m);

//Число обусловленности с нормой третьего типа
float condition_number3(Matrix *m);

//Определитель матрицы
float determinant(Matrix *m);

//Обратная матрица
Matrix *inversed_matrix(Matrix *m);
