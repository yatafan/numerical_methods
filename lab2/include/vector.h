#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef struct{
	uint8_t size;
	float *data;
} Vector;

//Выделяет память под новый вектор
Vector *new_vector(uint8_t _size);

//Возвращает копию вектора
Vector *copy_vector(Vector *v);

//Освобождает память из-под вектора
void free_vector(Vector *v);

//Выводит вектор в терминал
void print_vector(Vector *v);

//Умножает вектор на константу, возвращает результат умножения
Vector *multiply_vector_by_number(const float num, Vector *v);

//Складывает два вектора, возвращает результат сложения
Vector *add_vectors(Vector *v1, Vector *v2);

//Вычитает два вектора, возвращает результат вычитания
Vector *subtruct_vectors(Vector *v1, Vector *v2);

//Находит Евклидову норму вектора
float euclidean_norm(Vector *v);
