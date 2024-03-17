#include "../include/vector.h"

Vector *new_vector(uint8_t _size){
	Vector *vec = malloc(sizeof(Vector));

	vec->size = _size;
	vec->data = malloc(_size * sizeof(float));
	
	return vec;
}

Vector *copy_vector(Vector *v){
	Vector *copy_v = new_vector(v->size);
	
	for(uint8_t i = 0; i < v->size; i++){
		copy_v->data[i] = v->data[i];
	}

	return copy_v;
}

void free_vector(Vector *v){
	free(v->data);
	free(v);
}

void print_vector(Vector *v){
	for(uint8_t i = 0; i < v->size; i++){
		printf("%f\t", v->data[i]);
	}
	printf("\n");
}

Vector *multiply_vector_by_number(const float num, Vector *v){
	Vector *result = new_vector(v->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		result->data[i] = num * v->data[i];
	}

	return result;
}

Vector *add_vectors(Vector *v1, Vector *v2){
	Vector *result = new_vector(v1->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		result->data[i] = v1->data[i] + v2->data[i];
	}

	return result;
}

Vector *subtruct_vectors(Vector *v1, Vector *v2){
	Vector *result = new_vector(v1->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		result->data[i] = v1->data[i] - v2->data[i];
	}

	return result;
}

float euclidean_norm(Vector *v){
	float norm = 0;
	
	for(uint8_t i = 0; i < v->size; i++){
		norm += v->data[i] * v->data[i];
	}

	norm = sqrtf(norm);

	return norm;
}
