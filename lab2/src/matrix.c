#include "../include/matrix.h"

Matrix *new_matrix(uint8_t _size){
	Matrix *m = malloc(sizeof(Matrix));

	m->size = _size;
	m->data = malloc(_size * sizeof(float*));

	for(uint8_t i = 0; i < _size; i++){
		m->data[i] = malloc(_size * sizeof(float));
	}

	return m;
}

Matrix *copy_matrix(Matrix *m){
	Matrix *copy_m = new_matrix(m->size);

	copy_m->size = m->size;

	for(uint8_t i = 0; i < m->size; i++){
		for(uint8_t j = 0; j < m->size; j++){
			copy_m->data[i][j] = m->data[i][j];
		}
	}

	return copy_m;
}

Matrix *identity_matrix(uint8_t _size){
	Matrix *m = new_matrix(_size);
	
	for(uint8_t i = 0; i < _size; i++){
		for(uint8_t j = 0; j < _size; j++){
			m->data[i][j] = (i == j); // 0 или 1
		}
	}

	return m;
}

void free_matrix(Matrix *m){
	for(uint8_t i = 0; i < m->size; i++){
		free(m->data[i]);
	}
	free(m->data);
	free(m);
}

void print_matrix(Matrix *m){
	for(uint8_t i = 0; i < m->size; i++){
		for(uint8_t j = 0; j < m->size; j++){
			printf("%f\t", m->data[i][j]);
		}
		printf("\n");
	}
}

Matrix *transposed_matrix(Matrix *m){
	Matrix *transp_matrix = new_matrix(m->size);

	for(uint8_t i = 0; i < m->size; i++){
		for(uint8_t j = 0; j < m->size; j++){
			transp_matrix->data[i][j] = m->data[j][i];
		}
	}

	return transp_matrix;
}

Matrix *matrix_diagonal(Matrix *m){
	Matrix *result = new_matrix(m->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		for(uint8_t j = 0; j < result->size; j++){
			if(i == j){
				result->data[i][j] = m->data[i][j];
			}
			else{
				result->data[i][j] = 0.0;
			}
		}
	}

	return result;
}

Matrix *lower_triangle(Matrix *m){
	Matrix *result = new_matrix(m->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		for(uint8_t j = 0; j < result->size; j++){
			if(j < i){
				result->data[i][j] = m->data[i][j];
			}
			else{
				result->data[i][j] = 0.0;
			}
		}
	}

	return result;
}

Matrix *upper_triangle(Matrix *m){
	Matrix *result = new_matrix(m->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		for(uint8_t j = 0; j < result->size; j++){
			if(j > i){
				result->data[i][j] = m->data[i][j];
			}
			else{
				result->data[i][j] = 0.0;
			}
		}
	}

	return result;
}

Matrix *matrix_multiply(Matrix *m1, Matrix *m2){
	Matrix *result = new_matrix(m1->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		float sum = 0;
		for(uint8_t j = 0; j < result->size; j++){
			for(uint8_t k = 0; k < result->size; k++){
				sum += m1->data[i][k] * m2->data[k][j];
			}

			result->data[i][j] = sum;
			sum = 0;
		}
	}

	return result;
}

Matrix *multiply_on_number(const float num, Matrix *m){
	Matrix *result = new_matrix(m->size);

	for(uint8_t i = 0; i < result->size; i++){
		for(uint8_t j = 0; j < result->size; j++){
			result->data[i][j] = num * m->data[i][j];
		}
	}

	return result;
}

Matrix *matrix_subtruct(Matrix *m1, Matrix *m2){
	Matrix *result = new_matrix(m1->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		for(uint8_t j = 0; j < result->size; j++){
			result->data[i][j] = m1->data[i][j] - m2->data[i][j];
		}
	}

	return result;
}

Matrix *matrix_add(Matrix *m1, Matrix *m2){
	Matrix *result = new_matrix(m1->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		for(uint8_t j = 0; j < result->size; j++){
			result->data[i][j] = m1->data[i][j] + m2->data[i][j];
		}
	}

	return result;
}

float matrix_norm1(const float max_eigenvalue){
	return sqrtf(max_eigenvalue);
}

float matrix_norm2(Matrix *m){
	float norm = 0;
	float sum = 0;
	
	for(uint8_t i = 0; i < m->size; i++){
		for(uint8_t j = 0; j < m->size; j++){
			sum += fabsf(m->data[i][j]);
		}
		if(sum > norm){ norm = sum; }
		sum = 0;
	}
	
	return norm;
}

float matrix_norm3(Matrix *m){
	float norm = 0;
	float sum = 0;
	
	for(uint8_t i = 0; i < m->size; i++){
		for(uint8_t j = 0; j < m->size; j++){
			sum += fabsf(m->data[j][i]);
		}
		if(sum > norm){ norm = sum; }
		sum = 0;
	}
	
	return norm;
}

float condition_number1(const float eigen1, const float eigen2){
	return matrix_norm1(eigen1) * matrix_norm1(eigen2);
}

float condition_number2(Matrix *m){
	Matrix *inv_m = inversed_matrix(m);
	
	float cond_number = matrix_norm2(m) * matrix_norm2(inv_m);

	free_matrix(inv_m);
	return cond_number;
}

float condition_number3(Matrix *m){
	Matrix *inv_m = inversed_matrix(m);
	
	float cond_number = matrix_norm3(m) * matrix_norm3(inv_m);

	free_matrix(inv_m);
	return cond_number;
}

float determinant(Matrix *m){
	Matrix *copy_m = copy_matrix(m);
	float det;
	
	if(copy_m->size == 2){
		det = (copy_m->data[0][0] * copy_m->data[1][1]) - (copy_m->data[0][1] * copy_m->data[1][0]);
	}
	else{
		float first_elem = copy_m->data[0][0];

		for(uint8_t i = 1; i < copy_m->size; i++){
			float main_elem = copy_m->data[i][0];

			for(uint8_t j = 0; j < copy_m->size; j++){
				copy_m->data[i][j] -= copy_m->data[0][j] * (main_elem / first_elem);
			}
		}

		Matrix *minor = new_matrix(copy_m->size - 1);
		
		for(uint8_t i = 1; i < copy_m->size; i++){
			for(uint8_t j = 1; j < copy_m->size; j++){
				minor->data[i-1][j-1] = copy_m->data[i][j];
			}
		}

		det = first_elem * determinant(minor);

		free_matrix(minor);
	}

	free_matrix(copy_m);
	return det;
}

Matrix *inversed_matrix(Matrix *m){
	float det_m = fabsf(determinant(m));

	if(det_m == 0){ return NULL; }

	Matrix *complements_matrix = new_matrix(m->size);
	
	for(uint8_t i = 0; i < m->size; i++){
		for(uint8_t j = 0; j < m->size; j++){
			float algebraic_complement = powf(-1, (i+j) % 2);
			uint8_t mx, my; //Индексы элементов в миноре

			Matrix *minor = new_matrix(m->size - 1);
			
			for(uint8_t x = 0; x < m->size; x++){
				if(x == i){ continue; }
				mx = x;
				if(x > i){ mx -= 1; }

				for(uint8_t y = 0; y < m->size; y++){
					if(y == j){ continue; }

					my = y;
					if(y > j){ my -= 1; }

					minor->data[mx][my] = m->data[x][y];
				}
			}

			complements_matrix->data[i][j] = algebraic_complement * determinant(minor);

			free_matrix(minor);
		}
	}

	Matrix *inversed_m = transposed_matrix(complements_matrix);
	free_matrix(complements_matrix);

	for(uint8_t i = 0; i < inversed_m->size; i++){
		for(uint8_t j = 0; j < inversed_m->size; j++){
			inversed_m->data[i][j] /= det_m;
		}
	}

	return inversed_m;
}
