#include "../include/slau.h"

Vector *multiply_matrix_by_vector(Matrix *m, Vector *v){
	Vector *result = new_vector(v->size);
	
	for(uint8_t i = 0; i < result->size; i++){
		float sum = 0;
		for(uint8_t j = 0; j < result->size; j++){
			sum += m->data[i][j] * v->data[j];
		}

		result->data[i] = sum;
		sum = 0;
	}

	return result;
}

Vector *simple_iterations(Matrix *m, Vector *x0, Vector *b, const float t, const float accuracy, unsigned int *iterations){
	Matrix *E = identity_matrix(m->size);
	Matrix *tA = multiply_on_number(t, m);
	Matrix *main_matrix = matrix_subtruct(E, tA);
	
	Vector *tb = multiply_vector_by_number(t, b);

	Vector *X = new_vector(x0->size);
	Vector *X0 = copy_vector(x0);
	
	//printf("матрица E-tA:\n");
	//print_matrix(main_matrix);

	bool terminate = false;

	while(!terminate){
		Vector *first_term = multiply_matrix_by_vector(main_matrix, X0); //Первое слагаемое
		free_vector(X);
		X = add_vectors(first_term, tb);
		free_vector(first_term);

		Vector *accuracy_vector = subtruct_vectors(X, X0); //Для проверки точности
		if(euclidean_norm(accuracy_vector) < accuracy){ terminate = true; }
		free_vector(accuracy_vector);

		free_vector(X0);
		X0 = copy_vector(X);

		(*iterations)++;
	}

	free_matrix(E);
	free_matrix(tA);
	free_matrix(main_matrix);

	free_vector(tb);
	free_vector(X0);
	return X;
}

Vector *jacobi_method(Matrix *m, Vector *x0, Vector *b, const float accuracy, unsigned int *iterations){
	Matrix *L = lower_triangle(m);
	Matrix *U = upper_triangle(m);
	Matrix *D = matrix_diagonal(m);
	
	Vector *X = new_vector(x0->size);
	Vector *X0 = copy_vector(x0);

	//Нахождение обратной матрицы методом Гаусса
	Matrix *inv_D = identity_matrix(m->size);

	for(uint8_t i = 0; i < D->size; i++){
		inv_D->data[i][i] /= D->data[i][i];
	}

	Matrix *LU = matrix_add(L, U); // L + U
	
	bool terminate = false;

	while(!terminate){
		free_vector(X);
		Vector *term = multiply_matrix_by_vector(LU, X0); //Слагаемое
		Vector *multiplier = subtruct_vectors(b, term); //Множитель

		X = multiply_matrix_by_vector(inv_D, multiplier);

		free_vector(term);
		free_vector(multiplier);

		Vector *accuracy_vector = subtruct_vectors(X, X0);
		if(euclidean_norm(accuracy_vector) < accuracy){ terminate = true; }
		free_vector(accuracy_vector);
		
		free_vector(X0);
		X0 = copy_vector(X);

		(*iterations)++;
	}

	free_vector(X0);

	free_matrix(L);
	free_matrix(U);
	free_matrix(D);
	free_matrix(inv_D);
	free_matrix(LU);
	return X;
}

Vector *gauss_seidel_method(Matrix *m, Vector *x0, Vector *b, const float accuracy, unsigned int *iterations){
	Matrix *L = lower_triangle(m);
	Matrix *U = upper_triangle(m);
	Matrix *D = matrix_diagonal(m);
	
	Vector *X = new_vector(x0->size);
	Vector *X0 = copy_vector(x0);
	
	Matrix *LD = matrix_add(L, D); // L + D
	Matrix *inv_LD = identity_matrix(LD->size);

	//Нахождение обратной матрицы методом Гаусса
	for(uint8_t i = 0; i < LD->size; i++){
		float main_elem = LD->data[i][i];
		
		for(uint8_t j = i + 1; j < LD->size; j++){
			float koef = main_elem / LD->data[j][i];
			
			//Приравнивание коэффициентов
			for(uint8_t k = 0; k < LD->size; k++){
				LD->data[j][k] *= koef;
				inv_LD->data[j][k] *= koef;
			}
			
			//Вычитание строк
			for(uint8_t k = 0; k < LD->size; k++){
				LD->data[j][k] -= LD->data[i][k];
				inv_LD->data[j][k] -= inv_LD->data[i][k];
			}
		}
		
		//Приведение к единице
		for(uint8_t j = 0; j < LD->size; j++){
			LD->data[i][j] /= main_elem;
			inv_LD->data[i][j] /= main_elem;
		}
	}
	
	Vector *term = multiply_matrix_by_vector(inv_LD, b); //Слагаемое
	Matrix *multiplier = matrix_multiply(inv_LD, U); //Множитель
	
	bool terminate = false;

	while(!terminate){
		free_vector(X);
		Vector *term2 = multiply_matrix_by_vector(multiplier, X0);

		X = subtruct_vectors(term, term2);
		
		free_vector(term2);

		Vector *accuracy_vector = subtruct_vectors(X, X0);
		if(euclidean_norm(accuracy_vector) < accuracy){ terminate = true; }
		free_vector(accuracy_vector);

		free_vector(X0);
		X0 = copy_vector(X);

		(*iterations)++;
	}
	

	free_vector(X0);
	free_vector(term);
	
	free_matrix(multiplier);
	free_matrix(inv_LD);
	free_matrix(LD);
	free_matrix(L);
	free_matrix(U);
	free_matrix(D);
	return X;
}
