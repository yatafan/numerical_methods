#include "include/slau.h"

#define N 4 //Размерность матрицы

#define EIGEN1 /*825.05843*/ 836.949 //Наибольшее собственное значение матрицы для нормы первого типа
#define EIGEN2 /*0.21669*/ 0.055 //Наибольшее собственное значение обратной матрицы для нормы первого типа

#define T 0.0725 //Тау

#define ACC1 0.01
#define ACC2 0.001
#define ACC3 0.0001

int main(){
	float static_matrix[N][N] = {{10.0, 4.0, 2.0, 3.0},
								{8.0, 12.0, 1.0, 2.0},
								{1.0, 2.0, 6.0, 2.0},
								{1.0, 10.0, 11.0, 22.0}};

/*	float static_matrix[N][N] = {{12.0, 7.0, 4.0, 1.0},
								{5.0, 20.0, 8.0, 7.0},
								{1.0, 4.0, 14.0, 8.0},
								{0.0, 1.0, 4.0, 5.0}};
*/
	Matrix *matr = new_matrix(N);

	for(uint8_t i = 0; i < matr->size; i++){
		for(uint8_t j = 0; j < matr->size; j++){
			matr->data[i][j] = static_matrix[i][j];
		}
	}
	
	//Матрица с количеством итераций для каждого из методов и для каждой точности
	//Строка соотвествует методу, а столбец точности
	//Методы идут в порядке: простых итераций, Якоби, Гаусса-Зейделя
	unsigned int method_iterations[3][3];
	
	//Матрица с разностями норм точного решения и начального приближения
	float difference_norms[3][3];
	Vector *difference; //Разность точного решения и начального приближения

	Vector *b = new_vector(N);

	b->data[0] = 2.0;
	b->data[1] = 5.0;
	b->data[2] = 10.0;
	b->data[3] = 0.0;
/* 	
	b->data[0] = 18.0;
	b->data[1] = 36.0;
	b->data[2] = 26.0;
	b->data[3] = 12.0;
*/
	Vector *x0 = new_vector(N);

	x0->data[0] = 1.0;
	x0->data[1] = 2.0;
	x0->data[2] = 3.0;
	x0->data[3] = 4.0;
 	
/*
	x0->data[0] = 1.0;
	x0->data[1] = 1.0;
	x0->data[2] = 1.0;
	x0->data[3] = 1.0;
*/
	Vector *exact_solution = new_vector(N); //Точное решение системы

	exact_solution->data[0] = -5.0 / 286.0;
	exact_solution->data[1] = 265.0 / 572.0;
	exact_solution->data[2] = 272.0 / 143.0;
	exact_solution->data[3] = -166.0 / 143.0;
/*
	exact_solution->data[0] = 111.0 / 155.0;
	exact_solution->data[1] = 119.0 / 155.0;
	exact_solution->data[2] = 173.0 / 310.0;
	exact_solution->data[3] = 9.0 / 5.0;
*/	
	printf("исходная матрица:\n");
	print_matrix(matr);
	printf("\n");

	printf("число обусловленности (тип 1): %f\n", condition_number1(EIGEN1, EIGEN2));
	printf("число обусловленности (тип 2): %f\n", condition_number2(matr));
	printf("число обусловленности (тип 3): %f\n", condition_number3(matr));
	printf("\n");
	
	unsigned int iterations = 0;
	Vector *X;

	//-----= Метод простых итераций =------	
	printf("-=Метод простых итераций=-\n\n");

	float norm = matrix_norm1(0.96776); //matrix_norm1(0.935);
	
	printf("Сходимость: ");
	if(norm < 1){
		printf("да\n\n");
	}
	else{
		printf("нет\n\n");
	}
	
	printf("Точность %g\n", ACC1);

	X = simple_iterations(matr, x0, b, T, ACC1, &iterations);
	method_iterations[0][0] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[0][0] = euclidean_norm(difference);
	
	printf("Вычисленное решение:\n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	printf("Точность %g\n", ACC2);

	X = simple_iterations(matr, x0, b, T, ACC2, &iterations);
	method_iterations[0][1] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[0][1] = euclidean_norm(difference);

	printf("Вычисленное решение:\n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	printf("Точность %g\n", ACC3);

	X = simple_iterations(matr, x0, b, T, ACC3, &iterations);
	method_iterations[0][2] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[0][2] = euclidean_norm(difference);

	printf("Вычисленное решение:\n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	//-----= Метод Якоби =-----
	printf("-=Метод Якоби=-\n\n");

	printf("Достаточное условие: ");
	
	bool sufficient_condition = true; //Достаточное условие
	for(uint8_t i = 0; i < matr->size; i++){
		float diagonal_elem;
		float sum = 0;

		for(uint8_t j = 0; j < matr->size; j++){
			if(i == j){
				diagonal_elem = matr->data[i][j];
				continue;
			}

			sum += matr->data[i][j];
		}

		if(diagonal_elem < sum){
			sufficient_condition = false;
			break;
		}
		sum = 0;
	}

	if(sufficient_condition){
		printf("да\n");
	}
	else{
		printf("нет\n");
	}

	printf("Сходимость: да\n\n");
	
	printf("Точность %g\n", ACC1);

	X = jacobi_method(matr, x0, b, ACC1, &iterations);
	method_iterations[1][0] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[1][0] = euclidean_norm(difference);

	printf("Вычисленное решение:\n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	printf("Точность %g\n", ACC2);

	X = jacobi_method(matr, x0, b, ACC2, &iterations);
	method_iterations[1][1] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[1][1] = euclidean_norm(difference);

	printf("Вычисленное решение:\n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");
	
	printf("Точность %g\n", ACC3);

	X = jacobi_method(matr, x0, b, ACC3, &iterations);
	method_iterations[1][2] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[1][2] = euclidean_norm(difference);

	printf("Вычисленное решение:\n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	//-----=Метод Гаусса-Зейделя=-----
	printf("-=Метод Гаусса-Зейделя=-\n\n");

	printf("Сходимость: да\n\n");

	printf("Точность: %g\n", ACC1);

	X = gauss_seidel_method(matr, x0, b, ACC1, &iterations);
	method_iterations[2][0] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[2][0] = euclidean_norm(difference);

	printf("Вычисленное решение: \n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);
	
	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	printf("Точность: %g\n", ACC2);

	X = gauss_seidel_method(matr, x0, b, ACC2, &iterations);
	method_iterations[2][1] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[2][1] = euclidean_norm(difference);

	printf("Вычисленное решение: \n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	printf("Точность: %g\n", ACC3);

	X = gauss_seidel_method(matr, x0, b, ACC3, &iterations);
	method_iterations[2][2] = iterations;

	difference = subtruct_vectors(exact_solution, X);
	difference_norms[2][2] = euclidean_norm(difference);

	printf("Вычисленное решение: \n");
	print_vector(X);
	free_vector(X);
	free_vector(difference);

	printf("Количество итераций: %u\n", iterations);
	iterations = 0;
	printf("\n");

	//Запись данных для графиков в файл
	FILE *file = fopen("data.txt", "w");
	
	for(uint8_t i = 0; i < 3; i++){
		for(uint8_t j = 0; j < 3; j++){
			fprintf(file, "%u ", method_iterations[i][j]);
		}
	}

	fprintf(file, "\n");

	for(uint8_t i = 0; i < 3; i++){
		for(uint8_t j = 0; j < 3; j++){
			fprintf(file, "%f ", difference_norms[i][j]);
		}
	}

	fclose(file);

	//Освобождение памяти
	free_matrix(matr);
	free_vector(b);
	free_vector(x0);
	free_vector(exact_solution);
	return 0;
}
