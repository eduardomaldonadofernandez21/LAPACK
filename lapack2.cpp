#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <mkl.h>
#include <random>

/*
1. Repetir la factorización LU y el cálculo de la matriz inversa realizado en la
práctica anterior utilizando la rutina LAPACKE_dgesv().
2. Realizar una comparación entre la operación sobre matrices generales y
sobre matrices banda.
a) Crear matrices A (tridiagonal) y B, y rellenarlas con valores aleatorios (media
0 y varianza 1)
b) Codificar la matriz A en forma compacta A_banda, añadiendo una fila auxiliar
nula al principio.
c) Resolver a partir de la matriz general (LAPACKE_dgesv, comprobar el código
de error).
d) Resolver a partir de la matriz banda (LAPACKE_dgbsv, comprobar el código de
error).
e) Comparar los tiempos c y d promediando entre diferentes ejecuciones.
3. [OPTATIVO] Realizar la comparativa de tiempos en función del ancho de
la banda.
*/

int showMat(int* mat, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			printf("%i ", mat[col * i + j]);
		}
		printf("\n");
	}
	printf("\n");
	return 0;
}

int showMat(double* mat, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			printf("%.4f ", mat[col * i + j]);
		}
		printf("\n");
	}
	printf("\n");
	return 0;
}

double* pClone(double* mat, double* c, int size) {
	for (int i = 0; i < size; i++) {
		c[i] = mat[i];
	}
	return c;
}

double* generateMatrix(int N) {
	double* matrix = (double*)mkl_malloc((N*N) * sizeof(double), 64);
	std::default_random_engine generador;
	std::normal_distribution<double> aleatorio(0.0, 1.0);

	for (int i = 0; i < N*N; i++) {
		matrix[i] = aleatorio(generador);
	}
	return matrix;
}

int* generatePivot(int N) {
	int* pivot = (int*)mkl_malloc(N * sizeof(int), 32);
	return pivot;
}

//Metodo que generar una matriz tridiagonal NxN
double* generateTridiagonal(int N) { 
	double* matrix = (double*)mkl_malloc((N*N) * sizeof(double), 64);
	std::default_random_engine generador;
	std::normal_distribution<double> aleatorio(0.0, 1.0);

	for (int i = 0; i < N*N; i++) {
		if (i % (N+1) == 0 || (i-1) % (N+1) == 0 || (i+1) % (N + 1) == 0) {
			matrix[i] = aleatorio(generador);
		}
		else {
			matrix[i] = 0;
		}
	}
	return matrix;
}

double* getMainDiagonal(double* m, int N) {
	double* diagonal = (double*)mkl_malloc(N * sizeof(double), 64);
	int size = 0;
	while (N > size) {
		for (int i = 0; i < N * N; i++) {
			if (i % (N + 1) == 0) {
				diagonal[size] = m[i];
				size = size + 1;
			}
		}
	}
	return diagonal;
}

double* getInfDiagonal(double* m, int N) {
	double* diagonal = (double*)mkl_malloc((N) * sizeof(double), 64);
	int size = 0;
	while (N > size) {
		for (int i = 0; i < N * N; i++) {
			if ((i + 1) % (N + 1) == 0) {
				diagonal[size] = m[i];
				size = size + 1;
			}
		}
	}
	diagonal[N - 1] = 0;
	return diagonal;
}

double* getSupDiagonal(double* m, int N) {
	double* diagonal = (double*)mkl_malloc((N) * sizeof(double), 64);
	diagonal[0] = 0;
	int size = 1;
	while (N > size) {
		for (int i = 0; i < N * N; i++) {
			if ((i - 1) % (N + 1) == 0) {
				diagonal[size] = m[i];
				size = size + 1;
			}
		}
	}
	return diagonal;
}

double* codificaMatrix(double* m, int N) {
	double* matrix = (double*)mkl_malloc((N * N) * sizeof(double), 64);
	double* mainDig = getMainDiagonal(m, N);
	double* supDig = getSupDiagonal(m, N);
	double* infDig = getInfDiagonal(m, N);
	printf("Diagonal principal\n");
	showMat(mainDig, 1, 4);
	printf("Diagonal superior\n");
	showMat(supDig, 1, 4);
	printf("Diagonal inferior\n");
	showMat(infDig, 1, 4);
	int j = 0;
	int cont = 0;
	while (N-1 > cont) {
		for (int i = 0; i < N * N; i++) { //Recorremos la nueva matriz codificada
			if (i < N) {
				matrix[i] = 0;
			}else {
				//Añadimos la diagonal superior
				if (i < (N*2)-1) {
					matrix[i] = supDig[j];
					j++;
				}
				if (i == (N * 2)-1) {
					matrix[i] = supDig[j];
					j = 0;
					cont++;
				}
				//Añadimos la diagonal principal
				if (i < ((N * N) - N-1) && i >= N*2) {
					matrix[i] = mainDig[j];
					j++;
				}
				if (i == ((N * N)-N-1)) {
					matrix[i] = mainDig[j];
					j = 0;
					cont++;
				}
				//Añadimos la diagonal inferior
				if (i > (N*N)-N-1) {
					matrix[i] = infDig[j];
					j++;
				}
				if ((N * N) - 1 == i) {
					matrix[i] = 0;
					cont++;
				}
			}
		}
		
	}
	return matrix;
}



int main(int argc, char* argv[]) {
	double A[36] = { 9.0,3.0,10.0,8.0,7.0,8.0  ,10.0,6.0,5.0,10.0,8.0,1.0,  2.0,10.0,9.0,7.0,8.0,3.0  ,10.0,10.0,2.0,1.0,4.0,1.0,  7.0,2.0,5.0,9.0,7.0,1.0,  1.0,10.0,10.0,10.0,2.0,9.0 };
	double B[36] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
	lapack_int* piv = new lapack_int[6];
	double c[36];
	memset(c, 0, 36 * sizeof(int));
	pClone(A, c, 36);

	LAPACKE_dgesv(LAPACK_ROW_MAJOR, 6, 6, A, 6, piv, B, 6);
	printf("Ejercicio 2:");
	printf("LU calculado con LAPACKE_dgesv driver\n");
	showMat(A, 6, 6);
	
	printf("Matriz inversa calculado con LAPACKE_dgesv driver\n");
	//int p[6] = { 0,0,0,0,0,0 };
	lapack_int* p = new lapack_int[6];
	double I[36] = { 1.0,0.0,0.0,0.0,0.0,0.0,  0.0,1.0,0.0,0.0,0.0,0.0, 0.0,0.0,1.0,0.0,0.0,0.0, 0.0,0.0,0.0,1.0,0.0,0.0, 0.0,0.0,0.0,0.0,1.0,0.0, 0.0,0.0,0.0,0.0,0.0,1.0 };
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, 6, 6, c, 6, p, I, 6);
	showMat(I, 6, 6);

	//Crear matrices A (tridiagonal) y B, y rellenarlas con valores aleatorios (media 0 y varianza 1)
	printf("a): \n");
	double* matA = generateTridiagonal(4);
	double* matB = generateMatrix(4);
	printf("Matrix A tridiagonal 4x4: \n");
	showMat(matA, 4,4);
	printf("Matrix B general 4x4: \n");
	showMat(matB, 4, 4);
	//Codificar la matriz A en forma compacta A_banda, añadiendo una fila auxiliar nula al principio
	printf("b): \n");
	double* dig = codificaMatrix(matA, 4);
	printf("Matriz A codificada: \n");
	showMat(dig, 4, 4);

	//Resolver a partir de la matriz general (LAPACKE_dgesv, comprobar el código de error).
	lapack_int* pivot = new lapack_int[6];
	double be[16] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, 6, 6, matB, 6, pivot, be, 6);

	mkl_free(dig);
	mkl_free(matA);
	mkl_free(matB);
	getchar();
	return 0;
	
}