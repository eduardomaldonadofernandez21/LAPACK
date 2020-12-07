#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <mkl.h>

int showMat(double* mat, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			printf("%.4f ", mat[col*i+j]);
		}
		printf("\n");
	}
	printf("\n");
	return 0;
}

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

int getDet(double* mat, int row, double det) {
	for (int i = 0; i < row; i++) {
		det *= mat[row * i + i];
	}
	return det;
}

double* pClone(double* mat, double* c, int size) {
	for (int i = 0; i < size; i++) {
		c[i] = mat[i];
	}
	return c;
}

int main() {
	double arr[36] = { 9,3,10,8,7,8,10,6,5,10,8,1,2,10,9,7,8,3,10,10,2,1,4,1,7,2,5,9,7,1,1,10,10,10,2,9 };
	double c[36];
	memset(c, 0, 36 * sizeof(int));
	pClone(arr, c, 36);

	int p[6] = { 0,0,0,0,0,0 };
	// a) Factorización LU
	showMat(arr, 6, 6);
	showMat(p, 1, 6);
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 6, 6, arr, 6, p);
	printf("LU\n");
	showMat(arr, 6, 6);
	showMat(p, 1, 6);
	// b) Determinante
	double det = 1;
	printf("%f\n", det);
	det = getDet(arr, 6, det);
	printf("%f\n", det);
	pClone(arr, c, 36);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, 6, c, 6, p);
	//d) Calcular la inversa usando la rutina _dgetri()
	printf("INV\n");
	showMat(c, 6, 6);
	// c) Matriz inversa a partir de resolver el sistema AX = I
	double I[36] = { 1,0,0,0,0,  0,1,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1 };
	double x[36] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',6, 6, arr, 6, p, I, 6);
	printf("INV a partir del sistema A*X = I\n");
	showMat(I, 6, 6);
	getchar();
	return 0;
}