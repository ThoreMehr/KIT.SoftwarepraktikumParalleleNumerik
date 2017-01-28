/*
 * Sparse matrix format is array of 5 * N (5 entries per node of u)
 * Order: node itself, left neighbor, right neighbor, top neighbor, bottom neighbor
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"

void decompose(double* mat_L, double* mat_U, int* iterations) {
	*iterations = 0;
	
	// Initialize arrays with zero
	for (int i = 0; i < 5 * N; i++) {
		mat_L[i] = mat_U[i] = 0;
	}
	
	for (int j = 0; j < N; j++) {
		double sum = 0;
		// Access mat_L[j,k] for k < j
		int x = j % D;
		int y = j / D;
		if (x > 0) sum += mat_L[5 * j + 1] * mat_L[5 * j + 1];
		if (y > 0) sum += mat_L[5 * j + 3] * mat_L[5 * j + 3];
		mat_U[5 * j] = mat_L[5 * j] = sqrt(4 - sum);
		
		// Look at i > j for which a_ij != 0
		if (x < D - 1) {
			int i = j + 1;
			mat_U[5 * j + 2] = mat_L[5 * i + 1] = -1 / mat_L[5 * j];
		}
		if (y < D - 1) {
			int i = j + D;
			mat_U[5 * j + 4] = mat_L[5 * i + 3] = -1 / mat_L[5 * j];
		}
		
		(*iterations)++;
	}
}

void printSparseMatrix(double* A) {
	for (int i = 0; i < N; i++) {
		int x = i % D;
		int y = i / D;
		for (int j = 0; j < N; j++) {
			double value = 0;
			if (j == i) {
				value = A[5 * i];
			} else if (j == i - 1 && x > 0) {
				value = A[5 * i + 1];
			} else if (j == i + 1 && x < D - 1) {
				value = A[5 * i + 2];
			} else if (j == i - D && y > 0) {
				value = A[5 * i + 3];
			} else if (j == i + D && y < D - 1) {
				value = A[5 * i + 4];
			}
			printf("%10.6f", value);
		}
		printf("\n");
	}
}

int main() {
	double* mat_L = (double*) malloc(5 * N * sizeof(double));
	double* mat_U = (double*) malloc(5 * N * sizeof(double));
	int iterations;
	
	decompose(mat_L, mat_U, &iterations);
	
	printf("Used %d iterations\n\n", iterations);
	
	if (SHOW_RESULTS) {
		printf("Matrix L:\n");
		printSparseMatrix(mat_L);
		
		printf("\nMatrix U:\n");
		printSparseMatrix(mat_U);
	}
	
	free(mat_L);
	free(mat_U);
	
	return 0;
}

