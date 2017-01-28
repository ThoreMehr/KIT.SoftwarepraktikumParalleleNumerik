/*
 * Sparse matrix format is array of 5 * N (5 entries per node of u)
 * Order: node itself, left neighbor, right neighbor, top neighbor, bottom neighbor
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"

void swap(double** a, double** b) {
	double* c = *a;
	*a = *b;
	*b = c;
}

double max(double a, double b) {
	return a > b ? a : b;
}

void decompose(double* mat_L, double* mat_U, int* iterations) {
	*iterations = 0;
	
	double* mat_altL = (double*) malloc(5 * N * sizeof(double));
	double* mat_altU = (double*) malloc(5 * N * sizeof(double));
	
	// Initialize arrays with zero
	for (int i = 0; i < 5 * N; i++) {
		mat_altL[i] = mat_altU[i] = mat_L[i] = mat_U[i] = i % 5 == 0;
	}
	
	for (int m = 0;; m++) {
		double delta = 0;
		
		for (int i = 0; i < N; i++) {
			int x = i % D;
			int y = i / D;
			int js[3] = { i };
			int jCount = 1;
			if (x < D - 1) js[jCount++] = i + 1;
			if (y < D - 1) js[jCount++] = i + D;
			
			for (int idx = 0; idx < jCount; idx++) {
				int j = js[idx];
				if (i == j) {
					double sum = 0;
					if (x > 0) sum += mat_U[5 * (i - 1) + 2] * mat_U[5 * (i - 1) + 2];
					if (y > 0) sum += mat_U[5 * (i - D) + 4] * mat_U[5 * (i - D) + 4];
					double sm = 4 - sum;
					mat_altU[5 * i] = mat_altL[5 * i] = sqrt(sm);
					delta = max(fabsf(mat_altU[5 * i] - mat_U[5 * i]), delta);
				} else {
					double sm = -1;
					mat_altL[5 * j + 1 + 2 * (j - i == D)] = mat_altU[5 * i + 2 + 2 * ((j - i) == D)] = sm / mat_U[5 * i];
					delta = max(fabsf(mat_altU[5 * i + 2 + 2 * ((j - i) == D)] - mat_U[5 * i + 2 + 2 * ((j - i) == D)]), delta);
				}
			}
			
		}
		
		swap(&mat_altL, &mat_L);
		swap(&mat_altU, &mat_U);
		
		(*iterations)++;
		
		if (delta < EPSILON && *iterations % 2 == 0) {
			break;
		}
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
	
	printf("Calculate LU decomposition of A for N=%d\n", N);
	
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

