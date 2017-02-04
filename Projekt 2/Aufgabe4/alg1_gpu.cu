/*
 * Sparse matrix format is array of 5 * N (5 entries per node of u)
 * Order: node itself, left neighbor, right neighbor, top neighbor, bottom neighbor
 */

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"


__global__
void iterate(double* srcU, double* dstU, char* smallError) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int offset = index % 3;
	int i = index / 3;
	int x = i % D;
	int y = i / D;
	if (i < N && (offset == 0 || x < D - 1 && offset == 1 || y < D - 1 && offset == 2)) {
		double value = 0;
		
		if (offset == 0) {
			if (x > 0) value += srcU[5 * (i - 1) + 2] * srcU[5 * (i - 1) + 2];
			if (y > 0) value += srcU[5 * (i - D) + 4] * srcU[5 * (i - D) + 4];
			value = sqrt(4 - value);
		} else {
			value = -1 / srcU[5 * i];
		}
		
		dstU[5 * i + 2 * offset] = value;
		
		if (fabs(value - srcU[5 * i + 2 * offset]) >= EPSILON) {
			*smallError = 0;
		}
	}
}

__global__
void initMatrix(double* A) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 5 * N) {
		A[i] = (i % 5 == 0);
	}
}

__global__
void transpose(double* A, double* B) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < 5 * N) {
		int offset = i % 5;
		if (offset == 0) {
			A[i] = B[i];
		} else {
			int j = i / 5 - D * (offset == 3) + D * (offset == 4)
				  - (offset == 1) + (offset == 2);
			A[i] = B[5 * j + 2 - (offset - 1) % 2 + 2 * ((offset - 1) / 2)];
		}
	}
}

void decompose(double* mat_L, double* mat_U, int* iterations, int blockSize) {
	*iterations = 0;
	
	int gridSize3N = (3 * N + blockSize - 1) / blockSize;
	int gridSize5N = (5 * N + blockSize - 1) / blockSize;
	
	double* d_U[2];
	cudaMalloc((void**) &d_U[0], 5 * N * sizeof(double));
	cudaMalloc((void**) &d_U[1], 5 * N * sizeof(double));
	
	double* d_L;
	cudaMalloc((void**) &d_L, 5 * N * sizeof(double));
	
	char* d_smallError;
	cudaMalloc((void**) &d_smallError, 1);
	cudaMemset(d_smallError, 0, 1);
	
	// Initialize matrices with identity
	initMatrix<<<gridSize5N, blockSize>>>(d_U[0]);
	initMatrix<<<gridSize5N, blockSize>>>(d_U[1]);
	
	for (int m = 0;; m++) {
		cudaMemset(d_smallError, 1, 1);
		
		iterate<<<gridSize3N,blockSize>>>(d_U[m % 2 == 0], d_U[m % 2 == 1], d_smallError);
		
		(*iterations)++;
		
		char smallError;
		cudaMemcpy(&smallError, d_smallError, 1, cudaMemcpyDeviceToHost);
		if (smallError) break;
	}
	
	transpose<<<gridSize5N, blockSize>>>(d_L, d_U[(*iterations) % 2 == 0]);
	
	cudaMemcpy(mat_U, d_U[(*iterations) % 2 == 0], 5 * N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(mat_L, d_L, 5 * N * sizeof(double), cudaMemcpyDeviceToHost);
	
	cudaFree(d_U[0]);
	cudaFree(d_U[1]);
	cudaFree(d_L);
	cudaFree(d_smallError);
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
	
	cudaSetDevice(CUDA_DEVICE);
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(& prop, device);
	int blockSize = prop.warpSize;
	
	printf("Run on %s (device %d) with blocksize %d\n",
			prop.name, device, blockSize);
	
	printf("l = %d\nd = %d\nn = %d\n\n", L, D, N);

	
	decompose(mat_L, mat_U, &iterations, blockSize);
	
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

