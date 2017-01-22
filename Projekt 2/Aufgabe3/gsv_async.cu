/*
	Just run sh compileRun.sh
	Use config.h in order to adjust problem size
 */

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#include "config.h"





__device__
float func(float x, float y) {
	return 32 * (x * (1 - x) + y * (1 - y));
}

__global__
void initBase(float *base, float h) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		int x = i % D;
		int y = i / D;
		float f = func(h * x + h, h * y + h);
		base[i] = h * h * f;
	}
}

__global__
void calculate(float *uHistory, float *base, char *smallError, int sourceTime, int time, int lastTime, int k) {
	int i = 2 * (blockIdx.x * blockDim.x + threadIdx.x);
	if (i < N) {
		int x = i % D;
		int y = i / D;
		int diagIdx = (x + y) / 2;
		
		if (diagIdx < k) {
			float sum = base[i];
			if (y > 0) sum += uHistory[i - D + sourceTime];
			if (y < D - 1) sum += uHistory[i + D + sourceTime];
			if (x > 0) sum += uHistory[i - 1 + sourceTime];
			if (x < D - 1) sum += uHistory[i + 1 + sourceTime];
			sum /= 4;
			
			if (fabsf(sum - uHistory[i + lastTime]) >= EPSILON) {
				smallError[(k - diagIdx + D) % D] = 0;
			}
			
			uHistory[i + time] = sum;
		}
	}
	
	__syncthreads();
	
	i++;
	sourceTime = time;
	if (i < N) {
		int x = i % D;
		int y = i / D;
		int diagIdx = (x + y) / 2;
		
		if (diagIdx < k) {
			float sum = base[i];
			if (y > 0) sum += uHistory[i - D + sourceTime];
			if (y < D - 1) sum += uHistory[i + D + sourceTime];
			if (x > 0) sum += uHistory[i - 1 + sourceTime];
			if (x < D - 1) sum += uHistory[i + 1 + sourceTime];
			sum /= 4;
			
			if (fabsf(sum - uHistory[i + lastTime]) >= EPSILON) {
				smallError[(k - diagIdx + D) % D] = 0;
			}
			
			uHistory[i + time] = sum;
		}
	}
}

__global__
void fetchU(float *uHistory, float *u, int k) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		int x = i % D;
		int y = i / D;
		int diagIdx = (x + y) / 2;
		u[i] = uHistory[i + ((k + 1 + diagIdx) % D) * N];
	}
}



void solve(float h, float *u, int *iterations, int blockSize) {
	*iterations = 0;
	int halfN = (N + 1) / 2;
	int gridSizeN = (N + blockSize - 1) / blockSize;
	int gridSizeHalfN = (halfN + blockSize - 1) / blockSize;
	
	// Allocate memory
	float *base_d;
	cudaMalloc((void**) &base_d, N * sizeof(float));
	cudaMemset(base_d, 0, N * sizeof(float));
	initBase<<<gridSizeN, blockSize>>>(base_d, (float) h);
	
	float *uHistory_d;
	cudaMalloc((void**) &uHistory_d, D * N * sizeof(float));
	cudaMemset(uHistory_d, 0, D * N * sizeof(float));
	
	char *smallError_d;
	cudaMalloc((void**) &smallError_d, D);
	cudaMemset(smallError_d, 0, D);
	
	// Calculate u
	for (int k = 1; ; k++) {
		int time = (k % D) * N;
		int lastTime = ((k - 1 + D) % D) * N;
		
		cudaMemset(smallError_d + (k % D), 1, 1);
		
		// Black fields
		calculate<<<gridSizeHalfN, blockSize>>>(uHistory_d, base_d, smallError_d, lastTime, time, lastTime, k);
		
		// White fields
		//calculate<<<gridSizeHalfN, blockSize>>>(uHistory_d, base_d, smallError_d, time, time, lastTime, 1, k);
		
		(*iterations)++;
		
		int smallError;
		cudaMemcpy(&smallError, smallError_d + ((k + 1) % D), 1, cudaMemcpyDeviceToHost);
		if (smallError) break;
	}
	
	// Fetch result
	float* u_d;
	cudaMalloc((void**) &u_d, N * sizeof(float));
	fetchU<<<gridSizeN, blockSize>>>(uHistory_d, u_d, *iterations);
	cudaMemcpy(u, u_d, N * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(u_d);
	
	// Release memory
	cudaFree(base_d);
	cudaFree(uHistory_d);
	cudaFree(smallError_d);
}





float analyticU(float x, float y) {
	return 16 * x * (1 - x) * y * (1 - y);
}

int main() {
	int i, j;
	
	float u[N];
	float h = 1.f / (D + 1);
	
	cudaSetDevice(CUDA_DEVICE);
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(& prop, device);
	int blockSize = prop.warpSize;
	
	printf("Run on %s (device %d) with blocksize %d\n",
			prop.name, device, blockSize);
	
	printf("l = %d\nd = %d\nn = %d\n\n", L, D, N);
	
	int it;
	solve(h, u, &it, blockSize);
	
	if (SHOW_RESULTS) {
		printf("\nResult:\n");
		for (i = 0; i < D; i++) {
			for (j = 0; j < D; j++) {
				printf("%8.4f", u[j + D * i]);
			}
			printf("\n");
		}
		
		printf("\nAnalytic:\n");
		for (i = 0; i < D; i++) {
			for (j = 0; j < D; j++) {
				printf("%8.4f", analyticU(j * h + h, i * h + h));
			}
			printf("\n");
		}
		printf("\n");
	}
	
	float maxError = 0.f;
	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++) {
			float error = analyticU(j * h + h, i * h + h) - u[j + D * i];
			error = error > 0 ? error : -error;
			maxError = error > maxError ? error : maxError;
		}
	}
	printf("Max error: %4.8f\n", maxError);
	printf("Iterations: %d\n", it);
	
	return 0;
}
