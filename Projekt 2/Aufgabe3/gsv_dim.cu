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
void calculate(float *uHistory, float *base, char *smallError, int sourceTime, int time, int lastTime, int offset, int k) {
	extern __shared__ float values[];
	
	int lx = threadIdx.x;
	int ly = threadIdx.y;
	int x = 2 * (blockIdx.x * blockDim.x + lx);
	int y = blockIdx.y * blockDim.y + ly;
	x += (y + offset) % 2;
	int i = x + y * D;
	
	int width = 2 * blockDim.x + 2;
	int height = blockDim.y + 2;
	int sx = 2 * blockIdx.x * blockDim.x - 1;
	int sy = blockIdx.y * blockDim.y - 1;
	int myIndex = lx + ly * blockDim.x;

	
	int idx = 4 * myIndex + 1 - offset;
	if (idx < width * height) {
		int cx = sx + idx % width;
		int cy = sy + idx / width;
		if (cx <= D && cy <= D && idx < width * height) {
			if (cx >= 0 && cy >= 0 && cx < D && cy < D) {
				values[idx] = uHistory[cx + cy * D + sourceTime];
			} else {
				values[idx] = 0;
			}
		}
		
		idx += 2;
		if (idx < width * height) {
			cx = sx + idx % width;
			cy = sy + idx / width;
			if (cx <= D && cy <= D) {
				if (cx >= 0 && cy >= 0 && cx < D && cy < D) {
					values[idx] = uHistory[cx + cy * D + sourceTime];
				} else {
					values[idx] = 0;
				}
			}
		}
	}
	
	
	/*if (myIndex == 0) {
		for (int idx = 1 - offset; idx < 4 * blockDim.x * blockDim.y; idx += 2) {
			if (idx < width * height) {
				int cx = sx + idx % width;
				int cy = sy + idx / width;
				if (cx <= D && cy <= D) {
					if (cx >= 0 && cy >= 0 && cx < D && cy < D) {
						values[idx] = uHistory[cx + cy * D + sourceTime];
					} else {
						values[idx] = 0;
					}
				}
			}
		}
	}*/
	
	/*if (myIndex == 0) {
		for (int idx = 1 - offset; idx < width * height; idx += 2) {
			int cx = sx + idx % width;
			int cy = sy + idx / width;
			if (cx >= 0 && cy >= 0 && cx < D && cy < D) {
				values[idx] = uHistory[cx + cy * D + sourceTime];
			} else {
				values[idx] = 0;
			}
		}
	}*/
	
	__syncthreads();
	
	if (x < D && y < D) {
		int diagIdx = (x + y) / 2;
		int s = x - sx + (y - sy) * width;
		
		if (diagIdx < k) {
			float sum = base[i];
			sum += values[s - width];
			sum += values[s + width];
			sum += values[s - 1];
			sum += values[s + 1];
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



void solve(float h, float *u, int *iterations, int blockSize, int blockSizeX, int blockSizeY) {
	*iterations = 0;
	
	int gridSizeN = (N + blockSize - 1) / blockSize;
	int halfD = (D + 1) / 2;
	dim3 dimBlock(blockSizeX, blockSizeY, 1);
	dim3 dimGrid((halfD + blockSizeX - 1) / blockSizeX, (D + blockSizeY - 1) / blockSizeY, 1);
	int parcelSize = sizeof(float) * (2 * blockSizeX + 2) * (blockSizeY + 2);
	
	// Allocate memory
	float *base_d;
	cudaMalloc((void**) &base_d, N * sizeof(float));
	cudaMemset(base_d, 0, N * sizeof(float));
	initBase<<<gridSizeN, blockSize>>>(base_d, (float) h);
	
	float *uHistory_d;
	cudaMalloc((void**) &uHistory_d, D * N * sizeof(float));
	cudaMemset(uHistory_d, 0, D * N * sizeof(float));
	
	char *smallError_d;
	cudaMalloc((void**) &smallError_d, D * sizeof(char));
	cudaMemset(smallError_d, 0, D);
	
	// Calculate u
	for (int k = 1; ; k++) {
		int time = (k % D) * N;
		int lastTime = ((k - 1 + D) % D) * N;
		
		cudaMemset(smallError_d + (k % D), 1, 1);
		
		// Black fields
		calculate<<<dimGrid, dimBlock, parcelSize>>>(uHistory_d, base_d, smallError_d, lastTime, time, lastTime, 0, k);
		
		// White fields
		calculate<<<dimGrid, dimBlock, parcelSize>>>(uHistory_d, base_d, smallError_d, time, time, lastTime, 1, k);
		
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
	int blockSizeX = 4;
	int blockSizeY = 8;
	
	printf("Run on %s (device %d) with blocksize %d and %dx%d\n",
			prop.name, device, blockSize, blockSizeX, blockSizeY);
	
	printf("l = %d\nd = %d\nn = %d\n\n", L, D, N);
	
	int it;
	solve(h, u, &it, blockSize, blockSizeX, blockSizeY);
	
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