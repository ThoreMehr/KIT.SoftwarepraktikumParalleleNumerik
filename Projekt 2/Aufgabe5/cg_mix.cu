/*
 * Just run sh compileRun.sh
 * Use config.h in order to adjust problem size
 */
 
#include <cuda.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"


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


__global__
void iterateILU(double* srcU, double* dstU, char* smallError) {
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
		
		if (fabs(value - srcU[5 * i + 2 * offset]) >= EPSILON_ILU) {
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

double func(double x, double y) {
	return 8 * M_PI * M_PI * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

__global__
void calculateGSV(double* A, double *uHistory, double *base, char *smallError, int sourceTime, int time, int lastTime, int offset, int k) {
	int i = 2 * (blockIdx.x * blockDim.x + threadIdx.x) + offset;
	if (i < N) {
		int x = i % D;
		int y = i / D;
		int diagIdx = (x + y) / 2;
		
		if (diagIdx < k) {
			double sum = base[i];
			if (y > 0) sum -= A[5 * i + 3] * uHistory[i - D + sourceTime];
			if (y < D - 1) sum -= A[5 * i + 4] * uHistory[i + D + sourceTime];
			if (x > 0) sum -= A[5 * i + 1] * uHistory[i - 1 + sourceTime];
			if (x < D - 1) sum -= A[5 * i + 2] * uHistory[i + 1 + sourceTime];
			sum /= A[5 * i];
			
			if (fabsf(sum - uHistory[i + lastTime]) >= EPSILON_GSV) {
				smallError[(k - diagIdx + D) % D] = 0;
			}
			
			uHistory[i + time] = sum;
		}
	}
}

__global__
void fetchU(double *uHistory, double *u, int k) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		int x = i % D;
		int y = i / D;
		int diagIdx = (x + y) / 2;
		u[i] = uHistory[i + ((k + 1 + diagIdx) % D) * N];
	}
}

void decompose(double* d_U_, double* d_L, char* d_smallError, int* iterations, int blockSize) {
	*iterations = 0;
	
	int gridSize3N = (3 * N + blockSize - 1) / blockSize;
	int gridSize5N = (5 * N + blockSize - 1) / blockSize;
	
	double* d_U[2];
	d_U[0] = d_U_;
	cudaMalloc((void**) &d_U[1], 5 * N * sizeof(double));
	
	cudaMemset(d_smallError, 0, 1);
	
	// Initialize matrices with identity
	initMatrix<<<gridSize5N, blockSize>>>(d_U[0]);
	initMatrix<<<gridSize5N, blockSize>>>(d_U[1]);
	
	for (int m = 0;; m++) {
		cudaMemset(d_smallError, 1, 1);
		
		iterateILU<<<gridSize3N,blockSize>>>(d_U[m % 2 == 0], d_U[m % 2], d_smallError);
		
		(*iterations)++;
		
		char smallError;
		cudaMemcpy(&smallError, d_smallError, 1, cudaMemcpyDeviceToHost);
		if (smallError && *iterations % 2 == 0) break;
	}
	
	transpose<<<gridSize5N, blockSize>>>(d_L, d_U[(*iterations) % 2 == 0]);
	
	cudaFree(d_U[1]);
}

void solveGSV(double* d_A, double* d_u, double* d_b, double* d_uHistory, char* d_smallError, int *iterations, int blockSize) {
	*iterations = 0;
	int halfN = (N + 1) / 2;
	int gridSizeN = (N + blockSize - 1) / blockSize;
	int gridSizeHalfN = (halfN + blockSize - 1) / blockSize;
	
	cudaMemset(d_smallError, 0, D);
	
	// Calculate u
	for (int k = 1; ; k++) {
		int time = (k % D) * N;
		int lastTime = ((k - 1 + D) % D) * N;
		
		cudaMemset(d_smallError + (k % D), 1, 1);
		
		// Black fields
		calculateGSV<<<gridSizeHalfN, blockSize>>>(d_A, d_uHistory, d_b, d_smallError, lastTime, time, lastTime, 0, k);
		
		// White fields
		calculateGSV<<<gridSizeHalfN, blockSize>>>(d_A, d_uHistory, d_b, d_smallError, time, time, lastTime, 1, k);
		
		(*iterations)++;
		
		if (k >= D) {
			char smallError;
			cudaMemcpy(&smallError, d_smallError + ((k + 1) % D), 1, cudaMemcpyDeviceToHost);
			if (smallError) break;
		}
	}
	
	// Fetch result
	fetchU<<<gridSizeN, blockSize>>>(d_uHistory, d_u, *iterations);
}

void solveBr(double* d_L, double* d_U, double* d_r, double* d_p, double* d_tmp, double* d_uHistory, char* d_smallError, int blockSize) {
	int it;
	
	solveGSV(d_L, d_tmp, d_r, d_uHistory, d_smallError, &it, blockSize);
	solveGSV(d_U, d_p, d_tmp, d_uHistory, d_smallError, &it, blockSize);
}

void solve(double *x, int *iterations, int blockSize) {
	*iterations = 0;
	
	// Allocate memory
	
	double *d_uHistory;
	cudaMalloc((void**) &d_uHistory, D * N * sizeof(double));
	cudaMemset(d_uHistory, 0, D * N * sizeof(double));
	
	char *d_smallError;
	cudaMalloc((void**) &d_smallError, D);
	cudaMemset(d_smallError, 0, D);
	
	double* d_U;
	cudaMalloc((void**) &d_U, 5 * N * sizeof(double));
	
	double* d_L;
	cudaMalloc((void**) &d_L, 5 * N * sizeof(double));
	
	int it;
	decompose(d_U, d_L, d_smallError, &it, blockSize);
	printf("%d iterations for ILU decomposition\n", it);
	
	double* d_r;
	cudaMalloc((void**) &d_r, N * sizeof(double));
	
	double* d_p;
	cudaMalloc((void**) &d_p, N * sizeof(double));
	
	double* d_tmp0;
	cudaMalloc((void**) &d_tmp0, N * sizeof(double));
	double* d_tmp1;
	cudaMalloc((void**) &d_tmp1, N * sizeof(double));
	
	double delta = 0, deltaHat;
	
	double* Br = (double*) malloc(N * sizeof(double));
	double* p = (double*) malloc(N * sizeof(double));
	double* r = (double*) malloc(N * sizeof(double));
	double* base = (double*) malloc(N * sizeof(double));
	int i, ix;
	
	#pragma omp parallel for private(i)
	for (i = 0; i < N; i++) {
		x[i] = 1;
		int x = i % D;
		int y = i / D;
		double f = func(H * x + H, H * y + H);
		base[i] = H * H * f;
	}
	
	//initR0<<<gridSizeN,blockSize>>>(d_r, d_base, d_x);
	#pragma omp parallel for private(i,ix)
	for (i = 0; i < N; i++) {
		int ix = i % D;
		double residuum = base[i];
		if (ix - 1 >= 0) residuum += x[i - 1];
		if (ix + 1 < D) residuum += x[i + 1];
		if (i - D >= 0) residuum += x[i - D]; 
		if (i + D < N) residuum += x[i + D];
		residuum -= 4 * x[i];
		r[i] = residuum;
	}
	
	cudaMemcpy(d_r, r, N * sizeof(double), cudaMemcpyHostToDevice);
	solveBr(d_L, d_U, d_r, d_p, d_tmp0, d_uHistory, d_smallError, blockSize);
	cudaMemcpy(p, d_p, N * sizeof(double), cudaMemcpyDeviceToHost);
	
	//cudaMemset(d_delta, 0, 2 * sizeof(double));
	//scalarProduct<<<gridSizeN,blockSize>>>(d_r, d_p, d_delta);
	
	
	delta = 0;
	#pragma omp parallel for private(i) reduction(+:delta)
	for (i = 0; i < N; i++) {
		delta += r[i] * p[i];
	}
	
	
	while (delta >= EPSILON * EPSILON) {
		
		/*calculateAx<<<gridSizeN,blockSize>>>(d_p, d_tmp0);
		scalarProduct<<<gridSizeN,blockSize>>>(d_p, d_tmp0, d_delta + 1);
		cudaMemcpy(&deltaHat, d_delta + 1, sizeof(double), cudaMemcpyDeviceToHost);
		deltaHat = delta / deltaHat;*/
		
		deltaHat = 0;
		
		#pragma omp parallel for private(i, ix) reduction(+:deltaHat)
		for (i = 0; i < N; i++) {
			ix = i % D;
			double v = 0;
			if (ix - 1 >= 0) v -= p[i - 1];
			if (ix + 1 < D) v -= p[i + 1];
			if (i - D >= 0) v -= p[i - D]; 
			if (i + D < N) v -= p[i + D];
			v += 4 * p[i];
			deltaHat += p[i] * v;
		}
		deltaHat = delta / deltaHat;
		
		
		/*addfv<<<gridSizeN,blockSize>>>(d_x, d_x, deltaHat, d_p);
		addfv<<<gridSizeN,blockSize>>>(d_r, d_r, -deltaHat, d_tmp0);*/
		
		#pragma omp parallel for private(i, ix)
		for (i = 0; i < N; i++) {
			ix = i % D;
			double v = 0;
			if (ix - 1 >= 0) v -= p[i - 1];
			if (ix + 1 < D) v -= p[i + 1];
			if (i - D >= 0) v -= p[i - D]; 
			if (i + D < N) v -= p[i + D];
			v += 4 * p[i];
			x[i] += deltaHat * p[i];
			r[i] -= deltaHat * v;
		}
		
		cudaMemcpy(d_r, r, N * sizeof(double), cudaMemcpyHostToDevice);
		//cudaMemset(d_delta, 0, 2 * sizeof(double));
		solveBr(d_L, d_U, d_r, d_tmp1, d_tmp0, d_uHistory, d_smallError, blockSize);
		//scalarProduct<<<gridSizeN,blockSize>>>(d_r, d_tmp1, d_delta);
		//cudaMemcpy(&newDelta, d_delta, sizeof(double), cudaMemcpyDeviceToHost);*/
		cudaMemcpy(Br, d_tmp1, N * sizeof(double), cudaMemcpyDeviceToHost);
		
		double newDelta = 0;
		#pragma omp parallel for private(i) reduction(+:newDelta)
		for (i = 0; i < N; i++) {
			newDelta += r[i] * Br[i];
		}
		
		//addfv<<<gridSizeN,blockSize>>>(d_p, d_tmp1, newDelta / delta, d_p);
		
		delta = newDelta / delta;
		#pragma omp parallel for private(i)
		for (i = 0; i < N; i++) {
			p[i] = Br[i] + delta * p[i];
		}
		
		delta = newDelta;
		(*iterations)++;
	}
	
	// Release memory
	cudaFree(d_uHistory);
	cudaFree(d_smallError);
	cudaFree(d_U);
	cudaFree(d_L);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_tmp0);
	cudaFree(d_tmp1);
	
	free(Br);
	free(p);
	free(r);
	free(base);
}



float analyticU(double x, double y) {
	return sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

int main(void) {
	int i, j;
	
	double u[N];
	
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
	solve(u, &it, blockSize);
	
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
				printf("%8.4f", analyticU(j * H + H, i * H + H));
			}
			printf("\n");
		}
		printf("\n");
	}
	
	double maxError = 0.0;
	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++) {
			double error = fabs(analyticU(j * H + H, i * H + H) - u[j + D * i]);
			maxError = fmax(error, maxError);
		}
	}
	printf("Max error: %4.8f\n", maxError);
	printf("Iterations: %d\n", it);
	
	return 0;
}
