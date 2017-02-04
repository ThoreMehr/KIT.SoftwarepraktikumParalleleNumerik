/*
 * Just run sh compileRun.sh
 * Use config.h in order to adjust problem size
 */
 
#include <cuda.h>
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

__device__
double func(double x, double y) {
	return 8 * M_PI * M_PI * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

__global__
void initBase(double *base, double h) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		int x = i % D;
		int y = i / D;
		double f = func(h * x + h, h * y + h);
		base[i] = h * h * f;
	}
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

__global__
void initX(double* x, double value) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		x[i] = value;
	}
}

__global__
void initR0(double* d_r, double* d_b, double* d_x) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		int ix = i % D;
		double residuum = d_b[i];
		
		if (ix - 1 >= 0) residuum += d_x[i - 1];
		if (ix + 1 < D) residuum += d_x[i + 1];
		if (i - D >= 0) residuum += d_x[i - D]; 
		if (i + D < N) residuum += d_x[i + D];
		residuum -= 4 * d_x[i];
		
		d_r[i] = residuum;
	}
}

__global__
void calculateAx(double* d_x, double* d_r) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		int ix = i % D;
		double residuum = 0;
		
		if (ix - 1 >= 0) residuum -= d_x[i - 1];
		if (ix + 1 < D) residuum -= d_x[i + 1];
		if (i - D >= 0) residuum -= d_x[i - D]; 
		if (i + D < N) residuum -= d_x[i + D];
		residuum += 4 * d_x[i];
		
		d_r[i] = residuum;
	}
}

#if __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

__global__
void scalarProduct(double* a, double* b, double* sum) {
	__shared__ double localSum[1];
	if (threadIdx.x == 0) {
		localSum[0] = 0;
	}
	
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		atomicAdd(localSum, a[i] * b[i]);
	}
	
	if (threadIdx.x == 0) {
		atomicAdd(sum, localSum[0]);
	}
}

__global__
void addfv(double* a, double* b, double f, double* c) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N) {
		a[i] = b[i] +  f * c[i];
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
	//printf("%d Ly=r iterations\n", it);
	solveGSV(d_U, d_p, d_tmp, d_uHistory, d_smallError, &it, blockSize);
	//printf("%d Up=y iterations\n", it);
}

__global__
void dotProduct(double* t, double* a, double* b) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < N) {
		t[i] = a[i] * b[i];
	}
}

__global__
void sumReduction(double* a, int n) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < n / 2) {
		a[i] += a[(n + 1) / 2 + i];
	}
}

double reductionScalarProduct(double* d_a, double* d_b, double* d_tmp, int blockSize) {
	int size = N;
	int gridSizeN = (N + blockSize - 1) / blockSize;
	dotProduct<<<gridSizeN,blockSize>>>(d_tmp, d_a, d_b);
	
	while (size > 1) {
		int gridSize = ((size + 1) / 2 + blockSize - 1) / blockSize;
		sumReduction<<<gridSize,blockSize>>>(d_tmp, size);
		size = (size + 1) / 2;
	}
	
	double result;
	cudaMemcpy(&result, d_tmp, sizeof(double), cudaMemcpyDeviceToHost);
	return result;
}
	

void solve(double *u, int *iterations, int blockSize) {
	*iterations = 0;
	int gridSizeN = (N + blockSize - 1) / blockSize;
	
	// Allocate memory
	double *d_base;
	cudaMalloc((void**) &d_base, N * sizeof(double));
	initBase<<<gridSizeN, blockSize>>>(d_base, H);
	
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
	
	double* d_x;
	cudaMalloc((void**) &d_x, N * sizeof(double));
	initX<<<gridSizeN,blockSize>>>(d_x, 1);
	
	double* d_r;
	cudaMalloc((void**) &d_r, N * sizeof(double));
	
	double* d_p;
	cudaMalloc((void**) &d_p, N * sizeof(double));
	
	double* d_tmp0;
	cudaMalloc((void**) &d_tmp0, N * sizeof(double));
	double* d_tmp1;
	cudaMalloc((void**) &d_tmp1, N * sizeof(double));
	
	double* d_delta;
	cudaMalloc((void**) &d_delta, 2 * sizeof(double));
	
	double delta = 0, newDelta, deltaHat;
	
	initR0<<<gridSizeN,blockSize>>>(d_r, d_base, d_x);
	
	solveBr(d_L, d_U, d_r, d_p, d_tmp0, d_uHistory, d_smallError, blockSize);
	
	//cudaMemset(d_delta, 0, 2 * sizeof(double));
	//scalarProduct<<<gridSizeN,blockSize>>>(d_r, d_p, d_delta);
	//cudaMemcpy(&delta, d_delta, sizeof(double), cudaMemcpyDeviceToHost);
	delta = reductionScalarProduct(d_r, d_p, d_tmp0, blockSize);
	
	while (delta >= EPSILON * EPSILON) {
		//cudaMemset(d_delta, 0, 2 * sizeof(double));
		
		calculateAx<<<gridSizeN,blockSize>>>(d_p, d_tmp0);
		
		deltaHat = reductionScalarProduct(d_p, d_tmp0, d_tmp1, blockSize);
		//scalarProduct<<<gridSizeN,blockSize>>>(d_p, d_tmp0, d_delta + 1);
		//cudaMemcpy(&deltaHat, d_delta + 1, sizeof(double), cudaMemcpyDeviceToHost);
		deltaHat = delta / deltaHat;
		
		addfv<<<gridSizeN,blockSize>>>(d_x, d_x, deltaHat, d_p);
		addfv<<<gridSizeN,blockSize>>>(d_r, d_r, -deltaHat, d_tmp0);
		
		solveBr(d_L, d_U, d_r, d_tmp1, d_tmp0, d_uHistory, d_smallError, blockSize);
		
		newDelta = reductionScalarProduct(d_r, d_tmp1, d_tmp0, blockSize);
		//scalarProduct<<<gridSizeN,blockSize>>>(d_r, d_tmp1, d_delta);
		//cudaMemcpy(&newDelta, d_delta, sizeof(double), cudaMemcpyDeviceToHost);
		
		addfv<<<gridSizeN,blockSize>>>(d_p, d_tmp1, newDelta / delta, d_p);
		
		delta = newDelta;
		(*iterations)++;
	}
	
	cudaMemcpy(u, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
	
	// Release memory
	cudaFree(d_base);
	cudaFree(d_uHistory);
	cudaFree(d_smallError);
	cudaFree(d_U);
	cudaFree(d_L);
	cudaFree(d_x);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_tmp0);
	cudaFree(d_tmp1);
	cudaFree(d_delta);
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
