/*
 * Compile with
 * gcc -Wall -g -fopenmp -o $(NAME_PART) $(FILE_NAME) -lm
 */
 
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#include "config.h"

double func(double x, double y) {
	return 32 * (x * (1 - x) + y * (1 - y));
}

double analyticU(double x, double y) {
	return 16 * x * (1 - x) * y * (1 - y);
}

/*
 * Solves Au = h²f for with A of size n*n, u, f of size n
 */
void solve(double h, double *f, int n, double *u, int *iterations) {
	int i, j, k;
	double sum;
	double base[N];
	double *uHistory = (double*) malloc(D * N * sizeof(double));
	*iterations = 0;
	int smallError[D];
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < D; j++) {
			uHistory[i + j * N] = 0;
		}
		base[i] = h * h * f[i];
	}
	
	for (i = 0; i < D; i++) {
		smallError[i] = 0;
	}
	
	for (k = 1; ; k++) {
		int time = (k % D) * N;
		int lastTime = ((k - 1 + D) % D) * N;
		smallError[k % D] = 1;
		
		// Black fields
		#pragma omp parallel for private(j, sum)
		for (j = 0; j < N; j += 2) {
			int x = j % D;
			int y = j / D;
			int diagIndex = (x + y) / 2;
			
			if (diagIndex < k) {
				sum = base[j];
				if (j - D >= 0) sum += uHistory[j - D + lastTime];
				if (j + D < n) sum += uHistory[j + D + lastTime];
				if (j - 1 >= 0 && (j - 1) / D == j / D) sum += uHistory[j - 1 + lastTime];
				if (j + 1 < n && (j + 1) / D == j / D) sum += uHistory[j + 1 + lastTime];
				sum /= 4;
				
				if (fabs(sum - uHistory[j + lastTime]) >= EPSILON) smallError[(k - diagIndex + D) % D] = 0;
				
				uHistory[j + time] = sum;
			}
		}
		
		// White fields
		#pragma omp parallel for private(j, sum)
		for (j = 1; j < N; j += 2) {
			int x = j % D;
			int y = j / D;
			int diagIndex = (x + y) / 2;
			
			if (diagIndex < k) {
				sum = base[j];
				if (j - D >= 0) sum += uHistory[j - D + time];
				if (j + D < n) sum += uHistory[j + D + time];
				if (j - 1 >= 0 && (j - 1) / D == j / D) sum += uHistory[j - 1 + time];
				if (j + 1 < n && (j + 1) / D == j / D) sum += uHistory[j + 1 + time];
				sum /= 4;
				
				if (fabs(sum - uHistory[j + lastTime]) >= EPSILON) smallError[(k - diagIndex + D) % D] = 0;
				
				uHistory[j + time] = sum;
			}
		}
		
		(*iterations)++;
		
		if (smallError[(k + 1) % D]) {
			break;
		}
	}
	
	for (j = 0; j < N; j++) {
		int x = j % D;
		int y = j / D;
		int diagIndex = (x + y) / 2;
		u[j] = uHistory[j + ((k + 1 + diagIndex) % D) * N];
	}
	
	free(uHistory);
}

int main(void) {
	int i, j;
	
	double f[N];
	double u[N];
	double h = 1.0 / (D + 1);
	
	printf("Run on system with %d processors and max %d threads\n", omp_get_num_procs(), omp_get_max_threads());
	
	printf("l = %d\nd = %d\nn = %d\n\n", L, D, N);
	
	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++) {
			f[i + D * j] = func(i * h + h, j * h + h);
		}
	}
	
	int it;
	solve(h, f, N, u, &it);
	
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
	
	double maxError = 0.0;
	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++) {
			double error = analyticU(j * h + h, i * h + h) - u[j + D * i];
			error = error > 0 ? error : -error;
			maxError = error > maxError ? error : maxError;
		}
	}
	printf("Max error: %4.8f\n", maxError);
	printf("Iterations: %d\n", it);
	
	return 0;
}
