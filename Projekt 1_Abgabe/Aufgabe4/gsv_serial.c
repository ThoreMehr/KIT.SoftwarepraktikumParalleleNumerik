/*
 * Compile with
 * gcc -Wall -g -fopenmp -o $(NAME_PART) $(FILE_NAME) -lm
 */

#include <omp.h>
#include <stdio.h>
#include <omp.h>

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
	*iterations = 0;
	
	for (i = 0; i < n; i++) {
		u[i] = 0.0;
	}
	
	for (k = 0; ; k++) {
		double maxError = 0.0;
		
		for (j = 0; j < n; j++) {
			sum = h * h * f[j];
			// for (i = 0; i < j; i++) {
				// sum -= A(i, j) * u[i];
			// }
			// for (i = j + 1; i < n; i++) {
				// sum -= A(i, j) * u[i];
			// }
			if (j - D >= 0) sum += u[j - D];
			if (j + D < n) sum += u[j + D];
			if (j - 1 >= 0 && (j - 1) / D == j / D) sum += u[j - 1];
			if (j + 1 < n && (j + 1) / D == j / D) sum += u[j + 1];
			sum /= 4;
			
			double error = sum - u[j];
			error = error > 0 ? error : -error;
			maxError = error > maxError ? error : maxError;
			
			u[j] = sum;
		}
		
		(*iterations)++;
		
		if (maxError < EPSILON) break;
	}
}

int main(void) {
	int i, j;
	
	double f[N];
	double u[N];
	double h = 1.0 / (D + 1);
	
	
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
