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
	return 8 * M_PI * M_PI * sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

double analyticU(double x, double y) {
	return sin(2 * M_PI * x) * sin(2 * M_PI * y);
}

/*
 * Solves Ax = b for with A of size n*n and structure as in FEM, x, b of size n
 */
void solve(double *b, double *x, int *iterations) {
	int i, ix;
	double r[N];
	double p[N];
	double delta = 0, deltaHat;
	
	#pragma omp parallel for private(i)
	for (i = 0; i < N; i++) {
		x[i] = 1;
	}
	
	#pragma omp parallel for private(i, ix)
	for (i = 0; i < N; i++) {
		ix = i % D;
		double residuum = b[i];
		
		if (ix - 1 >= 0) residuum += x[i - 1];
		if (ix + 1 < D) residuum += x[i + 1];
		if (i - D >= 0) residuum += x[i - D]; 
		if (i + D < N) residuum += x[i + D];
		residuum -= 4 * x[i];
		
		r[i] = residuum;
		p[i] = r[i];
	}
	
	#pragma omp parallel for private(i) reduction(+:delta)
	for (i = 0; i < N; i++) {
		delta += r[i] * r[i];
	}
	
	(*iterations) = 0;
	while (delta >= EPSILON * EPSILON) {
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
		
		double newDelta = 0;
		#pragma omp parallel for private(i, ix) reduction(+:newDelta)
		for (i = 0; i < N; i++) {
			ix = i % D;
			double v = 0;
			if (ix - 1 >= 0) v += p[i - 1];
			if (ix + 1 < D) v += p[i + 1];
			if (i - D >= 0) v += p[i - D]; 
			if (i + D < N) v += p[i + D];
			v -= 4 * p[i];
			r[i] += deltaHat * v;
			x[i] += deltaHat * p[i];
			newDelta += r[i] * r[i];
		}
		
		delta = newDelta / delta;
		#pragma omp parallel for private(i)
		for (i = 0; i < N; i++) {
			p[i] = r[i] + delta * p[i];
		}
		
		delta = newDelta;
		(*iterations)++;
	}
	
}

int main(void) {
	int i, j;
	
	double f[N];
	double u[N];
	
	printf("Run on system with %d processors and max %d threads\n", omp_get_num_procs(), omp_get_max_threads());
	
	printf("l = %d\nd = %d\nn = %d\n", L, D, N);
	
	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++) {
			f[i + D * j] = H * H * func(i * H + H, j * H + H);
		}
	}
	
	int it;
	solve(f, u, &it);
	
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
