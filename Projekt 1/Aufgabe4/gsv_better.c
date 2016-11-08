
#include <omp.h>
#include <stdio.h>
#include <omp.h>

#define L 7
#define D ((1<<L)-1)
#define N (D*D)
#define EPSILON 1e-6

double func(double x, double y) {
	return 32 * (x * (1 - x) + y * (1 - y));
}

double analyticU(double x, double y) {
	return 16 * x * (1 - x) * y * (1 - y);
}

void solve(double h, double *f, double *u) {
	int i, j, k;
	double sum;
	double lastU[N];
	
	for (i = 0; i < N; i++) {
		u[i] = 0.0;
	}
	
	for (k = 0; ; k++) {
		#pragma omp parallel for
		for (j = 0; j < N; j++) {
			lastU[j] = u[j];
		}
		
		double error = 0.0;
		#pragma omp parallel for private(sum) reduction(max:error)
		for (j = 0; j < N; j++) {
			sum = h * h * f[j];
			
			if (j - D >= 0) sum += lastU[j - D];
			if (j + D < N) sum += lastU[j + D];
			if (j - 1 >= 0 && (j - 1) / D == j / D) sum += lastU[j - 1];
			if (j + 1 < N && (j + 1) / D == j / D) sum += lastU[j + 1];
			sum /= 4;
			
			error = sum - u[j];
			error = error > 0 ? error : -error;
			
			u[j] = sum;
		}
		
		if (error < EPSILON) break;
	}
	
	printf("%d iterations\n", k);
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
	
	solve(h, f, u);
	
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
	
	double maxError = 0.0;
	for (i = 0; i < D; i++) {
		for (j = 0; j < D; j++) {
			double error = analyticU(j * h + h, i * h + h) - u[j + D * i];
			error = error > 0 ? error : -error;
			maxError = error > maxError ? error : maxError;
		}
	}
	printf("Max error: %4.8f", maxError);
	
	return 0;
}
