
#include <omp.h>
#include <stdio.h>

static long num_steps = 1000000000;

int main() {
	long i;
	double x, pi, sum, step;
	
	sum = 0.0;
	step = 1.0 / (double) num_steps;
	
	#pragma omp parallel for private(x)
	for (i = 1; i <= num_steps; i++) {
		x = (i-0.5) * step;
		#pragma omp atomic
		sum = sum + 4.0/(1.0+x*x);
	}
	pi = step * sum;
	
	printf("PI: %f\n", pi);
	
	return 0;
}