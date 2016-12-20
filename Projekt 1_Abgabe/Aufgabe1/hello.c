
/*
 * compile: gcc -Wall -g -fopenmp -o hello hello.c -lm
 * use: [set|export] OMP_NUM_THREADS=8
 */

#include <stdio.h>
#include <omp.h>

int main() {
	#pragma omp parallel
	{
		printf("Hello World, this is Thread%i\n", omp_get_thread_num());
	}
	return 0;
}

