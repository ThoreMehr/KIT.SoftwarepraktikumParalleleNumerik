
/*
 * Compilation:
 * 		w/ image output: gcc -Wall -o mandelbrot -D IMAGE_OUTPUT mandelbrot.c -lm
 * 		w/o image output: gcc -Wall -o mandelbrot mandelbrot.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#ifndef N
	#define N 5000 /* NxN grid size */
#endif
#ifndef maxiter
	#define maxiter 500
#endif

typedef struct {
	float re, im;
} complex;

float cabs1(complex x) {
	return sqrt(x.re*x.re+x.im*x.im);
}

complex add(complex x, complex y) {
	complex r;
	r.re=x.re+y.re; r.im=x.im+y.im;
	return r;
}

complex mult(complex x, complex y) {
	complex r;
	r.re = x.re*y.re-x.im*y.im;
	r.im = x.re*y.im+x.im*y.re;
	return r;
}

int main() {
	complex z, kappa;
	int i, j, k;
	int *T;
	
	T = (int*) malloc(sizeof(int)*N*N);
	float tp;
	printf("Starting calculation for N=%d... with %d threads", N,omp_get_max_threads());
	struct timeval time1,time2;
	struct timezone zone;	
	gettimeofday(&time1,&zone);
	#ifndef scheduleType
		#pragma omp parallel for private(j,z,kappa,k) shared(T)
	#else
		#ifndef chunksize
	 		#pragma omp parallel for private(j,z,kappa,k) shared(T) schedule(scheduleType)
		#else
			#pragma omp parallel for private(j,z,kappa,k) shared(T) schedule(scheduleType,chunksize)
		#endif
	#endif
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			z.re = kappa.re = (4.0*(i-N/2))/N;
			z.im = kappa.im = (4.0*(j-N/2))/N;
			
			for (k=0; ; k++) {
				if (cabs1(z)>2 || k==maxiter) {
					T[i*N+j]=(k*256)/maxiter;
					break;
				}
				z = add( mult(z,z), kappa);
			}
		}
	}
	gettimeofday(&time2,&zone);
	tp=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;	
	printf("time:%f s\n",tp);
	#ifdef IMAGE_OUTPUT
	printf("Writing simple image file...\n");
	
	FILE *f = fopen("output.ppm", "w");
	fprintf(f, "P3 %d %d 255 ", N, N);
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			fprintf(f, "%d %d %d ", T[i*N+j]%256, T[i*N+j]%256, T[i*N+j]%256);
		}
	}
	fclose(f);
	#endif
	
	free(T);
	
	return 0;
}
