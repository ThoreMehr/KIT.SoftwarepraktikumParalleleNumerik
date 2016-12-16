#include <omp.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
	long num_steps =1;
	int j;
	for(j =0 ;j<atoi(argv[1]);j++)
		num_steps*=10;
	long i;
	double x, pi, sum, step,t;
	
	struct timeval time1,time2;
	struct timezone zone;
	gettimeofday(&time1,&zone);
	sum = 0.0;
	step = 1.0 / (double) num_steps;
	for (i = 1; i <= num_steps; i++) {
		x = (i-0.5) * step;
		sum +=4.0/(1.0+x*x);
	}
	pi = step * sum;
	gettimeofday(&time2,&zone);
	printf("PI: %f\n", pi);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;	
	printf("time:%f s\n",t);

//para1
	gettimeofday(&time1,&zone);
	sum = 0.0;
	step = 1.0 / (double) num_steps;
	#pragma omp parallel for reduction(+:sum) private(x)
	for (i = 1; i <= num_steps; i++) {
		x = (i-0.5) * step;
		sum +=4.0/(1.0+x*x);
	}
	pi = step * sum;
	gettimeofday(&time2,&zone);
	printf("PI: %f\n", pi);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;	
	printf("time:%f s\n",t);

//para2
	gettimeofday(&time1,&zone);
	sum = 0.0;
	step = 1.0 / (double) num_steps;

	#pragma omp parallel for private(x)
	for (i = 1; i <= num_steps; i++) {
		x = (i-0.5) * step;
		#pragma omp critical
		{
			sum += 4.0/(1.0+x*x);
		}
	}
	pi = step * sum;
	
	gettimeofday(&time2,&zone);
	printf("PI: %f\n", pi);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;	
	printf("time:%f s\n",t);

	return 0;
}
