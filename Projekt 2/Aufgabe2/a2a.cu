#include <cuda.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
int main() {
	cudaSetDevice(0);
	srand(time(NULL));
	size_t size=N*sizeof(int);
	int* c_a=(int*)malloc(size);
	for(int i=0;i<N;i++){
		c_a[i]=rand();
	}
	int* c_b=(int*)malloc(size);
	struct timeval time1,time2;
	struct timezone zone;
	float t=0.0,r=0.0;

	gettimeofday(&time1,&zone);
	for(int i=0;i<N;i++){
		c_b[i]=c_a[i];
	}

	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("time:%f s, rate:%f GB/s\n",t,r);

	int* d_a;
	cudaMalloc(&d_a,size);
	int* d_b;
	cudaMalloc(&d_b,size);
	
	gettimeofday(&time1,&zone);

	cudaMemcpy(d_a,c_a,size,cudaMemcpyHostToDevice);

	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("time:%f s, rate:%f GB/s\n",t,r);

	gettimeofday(&time1,&zone);
	cudaMemcpy(d_b,d_a,size,cudaMemcpyDeviceToDevice);
	
	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("time:%f s, rate:%f GB/s\n",t,r);

	gettimeofday(&time1,&zone);

	cudaMemcpy(c_b,d_b,size,cudaMemcpyDeviceToHost);

	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("time:%f s, rate:%f GB/s\n",t,r);

	cudaFree(d_a);
	cudaFree(d_b);
	return 0;
}
