#include <cuda.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#define Mega 1000000

__global__ void Memcopy_add(int* A,int* B,int add){
	int i =threadIdx.x;
	B[i]=A[i]+add;
}

int main() {
	cudaSetDevice(0);
	srand(time(NULL));
	size_t size=__N__*sizeof(int);
	int* c_a=(int*)malloc(size);
	for(int i=0;i<__N__;i++){
		c_a[i]=rand()%1000;
	}
	printf("%d\n",c_a[1]);
	int* c_b=(int*)malloc(size);
	int* c_c=(int*)malloc(size);
	struct timeval time1,time2;
	struct timezone zone;
	float t=0.0,r=0.0;

	gettimeofday(&time1,&zone);
	for(int i=0;i<__N__;i++){
		c_b[i]=c_a[i];
	}
	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("RAM->RAM time:%f s, rate:%f GB/s\n",t,r);
	printf("%d\n",c_b[1]);
	int* d_a;
	cudaMalloc(&d_a,size);
	int* d_b;
	cudaMalloc(&d_b,size);
	
	gettimeofday(&time1,&zone);

	cudaMemcpy(d_a,c_a,size,cudaMemcpyHostToDevice);

	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("RAM->GPU time:%f s, rate:%f GB/s\n",t,r);

	gettimeofday(&time1,&zone);
	//cudaMemcpy(d_b,d_a,size,cudaMemcpyDeviceToDevice);
	Memcopy_add<<<1,__N__>>>(d_a,d_b,0);
	
	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("GPU->GPU time:%f s, rate:%f GB/s\n",t,r);
	gettimeofday(&time1,&zone);

	cudaMemcpy(c_c,d_b,size,cudaMemcpyDeviceToHost);

	gettimeofday(&time2,&zone);
	t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;
	r=(size/1000000000.0)/t;	
	printf("GPU->RAM time:%f s, rate:%f GB/s\n",t,r);
	printf("%d\n",c_c[1]);
	for(int i=0;i<__N__;i++){
		if(c_c[i]-2!=c_a[i]){
			printf("%d||%d\n",c_a[i],c_c[i]);
			printf("error at %d\n",i);
			break;
		}
	}
	cudaFree(d_a);
	cudaFree(d_b);
	return 0;
}
