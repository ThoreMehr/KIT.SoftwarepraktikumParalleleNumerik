#include <cuda.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>

__global__ void Memcpy(int* target,int* source){
	int i=threadIdx.x;
	target[i]=source[i];
}
__global__ void Memcpyadd(int* target,int*source,int add){
	int i=threadIdx.x;
	target[i]=source[i]+add;
}
#define Mega 1000000
#define Kilo 1000
int main(){
	srand(time(NULL));
	size_t N=100*Mega;
	size_t size=N*sizeof(int);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float t,r;
	
	int* c_a=(int*)malloc(size);
	for(int i=0;i<N;i++){
		c_a[i]=rand()%10000;
	}	
	int* d_a;
	cudaMalloc(&d_a,size);

	cudaEventRecord(start);

	//copy c_a to Device
	cudaMemcpy(d_a,c_a,size,cudaMemcpyHostToDevice);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t, start, stop);
	t/=1000.0;
	r=(size/(1000.0*Mega))/t;
	printf("copy to device time: %f,rate: %f GB/s\n",t,r);

	//c_b for testing	
	int* c_b=(int*)malloc(size);
	cudaMemcpy(c_b,d_a,size,cudaMemcpyDeviceToHost);

	int* d_b;
	cudaMalloc(&d_b,size);

	cudaEventRecord(start);
	//copy d_a in d_b
	cudaMemcpy(d_b,d_a,size,cudaMemcpyDeviceToDevice);

	//Memcpy<<<1,N>>>(d_b,d_a);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t, start, stop);
	t/=1000.0;
	r=(size/(1000.0*Mega))/t;
	printf("copy on device time: %f,rate: %f GB/s\n",t,r);

   	int* c_c=(int*)malloc(size);

	cudaEventRecord(start);
	//copy d_b from device to c_c
   	cudaMemcpy(c_c,d_b,size,cudaMemcpyDeviceToHost);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t, start, stop);
	t/=1000.0;
	r=(size/(1000.0*Mega))/t;
	printf("copy from device time: %f,rate: %f GB/s\n",t,r);
//checking correctness
   	for(int i=0;i<N;i++){
   		if(c_b[i]!=c_a[i]){
		printf("error in c_b at %d\n",i);
		break;
		}
   	}
   	for(int i=0;i<N;i++){
   		if(c_c[i]!=c_a[i]){
		printf("error in c_c at %d\n",i);
		break;
		}
   	}
}
