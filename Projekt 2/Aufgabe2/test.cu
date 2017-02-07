#include <cuda.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>

__global__ void Memcpy(int N,int* target,int* source){
	int i=blockDim.x*blockIdx.x+threadIdx.x;
	if (i<N) target[i]=source[i];
}
__global__ void Memcpyadd(int N,int* target,int*source,int add){
	int i=blockDim.x*blockIdx.x+threadIdx.x;
	if (i<N) target[i]=source[i]+add;
}
#define Mega 1000000
#define Kilo 1000
#define add 10
int main(){
	size_t N=__N__;
	size_t size=N*sizeof(int);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float t,r;
	
	int* c_a=(int*)malloc(size);
	#pragma omp parallel
	{
		srand(int(time(NULL)^omp_get_thread_num()));
		#pragma omp for
		for(int i=0;i<N;i++){
			c_a[i]=rand();
		}
	}	
	int* d_a;
	cudaMalloc(&d_a,size);
	int *c_d=(int*)malloc(size);
	cudaEventRecord(start);

	#pragma omp parallel for
	for(int i=0;i<N;i++){
		c_d[i]=c_a[i]+add;
		//c_d[i]=c_a[i];
	}

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&t, start, stop);
	t/=1000.0;
	r=(size/(1000.0*Mega))/t;
	printf("copy in Ram time: %f,rate: %f GB/s\n",t,r);
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
	dim3 BlockDim=dim3(1024,1,1);
	dim3 GridDim=dim3(N/1024+(((N %1024) == 0) ? 0 : 1),1,1);

	cudaEventRecord(start);
	//copy d_a in d_a
//	cudaMemcpy(d_b,d_a,size,cudaMemcpyDeviceToDevice);
	Memcpyadd<<<GridDim,BlockDim>>>(N,d_b,d_a,add);
	cudaDeviceSynchronize();
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
   		if(c_c[i]!=c_d[i]){
		printf("error in c_c at %d\n",i);
		break;
		}
   	}
}
