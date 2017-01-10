#include <cuda.h>
#include <stdio.h>

int main() {
	int count;
	cudaGetDeviceCount(&count);
	printf("Found %d CUDA devices\n", count);
	printf("=========================================\n");


	for (int device = 0; device < count; device++) {
		struct cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, device);
		
		printf("Information for %s (device %d):\n", prop.name, device);
		printf("  Total global  memory: %zd\n", prop.totalGlobalMem);
		printf("  Total const memory: %zd\n", prop.totalConstMem);
		printf("  Shared memory per block: %zd\n", prop.sharedMemPerBlock);
		printf("  Warp size: %d\n", prop.warpSize);
		printf("  Max threads per block: %d\n", prop.maxThreadsPerBlock);
		printf("  Max threads dimension: [ %d, %d, %d ]\n",
			prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
		printf("  Max grid size: [ %d, %d, %d ]\n",
			prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
		
		printf("-----------------------------------------\n");
	}
	
	return 0;
}
