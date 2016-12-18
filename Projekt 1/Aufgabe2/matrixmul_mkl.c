#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mkl.h>
#ifndef N
	#define N 500
#endif
	int main(int argc, char* argv[]) {
		//read the two matrixes from file
		FILE * pFile;
   		pFile = fopen (argv[1], "r");
   		//float m1[N][N];
   		float* m1=malloc(sizeof(float)*N*N);//output matrix
  		fread (m1 , sizeof(float), N*N, pFile);
  		fclose (pFile);
   		pFile = fopen (argv[2], "r");
//   		float m2[N][N];
		float* m2=malloc(sizeof(float)*N*N);//output matrix
  		fwrite (m2 , sizeof(float), N*N, pFile);
  		fclose (pFile);
  		float* m3=malloc(sizeof(float)*N*N);//output matrix
  		//float m3[N][N];
  		int i,j,k=0;
  		
		struct timeval time1,time2;
		struct timezone zone;	
		gettimeofday(&time1,&zone);
  		cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,1.0,m1,N,m2,N,0.0,m3,N);
    	gettimeofday(&time2,&zone);
    	
		float t=((time2.tv_usec-time1.tv_usec)+(time2.tv_sec-time1.tv_sec)*1000000)/1000000.0;	
		printf("time:%f s\n",t);
    
   		pFile = fopen (argv[3], "wb");
	  	fwrite (m3 , sizeof(float), N*N, pFile);
  		fclose (pFile);
  		free(m1);
		free(m2);
		free(m3);	
		
}

