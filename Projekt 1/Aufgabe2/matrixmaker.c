#include <stdlib.h>
#include <time.h>
#include <stdio.h>
	float randfloat(){
		srand(time(0));
	    double i = 0;
	    i = rand() % 1000 - 500; //Gives a number between -500 and +500;
	    return i / 100;
		}
	int main(int argc, char* argv[]) {
		printf("putting out to %s\n",argv[1]);
		long size=atoi(argv[2]);
		size*=size;
		float* matrix= malloc(sizeof(float)*size);
		int i=0;		
		for(i=0;i<size;i++){
			matrix[i]=randfloat();
		}
		FILE * pFile;
   		pFile = fopen (argv[1], "wb");
  		fwrite (matrix , sizeof(float), size, pFile);
  		fclose (pFile);
		
}

