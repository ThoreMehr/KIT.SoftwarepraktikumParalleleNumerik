#include "global.h"

int allocate () {
	int iX,iY,iZ, test=1;

	if(!(rawdata1 = (double *)calloc(max_x*max_y*max_z*19,sizeof(double)))) test=0;
	if(!(rawdata2 = (double *)calloc(max_x*max_y*max_z*19,sizeof(double)))) test=0;

	if(!(gridI = (double ****)calloc(max_x,sizeof(double ***)))) test=0;
	for(iX=0;iX<max_x;iX++) {
		if(!(gridI[iX] = (double ***)calloc(max_y,sizeof(double **)))) test=0;
		for(iY=0;iY<max_y;iY++) {
			if(!(gridI[iX][iY] = (double **)calloc(max_z,sizeof(double*)))) test=0;
			for(iZ=0;iZ<max_z;iZ++) {
				gridI[iX][iY][iZ] = rawdata1 + 19*(max_z* (max_y*iX+iY) + iZ);
			}
		}
	}

	if(!(gridO = (double ****)calloc(max_x,sizeof(double ***)))) test=0;
	for(iX=0;iX<max_x;iX++) {
		if(!(gridO[iX] = (double ***)calloc(max_y,sizeof(double **)))) test=0;
		for(iY=0;iY<max_y;iY++) {
			if(!(gridO[iX][iY] = (double **)calloc(max_z,sizeof(double*)))) test=0;
			for(iZ=0;iZ<max_z;iZ++) {
				gridO[iX][iY][iZ] = rawdata2 + 19*(max_z* (max_y*iX+iY) + iZ);
			}
		}
	}

	if(!(u_0 = (double **)calloc(max_x,sizeof(double *)))) test=0; 
	for(iX=0;iX<max_x;iX++) {
		if(!(u_0[iX] = (double *)calloc(max_y,sizeof(double)))) test=0;
	}

	if(!test) {
		printf("OUT OF MEMORY. test==BAD\n");
		exit(1);
	}

	return 0;
}
