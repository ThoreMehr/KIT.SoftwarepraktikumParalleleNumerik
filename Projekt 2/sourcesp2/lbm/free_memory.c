#include "global.h"

void free_memory() {
        int iX, iY;

        free(rawdata1);
        free(rawdata2);
	for(iX=0;iX<max_x;iX++) {
		for(iY=0;iY<max_y;iY++)
			free(gridI[iX][iY]);
		free(gridI[iX]);
	}
        for(iX=0;iX<max_x;iX++) {
		for(iY=0;iY<max_y;iY++)
			free(gridO[iX][iY]); 
		free(gridO[iX]);
	}
	free(gridI);
	free(gridO);
	free(u_0);
}
