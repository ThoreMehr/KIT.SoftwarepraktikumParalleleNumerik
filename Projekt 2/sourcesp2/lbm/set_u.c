#include "global.h"

void set_u () {

	int iX,iY;

	for(iX=0;iX<max_x;iX++) {
		for(iY=0;iY<max_y;iY++) {
			u_0[iX][iY] = max_u;
		}
	}
}
