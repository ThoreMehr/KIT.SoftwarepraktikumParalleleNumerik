#include "global.h"

void init() {

	int iX,iY,iZ;

	/* initialise distributions with incompressible equilibrium distributions and u=0, rho=rho_0 */
	for(iX=0;iX<max_x;iX++) {
		for(iY=0;iY<max_y;iY++) {
			for(iZ=0;iZ<max_z;iZ++) {
				gridI[iX][iY][iZ][0] = rho_0/3;
				gridI[iX][iY][iZ][1] = rho_0/18;
				gridI[iX][iY][iZ][2] = rho_0/18;
				gridI[iX][iY][iZ][3] = rho_0/18;
				gridI[iX][iY][iZ][4] = rho_0/36;
				gridI[iX][iY][iZ][5] = rho_0/36;
				gridI[iX][iY][iZ][6] = rho_0/36;
				gridI[iX][iY][iZ][7] = rho_0/36;
				gridI[iX][iY][iZ][8] = rho_0/36;
				gridI[iX][iY][iZ][9] = rho_0/36;
				gridI[iX][iY][iZ][10] = rho_0/18;
				gridI[iX][iY][iZ][11] = rho_0/18; 
				gridI[iX][iY][iZ][12] = rho_0/18;
				gridI[iX][iY][iZ][13] = rho_0/36;
				gridI[iX][iY][iZ][14] = rho_0/36;
				gridI[iX][iY][iZ][15] = rho_0/36;
				gridI[iX][iY][iZ][16] = rho_0/36;
				gridI[iX][iY][iZ][17] = rho_0/36;
				gridI[iX][iY][iZ][18] = rho_0/36; 
			}
		}
	}
}
