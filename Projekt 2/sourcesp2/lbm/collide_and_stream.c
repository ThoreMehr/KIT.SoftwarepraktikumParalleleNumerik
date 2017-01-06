#include "global.h"

void collide_and_stream (int x0, int x1, int y0, int y1, int z0, int z1)
{
   int iX, iY, iZ;
   double**** fTmp;
   double ro, rr, uxP, uxM, uyP, uyM, uzP, uzM;
   double rho, ux, uy, uz, square, c_u, f_eq, tau_inv=1/tau;
   //Boundary conditions at top z=z1
   for (iX=x0+1;iX<=x1-1;++iX){
   	for (iY=y0+1;iY<=y1-1;++iY){
		rho=1;
                ux = u_0[iX][iY];
                square = 1.5*ux*ux;
                f_eq = rho/36.0 * (1. - square);
		gridI[iX][iY][z1][3]  = 2.0*f_eq; // 3->12
		gridI[iX][iY][z1][8]  = f_eq;     // 8->17
		gridI[iX][iY][z1][18] = f_eq;      //18->9
		f_eq = rho/36.0 * (1. - 3.*ux + 2.0*square);
                gridI[iX][iY][z1][6]  = f_eq;  
                f_eq = rho/36.0 * (1. + 3.*ux + 2.0*square);
		gridI[iX][iY][z1][16] = f_eq;  
		
                
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // bulk: collide and stream
    for (iX=x0+1; iX<=x1-1; ++iX) {
        for (iY=y0+1; iY<=y1-1; ++iY) {
            for (iZ=z0+1; iZ<=z1-1; ++iZ) {
                // collide
                ro= gridI[iX][iY][iZ][0]  + gridI[iX][iY][iZ][2]  + gridI[iX][iY][iZ][3]  + \
                    gridI[iX][iY][iZ][8]  + gridI[iX][iY][iZ][9]  + gridI[iX][iY][iZ][11] + \
                    gridI[iX][iY][iZ][12] + gridI[iX][iY][iZ][17] + gridI[iX][iY][iZ][18];
                uxM=gridI[iX][iY][iZ][1]  + gridI[iX][iY][iZ][4]  + gridI[iX][iY][iZ][5]  + \
                    gridI[iX][iY][iZ][6]  + gridI[iX][iY][iZ][7];
		uxP=gridI[iX][iY][iZ][10] + gridI[iX][iY][iZ][13] + gridI[iX][iY][iZ][14] + \
                    gridI[iX][iY][iZ][15] + gridI[iX][iY][iZ][16];
		uyM=gridI[iX][iY][iZ][2]  + gridI[iX][iY][iZ][4]  + gridI[iX][iY][iZ][8]  + \
                    gridI[iX][iY][iZ][9]  + gridI[iX][iY][iZ][14];
		uyP=gridI[iX][iY][iZ][5]  + gridI[iX][iY][iZ][11] + gridI[iX][iY][iZ][13] + \
                    gridI[iX][iY][iZ][17] + gridI[iX][iY][iZ][18];
		uzM=gridI[iX][iY][iZ][3]  + gridI[iX][iY][iZ][6]  + gridI[iX][iY][iZ][8]  + \
                    gridI[iX][iY][iZ][16] + gridI[iX][iY][iZ][18];
		uzP=gridI[iX][iY][iZ][7]  + gridI[iX][iY][iZ][9]  + gridI[iX][iY][iZ][12] + \
                    gridI[iX][iY][iZ][15] + gridI[iX][iY][iZ][17];

		rho= ro+uxP+uxM;
                ux = (uxP-uxM)/rho;
                uy = (uyP-uyM)/rho;
                uz = (uzP-uzM)/rho;
		square = 1.5 * (ux*ux + uy*uy + uz*uz);
		ro=rho/3.0;
                f_eq = ro * (1. - square);
		gridI[iX][iY][iZ][0]  = gridI[iX][iY][iZ][0]  + (f_eq - gridI[iX][iY][iZ][0])*tau_inv;
                rr=ro;
                ro=ro/6.0;
                f_eq = ro * (1. - 3.*ux + 4.5*ux*ux - square);
		gridI[iX][iY][iZ][1]  = gridI[iX][iY][iZ][1]  + (f_eq - gridI[iX][iY][iZ][1])*tau_inv;
                f_eq = f_eq + rr*ux;	
                gridI[iX][iY][iZ][10] = gridI[iX][iY][iZ][10] + (f_eq - gridI[iX][iY][iZ][10])*tau_inv;
                f_eq = ro * (1. - 3.*uy + 4.5*uy*uy - square);
		gridI[iX][iY][iZ][2]  = gridI[iX][iY][iZ][2]  + (f_eq - gridI[iX][iY][iZ][2])*tau_inv;
		f_eq = f_eq + rr*uy;
                gridI[iX][iY][iZ][11] = gridI[iX][iY][iZ][11] + (f_eq - gridI[iX][iY][iZ][11])*tau_inv;
                f_eq = ro * (1. - 3.*uz + 4.5*uz*uz - square);
		gridI[iX][iY][iZ][3]  = gridI[iX][iY][iZ][3]  + (f_eq - gridI[iX][iY][iZ][3])*tau_inv;
		f_eq = f_eq + rr*uz;
		gridI[iX][iY][iZ][12] = gridI[iX][iY][iZ][12] + (f_eq - gridI[iX][iY][iZ][12])*tau_inv;
                rr=rr/2.0;
                ro=ro/2.0;
                c_u = ux + uy;
		f_eq = ro * (1. - 3.*c_u + 4.5*c_u*c_u - square);
		gridI[iX][iY][iZ][4]  = gridI[iX][iY][iZ][4]  + (f_eq - gridI[iX][iY][iZ][4])*tau_inv;
		f_eq = f_eq + rr*c_u;
		gridI[iX][iY][iZ][13] = gridI[iX][iY][iZ][13] + (f_eq - gridI[iX][iY][iZ][13])*tau_inv;
                c_u = ux - uy;
		f_eq = ro * (1. - 3.*c_u + 4.5*c_u*c_u - square);
		gridI[iX][iY][iZ][5]  = gridI[iX][iY][iZ][5]  + (f_eq - gridI[iX][iY][iZ][5])*tau_inv;
		f_eq = f_eq + rr*c_u;
		gridI[iX][iY][iZ][14] = gridI[iX][iY][iZ][14] + (f_eq - gridI[iX][iY][iZ][14])*tau_inv;
                c_u = ux + uz;
		f_eq = ro * (1. - 3.*c_u + 4.5*c_u*c_u - square);
		gridI[iX][iY][iZ][6]  = gridI[iX][iY][iZ][6]  + (f_eq - gridI[iX][iY][iZ][6])*tau_inv;
		f_eq = f_eq + rr*c_u;
		gridI[iX][iY][iZ][15] = gridI[iX][iY][iZ][15] + (f_eq - gridI[iX][iY][iZ][15])*tau_inv;
                c_u = ux - uz;
		f_eq = ro * (1. - 3.*c_u + 4.5*c_u*c_u - square);
		gridI[iX][iY][iZ][7]  = gridI[iX][iY][iZ][7]  + (f_eq - gridI[iX][iY][iZ][7])*tau_inv;
		f_eq = f_eq + rr*c_u;
		gridI[iX][iY][iZ][16] = gridI[iX][iY][iZ][16] + (f_eq - gridI[iX][iY][iZ][16])*tau_inv;
                c_u = uy + uz;
		f_eq = ro * (1. - 3.*c_u + 4.5*c_u*c_u - square);
		gridI[iX][iY][iZ][8]  = gridI[iX][iY][iZ][8]  + (f_eq - gridI[iX][iY][iZ][8])*tau_inv;
		f_eq = f_eq + rr*c_u;
		gridI[iX][iY][iZ][17] = gridI[iX][iY][iZ][17] + (f_eq - gridI[iX][iY][iZ][17])*tau_inv;
                c_u = uy - uz;
		f_eq = ro * (1. - 3.*c_u + 4.5*c_u*c_u - square);
		gridI[iX][iY][iZ][9]  = gridI[iX][iY][iZ][9]  + (f_eq - gridI[iX][iY][iZ][9])*tau_inv;
		f_eq = f_eq + rr*c_u;
		gridI[iX][iY][iZ][18] = gridI[iX][iY][iZ][18] + (f_eq - gridI[iX][iY][iZ][18])*tau_inv;
                // stream
		gridO[iX][iY][iZ][0]      = gridI[iX][iY][iZ][0];
                gridO[iX-1][iY][iZ][1]    = gridI[iX][iY][iZ][1];
                gridO[iX][iY-1][iZ][2]    = gridI[iX][iY][iZ][2];
                gridO[iX][iY][iZ-1][3]    = gridI[iX][iY][iZ][3];
                gridO[iX-1][iY-1][iZ][4]  = gridI[iX][iY][iZ][4];
                gridO[iX-1][iY+1][iZ][5]  = gridI[iX][iY][iZ][5];
                gridO[iX-1][iY][iZ-1][6]  = gridI[iX][iY][iZ][6];
                gridO[iX-1][iY][iZ+1][7]  = gridI[iX][iY][iZ][7];
                gridO[iX][iY-1][iZ-1][8]  = gridI[iX][iY][iZ][8];
                gridO[iX][iY-1][iZ+1][9]  = gridI[iX][iY][iZ][9];
		gridO[iX+1][iY][iZ][10]   = gridI[iX][iY][iZ][10];
                gridO[iX][iY+1][iZ][11]   = gridI[iX][iY][iZ][11];
                gridO[iX][iY][iZ+1][12]   = gridI[iX][iY][iZ][12];
                gridO[iX+1][iY+1][iZ][13] = gridI[iX][iY][iZ][13];
                gridO[iX+1][iY-1][iZ][14] = gridI[iX][iY][iZ][14];
                gridO[iX+1][iY][iZ+1][15] = gridI[iX][iY][iZ][15];
                gridO[iX+1][iY][iZ-1][16] = gridI[iX][iY][iZ][16];
                gridO[iX][iY+1][iZ+1][17] = gridI[iX][iY][iZ][17];
                gridO[iX][iY+1][iZ-1][18] = gridI[iX][iY][iZ][18];
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////
    // boundary: stream
    for (iX=x0+1; iX<=x1-1; ++iX){
	for (iY=y0+1; iY<=y1-1; ++iY){
                gridO[iX-1][iY][z0+1][7]  = gridI[iX][iY][z0][16];
                gridO[iX][iY-1][z0+1][9]  = gridI[iX][iY][z0][18];
                gridO[iX][iY][z0+1][12]   = gridI[iX][iY][z0][3];
                gridO[iX+1][iY][z0+1][15] = gridI[iX][iY][z0][6];
                gridO[iX][iY+1][z0+1][17] = gridI[iX][iY][z0][8];
	}
    }
    for (iX=x0+1; iX<=x1-1; ++iX){
	for (iY=y0+1; iY<=y1-1; ++iY){
                gridO[iX][iY][z1-1][3]    = gridI[iX][iY][z1][3];
                gridO[iX-1][iY][z1-1][6]  = gridI[iX][iY][z1][6];
                gridO[iX][iY-1][z1-1][8]  = gridI[iX][iY][z1][8];
                gridO[iX+1][iY][z1-1][16] = gridI[iX][iY][z1][16];
                gridO[iX][iY+1][z1-1][18] = gridI[iX][iY][z1][18];
	}
    }
    for (iX=x0+1; iX<=x1-1; ++iX){
	for (iZ=z0+1; iZ<=z1-1; ++iZ){
                gridO[iX-1][y0+1][iZ][5]  = gridI[iX][y0][iZ][14];
                gridO[iX][y0+1][iZ][11]   = gridI[iX][y0][iZ][2];
                gridO[iX+1][y0+1][iZ][13] = gridI[iX][y0][iZ][4];
                gridO[iX][y0+1][iZ+1][17] = gridI[iX][y0][iZ][8];
                gridO[iX][y0+1][iZ-1][18] = gridI[iX][y0][iZ][9];
	}
    }
    for (iX=x0+1; iX<=x1-1; ++iX){
	for (iZ=z0+1; iZ<=z1-1; ++iZ){
                gridO[iX][y1-1][iZ][2]    = gridI[iX][y1][iZ][11];
                gridO[iX-1][y1-1][iZ][4]  = gridI[iX][y1][iZ][13];
                gridO[iX][y1-1][iZ-1][8]  = gridI[iX][y1][iZ][17];
                gridO[iX][y1-1][iZ+1][9]  = gridI[iX][y1][iZ][18];
                gridO[iX+1][y1-1][iZ][14] = gridI[iX][y1][iZ][5];
	}
    }
    for (iY=y0+1; iY<=y1-1; ++iY){
	for (iZ=z0+1; iZ<=z1-1; ++iZ){
		gridO[x0+1][iY][iZ][10]   = gridI[x0][iY][iZ][1];
                gridO[x0+1][iY+1][iZ][13] = gridI[x0][iY][iZ][4];
                gridO[x0+1][iY-1][iZ][14] = gridI[x0][iY][iZ][5];
                gridO[x0+1][iY][iZ+1][15] = gridI[x0][iY][iZ][6];
                gridO[x0+1][iY][iZ-1][16] = gridI[x0][iY][iZ][7];
	}
    }
    for (iY=y0+1; iY<=y1-1; ++iY){
	for (iZ=z0+1; iZ<=z1-1; ++iZ){
                gridO[x1-1][iY][iZ][1]    = gridI[x1][iY][iZ][10];
                gridO[x1-1][iY-1][iZ][4]  = gridI[x1][iY][iZ][13];
                gridO[x1-1][iY+1][iZ][5]  = gridI[x1][iY][iZ][14];
                gridO[x1-1][iY][iZ-1][6]  = gridI[x1][iY][iZ][15];
                gridO[x1-1][iY][iZ+1][7]  = gridI[x1][iY][iZ][16];
	}
    }
    for (iX=x0+1; iX<=x1-1; iX++){
                gridO[iX][y1-1][z1-1][8]  = gridI[iX][y1][z1][17];
                gridO[iX][y1-1][z0+1][9]  = gridI[iX][y1][z0][18];
                gridO[iX][y0+1][z0+1][17] = gridI[iX][y0][z0][8];
                gridO[iX][y0+1][z1-1][18] = gridI[iX][y0][z1][9];
    }
    for (iY=y0+1; iY<=y1-1; iY++){
                gridO[x1-1][iY][z1-1][6]  = gridI[x1][iY][z1][15];
                gridO[x1-1][iY][z0+1][7]  = gridI[x1][iY][z0][16];
                gridO[x0+1][iY][z0+1][15] = gridI[x0][iY][z0][6];
                gridO[x0+1][iY][z1-1][16] = gridI[x0][iY][z1][7];
    }
    for (iZ=z0+1; iZ<=z1-1; iZ++){
                gridO[x1-1][y1-1][iZ][4]  = gridI[x1][y1][iZ][13];
                gridO[x1-1][y0+1][iZ][5]  = gridI[x1][y0][iZ][14];
                gridO[x0+1][y0+1][iZ][13] = gridI[x0][y0][iZ][4];
                gridO[x0+1][y1-1][iZ][14] = gridI[x0][y1][iZ][5];
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    fTmp = gridO; gridO=gridI; gridI=fTmp;
}

