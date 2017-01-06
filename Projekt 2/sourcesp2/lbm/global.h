#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* prototype declarations */
int  allocate();
void free_memory();
void init();
void collide_and_stream(int,int,int,int,int,int);
void set_u();

/* global variables */
int max_x, max_y, max_z;                    /* number of grid nodes in each direction */
double ****gridI;			    /* 4d matrix made up of 19 3D matrices */
double ****gridO;                           /* 4d matrix made up of 19 3D matrices */
double *rawdata1, *rawdata2;
double tau, nu;                             /* collision-relaxation time, shear viscosity */
double **u_0;                               /* values of u_0 */
double delta_x;                             /* grid fineness */
double max_u;                               /* maximal velocity */
double rho_0;

// static int c[57] = {
//              0, 0, 0,
//             -1, 0, 0,  0,-1, 0,  0, 0,-1,
//             -1,-1, 0, -1, 1, 0, -1, 0,-1,
//             -1, 0, 1,  0,-1,-1,  0,-1, 1,
// 
//              1, 0, 0,  0, 1, 0,  0, 0, 1,
//              1, 1, 0,  1,-1, 0,  1, 0, 1,
//              1, 0,-1,  0, 1, 1,  0, 1,-1
// };


