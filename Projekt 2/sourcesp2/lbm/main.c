#include "global.h"

int main() {
	/* Lid Driven Cavity */
	int t,t_max;

	/* initialisation */	
	gridI = NULL;
	gridO = NULL;
	
	/* read parameters */  
	//printf("max_x, max_y, max_z, t_max, pp_int, nu, max_u ?\n");
	scanf("%d %d %d %d %lf %lf", &max_x, &max_y, &max_z, &t_max, &nu, &max_u);
	//printf("max_x=%d, max_y=%d, max_z=%d, t_max=%d, ",max_x,max_y,max_z,t_max);
	//printf("nu=%lf, max_u=%lf, Re= %lf\n",nu,max_u,max_u/nu );
	
	/* Compute parameters */
	rho_0=1.0;
	delta_x = 1.0/((double)(max_x-1));
	tau = (6./delta_x * nu + 1.) / 2.;

	/* allocate memory */
	allocate();

	/* initialise mass and momentum (u=v=0, rho_0, f=f_eq) */
	init();

	/* assign velocity u_0 */	
	set_u();

	/* time loop */
	for(t=0;t<t_max;t++) {
		collide_and_stream(0, max_x-1, 0, max_y-1, 0, max_z-1);
		//printf("Zeitschritt: %d\n",t);
	}

	/* Testausgabe */
// 	for (x=1;x<max_x-1;x++){
// 		for(y=1;y<max_y-1;y++){
// 			for(z=1;z<max_z-1;z++){
// 				for(i=0;i<19;i++){
//                                         printf("%d %d %d %d", x,y,z,i);
// 					printf(" %f\n",gridI[x][y][z][i]);
// 				}
// 				printf("\n");
// 			}
// 		}
// 	} 


	/* free memory */
	free_memory();

	/* finish */
	return 0;
}

