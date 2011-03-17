#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const unsigned short dim = 2;	// Dimension of the state vector
const unsigned long bisect_limit = 20;
const unsigned long step_limit = 1000;	// Maximum number of steps to try on a single shot
const double delta_x = 0.1;		// Step size used in integration
char *outfile_name = "output.dat";	// Some output file

#include "./shoot.c"
#include "./pendulum.c"

double zeta[] = {0,1};	// The interval we wish to perform bisection on
void (*derivative)(double x, double *y, double *deriv_ret) = &f_pendulum;
int (*hit_func)(double x, double *y) = &hit_pendulum;
void (*init_func)(double *x, double *y, double zeta) = &init_pendulum;

int main(int argc, char *argv[]){
	unsigned long i;	// Dummy counter variable

	FILE *outfile = fopen(outfile_name,"w");
	if(outfile==NULL) {fprintf(stderr, "Unable to open %s for writing. Quitting.\n",outfile_name); return -1;}

	/* Set up and check that the interval is appropriate for bisection */
	double x_0, y_0[dim];
	int miss_sign[] = {0,0};
	struct hit_data HIT;
	for(i=0;i<2;i++){
		init_func(&x_0,y_0,zeta[i]);
		HIT = shoot(x_0,y_0,delta_x,step_limit,NULL,derivative,hit_func);
		if(HIT.hit) {break;}
		else {miss_sign[i] = HIT.miss_direction;}
	}
	if (miss_sign[0]==miss_sign[2]) {printf("Error! The miss directions from %lf and %lf are both %i.\n",zeta[0],zeta[1],miss_sign[0]); return 1;}	// Did we make a bad choice for the ends of the interval? This should only be possible on the first loop

	/* Start of interval bisection */
	for (i=0;i<bisect_limit;i++){
		printf("[%lf,%lf]\n",zeta[0],zeta[1]);
		init_func(&x_0,y_0,0.5*(zeta[0]+zeta[1]));
		HIT = shoot(x_0,y_0,delta_x,step_limit,NULL,derivative,hit_func);
		if(HIT.hit) {printf("Hit with parameter choice zeta = %lf\n",0.5*(zeta[0]+zeta[1])); break;}
		else if(HIT.miss_direction*(zeta[1]-zeta[0])>0) {zeta[1] = 0.5*(zeta[0]+zeta[1]);}
		else {zeta[0] = 0.5*(zeta[0]+zeta[1]);}
	}
	/* End of interval bisection */

	/* Perform a final shot using either the zeta that hit or the last zeta that we tried */
	init_func(&x_0,y_0,zeta[1]);
	shoot(x_0,y_0,delta_x,HIT.hit_steps,outfile,derivative,hit_func);
	if(HIT.hit==0) {printf("Failed to hit. Trajectory zeta = %lf has been written.\n",zeta[1]); return 1;}
	
	return 0;
}
