#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const unsigned short dim = 8;	// Dimension of the state vector
const unsigned long bisect_limit = 200;
const unsigned long step_limit = 10000;	// Maximum number of steps to try on a single shot
const double delta_x = 0.01;		// Step size used in integration
char *outfile_name = "output.dat";	// Some output file
char *outfileb_name = "output0.dat";	// Some output file
double zeta[] = {1.25,1.75};	// The interval we wish to perform bisection on. NB. This bears no relation to the variable "zeta" in the SUGRA system

#include "./shoot.c"
#include "./pendulum.c"
#include "./SUGRA.c"

void (*derivative)(double x, double *y, double *deriv_ret) = &f_SUGRA;
int (*hit_func)(double x, double *y) = &hit_AdS_small_zeta;
void (*init_func)(double *x, double *y, double zeta) = &init_Lifshitz;

int main(int argc, char *argv[]){
	unsigned long i;	// Dummy counter variable

	FILE *outfile = fopen(outfile_name,"w");
	if(outfile==NULL) {fprintf(stderr, "Unable to open %s for writing. Quitting.\n",outfile_name); return -1;}

	FILE *outfileb = fopen(outfileb_name,"w");
	if(outfile==NULL) {fprintf(stderr, "Unable to open %s for writing. Quitting.\n",outfileb_name); return -1;}

	/* Set up and check that the interval is appropriate for bisection */
	double x_0, y_0[dim];
	int miss_sign[] = {0,0};
	struct hit_data HIT;
	for(i=0;i<2;i++){
		init_func(&x_0,y_0,zeta[i]);
		HIT = shoot(x_0,y_0,delta_x,step_limit,NULL,derivative,hit_func);
		if(HIT.hit) {break;}
		else {miss_sign[i] = HIT.miss_direction; fprintf(stderr,"With zeta = %lf, missed by direction %d.\n",zeta[i],miss_sign[i]);}
	}
	if (miss_sign[0]==miss_sign[1]) {printf("Error! The miss directions from %lf and %lf are both %i.\n",zeta[0],zeta[1],miss_sign[0]); return 1;}	// Did we make a bad choice for the ends of the interval? This should only be possible on the first loop

	/* Start of interval bisection */
	for (i=0;i<bisect_limit;i++){
		if((zeta[1]==0.5*(zeta[1]+zeta[0]))||(zeta[0]==0.5*(zeta[1]+zeta[0]))) {printf("Bisection point is numerically the same as an end point!\n"); break;}
		printf("Interval number %ld:\t[%lf,%lf]\n",i,zeta[0],zeta[1]);
		init_func(&x_0,y_0,0.5*(zeta[0]+zeta[1]));
		HIT = shoot(x_0,y_0,delta_x,step_limit,NULL,derivative,hit_func);
		fprintf(stderr,"With zeta = %lf, missed by direction %d.\n",0.5*(zeta[0]+zeta[1]),HIT.miss_direction);
		if(HIT.hit) {printf("Hit with parameter choice zeta = %lf\n",0.5*(zeta[0]+zeta[1])); zeta[1] = 0.5*(zeta[1]+zeta[0]); break;}
		else if(HIT.miss_direction*(miss_sign[1]-miss_sign[0])>0) {zeta[1] = 0.5*(zeta[0]+zeta[1]);}
		else {zeta[0] = 0.5*(zeta[0]+zeta[1]);}
	}
	/* End of interval bisection */

	/* Perform a final shot using either the zeta that hit or the last zeta that we tried */
	init_func(&x_0,y_0,zeta[1]);
	shoot(x_0,y_0,delta_x,HIT.hit_steps,outfile,derivative,hit_func);
	if(HIT.hit==0) {
		printf("Failed to hit. Trajectory zeta = %lf has been written to %s.\n",zeta[1],outfile_name);
		init_func(&x_0,y_0,zeta[0]);
		shoot(x_0,y_0,delta_x,HIT.hit_steps,outfileb,derivative,hit_func);
		printf("Trajectory zeta = %lf has been written to %s,\n",zeta[0],outfileb_name);
		return 1;
	}

	fclose(outfile);
	fclose(outfileb);
	
	return 0;
}
