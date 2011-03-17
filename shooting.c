#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const unsigned short dim = 2;	// Dimension of the state vector
const unsigned long bisect_limit = 20;
char *outfile_name = "output.dat";	// Some output file

const double pi = M_PI;
struct hit_data {int hit; int miss_direction; double hit_time;};	// Data type to allow the return of data about hitting

struct hit_data shoot(double x_0, double *y_0, double delta_x, unsigned long x_steps, FILE *print, void (*derivative)(double x, double *y, double *deriv_ret), int (*hit_func)(double x, double *y));

/* Function declarations for derivatives */
void f(double x, double *y, double *deriv_ret);
void f_pendulum(double x, double *y, double *deriv_ret);
/* End of function declarations for derivatives */

/* Function declarations for initialisations */
void init_pendulum(double *x, double *y, double zeta);
/* End of function declarations for initialisations */

/* Function declarations for hit testing */
int hit_pendulum(double x, double *y);
/* End of function declarations for hit testing */

int main(int argc, char *argv[]){
	unsigned long i;	// Dummy counter variable


	FILE *outfile = fopen(outfile_name,"w");
	if(outfile==NULL) {fprintf(stderr, "Unable to open %s for writing. Quitting.\n",outfile_name); return -1;}

	
	double x_0, y_0[dim];
	double zeta[] = {0,1};	// The interval we wish to perform bisection on
	int miss_sign[] = {0,0};
	struct hit_data HIT;
	for(i=0;i<2;i++){
		init_pendulum(&x_0,y_0,zeta[i]);
		HIT = shoot(x_0,y_0,0.1,1000,NULL,&f_pendulum,&hit_pendulum);
		if(HIT.hit) {break;}
		else {miss_sign[i] = HIT.miss_direction;}
	}
	if (miss_sign[0]==miss_sign[2]) {printf("Error! The miss directions from %lf and %lf are both %i.\n",zeta[0],zeta[1],miss_sign[0]); return 1;}	// Did we make a bad choice for the ends of the interval? This should only be possible on the first loop
	for (i=0;i<bisect_limit;i++){
		printf("[%lf,%lf]\n",zeta[0],zeta[1]);
		init_pendulum(&x_0,y_0,0.5*(zeta[0]+zeta[1]));
		HIT = shoot(x_0,y_0,0.1,1000,NULL,&f_pendulum,&hit_pendulum);
		if(HIT.hit) {printf("Hit with parameter choice zeta = %lf\n",zeta[i]); break;}
		else if(HIT.miss_direction*(zeta[1]-zeta[0])>0) {zeta[1] = 0.5*(zeta[0]+zeta[1]);}
		else {zeta[0] = 0.5*(zeta[0]+zeta[1]);}
	}

	init_pendulum(&x_0,y_0,zeta[1]);
	shoot(x_0,y_0,0.1,1000,outfile,&f_pendulum,&hit_pendulum);
	if(HIT.hit==0) {printf("Failed to hit. Trajectory zeta = %lf has been written.\n",zeta[1]); return 1;}
	
	return 0;
}

struct hit_data shoot(double x_0, double *y_0, double delta_x, unsigned long x_steps,FILE *print, void (*derivative)(double x, double *y,double *deriv_ret), int (*hit_func)(double x, double *y)){
	// The shoot() function should perform a shot for given input parameters over the range of x specficied and conpare its shot against a 'hit' function
	struct hit_data HIT = {0,0,0.0};	// This gets returned
	unsigned long i,j;	// Dummy counter variables
	// Declare a variable to hold the state vector and allocate memory to it. Then initialise it.
	double *y = malloc(dim*sizeof(double));
	for (i=0;i<dim;i++) {y[i] = y_0[i];}
	double x = x_0;

	double *k_1 = malloc(dim*sizeof(double));
	double *k_2 = malloc(dim*sizeof(double));
	double *k_3 = malloc(dim*sizeof(double));
	double *k_4 = malloc(dim*sizeof(double));
	double *y_1 = malloc(dim*sizeof(double));

	if(print!=NULL){
		fprintf(print,"# x");
		for (i=0;i<dim;i++) {fprintf(print,"\ty_%lu",i);}
		fprintf(print,"\n");
		fprintf(print,"%lf",x);
		for (i=0;i<dim;i++) {fprintf(print,"\t%lf",y[i]);}
		fprintf(print,"\n");
	}

	for (j=0;j<x_steps;j++){
		/* Start of RK4 step */
		derivative(x,y,k_1);
		for (i=0;i<dim;i++) {y_1[i] = y[i] + 0.5*delta_x*k_1[i];}
		derivative(x + 0.5*delta_x,y_1,k_2);
		for (i=0;i<dim;i++) {y_1[i] = y[i] + 0.5*delta_x*k_2[i];}
		derivative(x + 0.5*delta_x,y_1,k_3);
		for (i=0;i<dim;i++) {y_1[i] = y[i] + delta_x*k_3[i];}
		derivative(x + delta_x,y_1,k_4);

		x = x + delta_x;
		for (i=0;i<dim;i++) {y[i] = y[i] + delta_x*(k_1[i] + 2*k_2[i] + 2*k_3[i] + k_4[i])/6;}
		/* End of RK4 step */
		
		if(print != NULL){
			fprintf(print,"%lf",x);
			for (i=0;i<dim;i++) {fprintf(print,"\t%lf",y[i]);}
			fprintf(print,"\n");
		}
		// Test whether we've hit our target yet
		if(hit_func(x,y)==0) {HIT.hit=1,HIT.miss_direction=0,HIT.hit_time=x; return HIT;};
	}

	// FIXME - do I need to free k_1,k_2,k_3,k_4,y_1,y at this point?
	
	// If I get to this point and haven't hit yet, need to know what direction I've missed by
	HIT.hit=0;
	HIT.miss_direction = hit_func(x,y);
	HIT.hit_time=0.0;
	return HIT;
}

/* Function definitions for derivatives */
void f(double x, double *y, double *deriv_ret){
	// Just some dummy placeholder function really
	deriv_ret[0] = 1.0;
	deriv_ret[1] = 1.0;
}

void f_pendulum(double x, double *y, double *deriv_ret){
	// The derivative appropriate to a pendumlum
	deriv_ret[0] = y[1];
	deriv_ret[1] = -sin(y[0]);
}
/* End of function definitions for derivatives */

/* Function definitions for initialisations */
void init_pendulum(double *x_0,double *y_0, double zeta){
	*x_0 = 0.0;
	y_0[0] = -pi + 0.3*cos(zeta*pi/2.0);
	y_0[1] = 0 + 0.3*sin(zeta*pi/2.0);
}
/* End of function definitions for initialisations */

/* Function declarations for hit testing */
int hit_pendulum(double x, double *y){
	float epsilon = 0.001;
	if((y[0]-pi)*(y[0]-pi) + y[1]*y[1] < epsilon*epsilon) {return 0;}
	else if(y[0] >= pi) {return +1;}
	else if(y[0] < pi) {return -1;}
	return 2;
}
/* Enf of function declarations for hit testing */
