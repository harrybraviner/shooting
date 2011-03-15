#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const unsigned short dim = 2;	// Dimension of the state vector
char *outfile_name = "output.dat";	// Some output file

const double pi = M_PI;
struct hit_data {int hit; double hit_time;};	// Data type to allow the return of data about hitting

struct hit_data shoot(double x_0, double *y_0, double delta_x, unsigned long x_steps, FILE *print, void (*derivative)(double x, double *y, double *deriv_ret));
/* Function declarations for derivatives */
void f(double x, double *y, double *deriv_ret);
void f_pendulum(double x, double *y, double *deriv_ret);
/* End of function declarations for derivatives */

int main(int argc, char *argv[]){

	FILE *outfile = fopen(outfile_name,"w");
	if(outfile==NULL) {fprintf(stderr, "Unable to open %s for writing. Quitting.\n",outfile_name); return -1;}

	// Example initial data for a pendulum
	double y_0[] = {-pi+0.1,0.0};
	double x_0 = 0.0;

	shoot(x_0,y_0,0.1,1000,outfile,&f_pendulum);	// Make a single shot and write the output
	
	return 0;
}

struct hit_data shoot(double x_0, double *y_0, double delta_x, unsigned long x_steps,FILE *print, void (*derivative)(double x, double *y,double *deriv_ret)){
	// The shoot() function should perform a shot for given input parameters over the range of x specficied and conpare its shot against a 'hit' function
	struct hit_data HIT = {0,0.0};	// This gets returned
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
		// FIXME - also insert calls to the 'hit' function
		// FIXME - the hit function is going to need some memory where it can store stuff between subsequent calls
	}

	// FIXME - do I need to free k_1,k_2,k_3,k_4,y_1,y at this point?

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
