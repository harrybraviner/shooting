struct hit_data {int hit; int miss_direction; double hit_time; unsigned long hit_steps;};	// Data type to allow the return of data about hitting

struct hit_data shoot(double x_0, double *y_0, double delta_x, unsigned long x_steps, FILE *print, void (*derivative)(double x, double *y, double *deriv_ret), int (*hit_func)(double x, double *y));

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
		
		/* If we have been passed a file stream, write the trajectory to it */
		if(print != NULL){
			fprintf(print,"%lf",x);
			for (i=0;i<dim;i++) {fprintf(print,"\t%lf",y[i]);}
			fprintf(print,"\n");
		}
		/* End of writing */

		// Test whether we've hit our target yet
		if(hit_func(x,y)==0) {HIT.hit=1,HIT.miss_direction=0,HIT.hit_time=x,HIT.hit_steps=j+1; return HIT;};
	}
	free(k_1),free(k_2),free(k_3),free(k_4),free(y_1);
	
	// If I get to this point and haven't hit yet, need to know what direction I've missed by
	HIT.hit=0;
	HIT.miss_direction = hit_func(x,y);
	HIT.hit_time=0.0;
	HIT.hit_steps=0;
	return HIT;
}

