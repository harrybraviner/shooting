/* pendulum.c

Contains the functions	f_pendulum (which gives the time derivative for this system, and hence defines the dynamical system)
			init_pendulum (which defines a 1-parameter family of initial conditions which we believe contains a point on the trajectory that passed through the 'target' point)
			hit_pendulum (which assesses whether or not we have 'hit', ie. come sufficiently close to, the target point, which in this case in (0,+pi))
*/

const double pi = M_PI;

/* Start of function declarations */
void f_pendulum(double x, double *y, double *deriv_ret);
void init_pendulum(double *x, double *y, double zeta);
int hit_pendulum(double x, double *y);
/* End of function declarations */

/* Start of function definitions */
void f_pendulum(double x, double *y, double *deriv_ret){
	// The derivative appropriate to a pendumlum
	deriv_ret[0] = y[1];
	deriv_ret[1] = -sin(y[0]);
}

void init_pendulum(double *x_0,double *y_0, double zeta){
	*x_0 = 0.0;
	y_0[0] = -pi + 0.3*cos(zeta*pi/2.0);
	y_0[1] = 0 + 0.3*sin(zeta*pi/2.0);
}

int hit_pendulum(double x, double *y){
	float epsilon = 0.001;
	if((y[0]-pi)*(y[0]-pi) + y[1]*y[1] < epsilon*epsilon) {return 0;}
	else if(y[0] >= pi) {return +1;}
	else if(y[0] < pi) {return -1;}
	return 2;
}
/* End of function definitions */
