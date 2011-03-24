/* SUGRA.c
FIXME - there's still stuff missing from here
Contains the functions	f_SUGRA (which gives the time derivative for this system, and hence defines the dynamical system)
			init_SUGRA_AdS 

The ordering of the fields in the state vector is \hat{\beta}, \partial_{\rho} \hat{\beta}, \zeta, \partial_{\rho} \zeta, \hat{F}, \partial_{\rho} \hat{F}, e^{-2\hat{h}}, \partial_{\rho} e^{-2\hat{h}}

*/

const double g2gamma2 = 1;	// g^2 \gamma^2   This turns out to be the only parameter we need to set.

const double zeta_init = 0.7;	// The value of \zeta of the AdS spacetime in the IR. This should be > (1 - 1/sqrt(6)) ~ 0.5918

const double SUGRA_epsilon = 0.001;	// The magnitude of the perturbation away from the fixed point

/* Start of function declarations */
void f_SUGRA(double x, double *y, double *deriv_ret);

/* Functions for AdS ---> AdS flows */
void init_AdS_large_zeta(double *x, double *y, double init_param);
/* End of functions for AdS ---> AdS flows */

/* End of function declarations */

/* Start of function definitions */
void f_SUGRA(double x, double *y, double *deriv_ret){
	double A = (-y[6] + (1 + 4*y[2] - y[2]*y[2] - (y[2] + 4*g2gamma2*y[6]*y[6])*y[0]*y[0] - 4*g2gamma2*y[2]*y[6]*y[6])/(4*sqrt(y[2])))/(1 + (y[7]/y[6])*(0.25*y[7]/y[6] - x*y[5] - y[4] - 2) + 2*(x*y[5] + y[4])*(x*y[5] + y[4]) - (y[3]/y[2])*(y[3]/y[2])/8 - (y[1] + 2*y[0])*(y[1] + 2*y[0])/(4*y[2]));	// This is e^{-2\hat{B}}
	double B = (1 + 4*y[2] - y[2]*y[2] - (3*y[2] + 4*g2gamma2*y[6]*y[6])*y[0]*y[0] + 4*g2gamma2*y[2]*y[6]*y[6])/(4*sqrt(y[2])) - 2*A*((y[1] + 2*y[0])*(y[1] + 2*y[0])/(4*y[2]) + 2 + x*y[5] + y[4] - y[7]/y[6]);	// This is \partial_{\rho} e^{-2\hat{B}}
	deriv_ret[0] = y[1];
	deriv_ret[1] = -2*y[1] + (y[1] + 2*y[0])*(y[7]/y[6] - x*y[5] - y[4] + y[3]/y[2] - 0.5*B/A) + sqrt(y[2])*(y[2] + 2*g2gamma2*y[6]*y[6])*y[0]/A;
	deriv_ret[2] = y[3];
	deriv_ret[3] = y[3]*y[3]/y[2] - y[3]*(2 + x*y[5] + y[4] + 0.5*B/A - y[7]/y[6]) + sqrt(y[2])*(1 - 4*y[2] + 3*y[2]*y[2] + (y[2] - 4*g2gamma2*y[6]*y[6])*y[0]*y[0])/(2*A) - (y[1] + 2*y[0])*(y[1] + 2*y[0]) + 2*g2gamma2*y[2]*sqrt(y[2])*y[6]*y[6]/A;
	deriv_ret[4] = y[4];
	deriv_ret[5] = -2*y[5]/x + x*(y[5] + y[4]/x)*(y[5] + y[4]/x) + (y[5] + y[4]/x)*(-2 + y[7]/y[6] - 0.5*B/A) + (1 + 4*y[2] -y[2]*y[2] + (y[2] + 12*g2gamma2*y[6]*y[6])*y[0]*y[0] + 4*y[2]*g2gamma2*y[6]*y[6])/(8*x*A*sqrt(y[2])) + (y[1] + 2*y[0])*(y[1] + 2*y[0])/(4*x*y[2]);
	deriv_ret[6] = y[7];
	deriv_ret[7] = -y[7]*(2 + x*y[5] + y[4] + 0.5*B/A - y[7]/y[6]) + 2*y[6]*y[6]/A - y[6]*(1 + 4*y[2] - y[2]*y[2] + (y[2] - 4*g2gamma2*y[6]*y[6])*y[0]*y[0])/(4*sqrt(y[2])*A) - y[6]*(y[1] + 2*y[0])*(y[1] + 2*y[0])/(2*y[2]) + 3*sqrt(y[2])*g2gamma2*y[6]*y[6]*y[6]/A;
}

void init_AdS_large_zeta(double *x, double *y, double init_param){
	// Set y[] to be the fixed point
	y[0] = 0;	// \hat{\beta}
	y[1] = 0;	// \partial_{\rho} \hat{\beta}
	y[2] = zeta_init;	// \zeta
	y[3] = 0;	// \partial_{\rho} \zeta
	y[4] = 1;	// \hat{F}
	y[5] = 1;	// \partial_{\rho} \hat{F}
	y[6] = (1 - 2*zeta_init + 2*zeta_init*zeta_init)/(2*sqrt(zeta_init));	// e^{-2\hat{h}}
	y[7] = 0;	// \partial_{\rho} e^{-2\hat{h}}
	// Declare the the two unstable eigenvectors
	double yp_5[8],yp_6[8];
	// The following is the eigenvector which is always unstable
	yp_5[0] = yp_5[1] = yp_5[4] = yp_5[5] = 0;
	yp_5[6] = 1;	// e^{-2\hat{h}}
	yp_5[2] = yp_5[6]*(4*zeta_init*sqrt(zeta_init)*(4*zeta_init - sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1)))/(1-2*zeta_init+2*zeta_init*zeta_init)*(1-3*zeta_init);	// \zeta
	yp_5[3]	= (-1.5 + 0.5*sqrt(3)*sqrt((-7*zeta_init*zeta_init + 22*zeta_init - 4 + 4*(1-zeta_init)*sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1))/(zeta_init*(2-zeta_init))))*yp_5[2];	// \partial_{\rho} \zeta
	yp_5[7]	= (-1.5 + 0.5*sqrt(3)*sqrt((-7*zeta_init*zeta_init + 22*zeta_init - 4 + 4*(1-zeta_init)*sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1))/(zeta_init*(2-zeta_init))))*yp_5[6];	// \partial_{\rho} e^{-2\hat{h}}
	// The following eigenvector is unstable for zeta_init > (1 - 1/sqrt(6)) ~ 0.5918
	yp_6[0] = yp_6[1] = yp_6[4] = yp_6[5] = 0;
	yp_6[2] = 1;	// \zeta
	yp_6[6] = yp_6[2]*(1-2*zeta_init+2*zeta_init*zeta_init)*(1-3*zeta_init)/(4*zeta_init*sqrt(zeta_init)*(4*zeta_init + sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1)));	// e^{\hat{h}}
	yp_6[3]	= (-1.5 + 0.5*sqrt(3)*sqrt((-7*zeta_init*zeta_init + 22*zeta_init - 4 - 4*(1-zeta_init)*sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1))/(zeta_init*(2-zeta_init))))*yp_6[2];	// \partial_{\rho} \zeta
	yp_6[7]	= (-1.5 + 0.5*sqrt(3)*sqrt((-7*zeta_init*zeta_init + 22*zeta_init - 4 - 4*(1-zeta_init)*sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1))/(zeta_init*(2-zeta_init))))*yp_6[6];	// \partial_{\rho} e^{-2\hat{h}}
	// Normalise these
	double yp_5_norm = 0, yp_6_norm = 0;
	short i;
	for(i=0;i<8;i++) {yp_5_norm += yp_5[i]*yp_5[i]; yp_6_norm += yp_6[i]*yp_6[i];}
	yp_5_norm = sqrt(yp_5_norm); yp_6_norm = sqrt(yp_6_norm);
	for(i=0;i<8;i++) {yp_5[i] = yp_5[i]/yp_5_norm; yp_6[i] = yp_6[i]/yp_6_norm;}
	// Now perturb the vector
	double pi = M_PI;
	for(i=0;i<8;i++) {y[i] += yp_5[i]*cos(init_param*pi) + yp_6[i]*sin(init_param*pi);}
}
// FIXME - write an initialisation function for the Lifshitz spacetime

// FIXME - write hitting functions
/* End of function definitions */
