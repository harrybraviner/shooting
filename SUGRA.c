/* SUGRA.c
FIXME - there's still stuff missing from here
Contains the functions	f_SUGRA (which gives the time derivative for this system, and hence defines the dynamical system)
			init_SUGRA_AdS 

The ordering of the fields in the state vector is \hat{\beta}, \partial_{\rho} \hat{\beta}, \zeta, \partial_{\rho} \zeta, \hat{F}, \partial_{\rho} \hat{F}, e^{-2\hat{h}}, \partial_{\rho} e^{-2\hat{h}}

*/

const double g2gamma2 = 1.1469463;	// g^2 \gamma^2   This turns out to be the only parameter we need to set. IT IS RELATED TO ZETA_INIT

const double zeta_init = 0.55;	// The value of \zeta of the AdS spacetime in the IR. This should be > (1 - 1/sqrt(6)) ~ 0.5918

const double SUGRA_init_e = 0.001;	// The magnitude of the perturbation away from the fixed point
const double SUGRA_hit_e = 0.01;	// The maximum deviation of any of the fields from the AdS fixed point that we are prepared to count as a hit

/* Start of function declarations */
void f_SUGRA(double x, double *y, double *deriv_ret);
double other_zeta(double zeta);

/* Functions for AdS ---> AdS flows */
void init_AdS_large_zeta(double *x, double *y, double init_param);
int hit_AdS_small_zeta(double x, double *y);
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
	deriv_ret[4] = y[5];
	deriv_ret[5] = -2*y[5]/x - x*(y[5] + y[4]/x)*(y[5] + y[4]/x) + (y[5] + y[4]/x)*(-2 + y[7]/y[6] - 0.5*B/A) + (1 + 4*y[2] -y[2]*y[2] + (y[2] + 12*g2gamma2*y[6]*y[6])*y[0]*y[0] + 4*y[2]*g2gamma2*y[6]*y[6])/(8*x*A*sqrt(y[2])) + (y[1] + 2*y[0])*(y[1] + 2*y[0])/(4*x*y[2]);
	deriv_ret[6] = y[7];
	deriv_ret[7] = -y[7]*(2 + x*y[5] + y[4] + 0.5*B/A - 2*y[7]/y[6]) + 2*y[6]*y[6]/A - y[6]*(1 + 4*y[2] - y[2]*y[2] + (y[2] - 4*g2gamma2*y[6]*y[6])*y[0]*y[0])/(4*sqrt(y[2])*A) - y[6]*(y[1] + 2*y[0])*(y[1] + 2*y[0])/(2*y[2]) + 3*sqrt(y[2])*g2gamma2*y[6]*y[6]*y[6]/A;
}

double other_zeta(double x){
	// There are typically two allowed values of \zeta for each value of g2gamma2. This function calculates the other one
	double other_zeta = cbrt((2*x*x - 2*x + 1)*sqrt((-48*pow(x,6) + 96*pow(x,5) - 76*pow(x,4) - 32*x*x*x + 128*x*x - 120*x +25)/(3*x*x - 4*x + 1))/(8*3*sqrt(3)*(3*x*x - 4*x + 1)) - (12*pow(x,4) - 4*x*x*x - 12*x*x + 18*x - 13)/(324*x - 108)) + (12*pow(x,4) - 16*x*x*x + 8*x*x + 4*x + 1)/((108*x*x - 144*x +36)*cbrt((2*x*x - 2*x + 1)*sqrt((-48*pow(x,6) + 96*pow(x,5) - 76*pow(x,4) - 32*x*x*x + 128*x*x - 120*x +25)/(3*x*x - 4*x + 1))/(8*3*sqrt(3)*(3*x*x - 4*x + 1)) - (12*pow(x,4) - 4*x*x*x - 12*x*x + 18*x - 13)/(324*x - 108))) - (x-2)/3;
	return other_zeta;
}

void init_AdS_large_zeta(double *x, double *y, double init_param){
	*x = 1.0;
	// Set y[] to be the fixed point
	y[0] = 0;	// \hat{\beta}
	y[1] = 0;	// \partial_{\rho} \hat{\beta}
	y[2] = zeta_init;	// \zeta
	y[3] = 0;	// \partial_{\rho} \zeta
	y[4] = 1;	// \hat{F}
	y[5] = 0;	// \partial_{\rho} \hat{F}
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
	for(i=0;i<8;i++) {y[i] += SUGRA_init_e*yp_5[i]*cos(init_param*pi) + SUGRA_init_e*yp_6[i]*sin(init_param*pi);}
}
// FIXME - write an initialisation function for the Lifshitz spacetime

int hit_AdS_small_zeta(double x, double *y){
	double zeta = other_zeta(zeta_init); // The (small) \zeta value of the spacetime that we're aiming for
	if((y[0] - 0)*(y[0] - 0)	// \hat{\beta}
	+ (y[1] - 0)*(y[1] - 0)		// \partial_{\rho} \hat{\beta}
	+ (y[2] - zeta)*(y[2] - zeta)	// \zeta
	+ (y[3] - 0)*(y[3] - 0)		// \partial_{\rho} \zeta
	+ (y[4] - 1)*(y[4] - 1)		// \hat{F}
	+ (y[5] - 0)*(y[5] - 0)		// \partial_{\rho} \hat{F}
	+ (y[6] - (1 - 2*zeta + 2*zeta*zeta)/(2*sqrt(zeta)))*(y[6] - (1 - 2*zeta + 2*zeta*zeta)/(2*sqrt(zeta)))	// e^{-2\hat{h}}
	+ (y[7] - 0)*(y[7] - 0)		// \partial_{\rho} e^{-2\hat{h}}
	< SUGRA_hit_e*SUGRA_hit_e){return 0;}
	else if(y[2] - zeta >= 0){return +1;}
	else {return -1;} // Our decision is purely based on whether \zeta is too large or too small
}
/* End of function definitions */
