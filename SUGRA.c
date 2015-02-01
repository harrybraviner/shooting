/* SUGRA.c
The ordering of the fields in the state vector is \hat{\beta}, \partial_{\rho} \hat{\beta}, \zeta, \partial_{\rho} \zeta, \hat{F}, \partial_{\rho} \hat{F}, e^{-2\hat{h}}, \partial_{\rho} e^{-2\hat{h}}

Setup instructions:
	Set the function pointers *derivative, *hit_func, *init_func in shooting.c
	Set the interval zeta[] in shooting.c
	Set SUGRA_init_e and SUGRA_hit_e in this file
For an AdS ---> AdS flow:
	Set zeta_init ( >1-1/sqrt(6), <1)
	Calculate and set g2gamma2 from zeta_init
For an AdS ---> (lower)Li flow:
	Set z_init ( >4.29, <7.01)
	Calculate and set g2gamma2 from z_init
	Calculate and set yp_3[] and yp_5[] in init_ls_Lifshitz()
For Li ---> AdS flow:
	Set zeta_init
	Calculate and set g2gamma2 from zeta_init
For Li --> Li flow:
	Set z_init (>7.01)
	Calculate and set g2gamma2 and z_target from z_init
	Calculate and set yp_3[] and yp_5[] in init_ls_Lifshitz()
For Li --> AdS --> AdS flow:
	NB: no idea if this will work!
	Set zeta_init (>0.592)
	Calculate and set g2gamma2 from zeta_init
*/

const double g2gamma2 = 2.049711539142868;	// g^2 \gamma^2   This turns out to be the only physical parameter we need to set.

const double zeta_init = 0.8;	// The value of \zeta of the AdS spacetime set up in the init_AdS_large_zetaxx functions. 
				// This should be > (1 - 1/sqrt(6)) ~ 0.5918 and < 1 for init_AdS_large_zeta56
				// This should be > 1 - sqrt(2/5) ~ 0.3675 for init_AdS_large_zeta35

const double z_init = 10.0;	// The value of z of the (lower sign) Lifshitz spacetime set up in init_ls_Lifshitz
const double z_target = 2.91000841281783;	// The value of z corresponding the (upper sign) Lifshitz spacetime that we expect to hit

const double SUGRA_init_e = 0.01;	// The magnitude of the perturbation away from the fixed point
const double SUGRA_init_e2 = 0.0002;	// The magnitude of the \hat{\beta} perturbation used in init_AdS_large_zeta56B()
const double SUGRA_hit_e = 0.0001;	// The maximum deviation of any of the fields from the AdS fixed point that we are prepared to count as a hit

/* Start of function declarations */
void f_SUGRA(double x, double *y, double *deriv_ret);
double other_zeta(double zeta);
double small_zeta_from_g2gamma2(double g2gamma2);

// Functions for AdS ---> AdS flows
void init_AdS_large_zeta56(double *x, double *y, double init_param);
int hit_AdS_small_zeta(double x, double *y);

// Functions for AdS ---> Li flows
void init_ls_Lifshitz(double *x, double *y, double init_param);
int hit_PFP_or_singularity(double x, double *y);

// Functions for Li ---> AdS flows
void init_AdS_large_zeta35(double *x,double *y, double init_param);
// Use hit_PFP_or_singularity() for these shots

// Functions for Li ---> Li flows
// Use init_ls_Lifshitz()
int hit_us_Li(double x, double *y);

// Functions for Li --> AdS --> AdS flows
void init_AdS_large_zeta56B(double *x, double *y, double init_param);
// Use hit_PFP_or_singularity() for these shots

// Function for shooting from upper sign Lifshitz to the 6D spacetime
void init_us_Lifshitz3(double *x, double *y, double init_param);

/* End of function declarations */

/* Start of function definitions */
void f_SUGRA(double x, double *y, double *deriv_ret){
	double A = (-y[6] + (1 + 4*y[2] - y[2]*y[2] - (y[2] + 4*g2gamma2*y[6]*y[6])*y[0]*y[0] - 4*g2gamma2*y[2]*y[6]*y[6])/(4*sqrt(y[2])))/(1 + (y[7]/y[6])*(0.25*y[7]/y[6] - x*y[5] - y[4] - 2) + 2*(x*y[5] + y[4]) - (y[3]/y[2])*(y[3]/y[2])/8 - (y[1] + 2*y[0])*(y[1] + 2*y[0])/(4*y[2]));	// This is e^{-2\hat{B}}
	double B = (1 + 4*y[2] - y[2]*y[2] - (3*y[2] + 4*g2gamma2*y[6]*y[6])*y[0]*y[0] + 4*g2gamma2*y[2]*y[6]*y[6])/(4*sqrt(y[2])) - 2*A*((y[1] + 2*y[0])*(y[1] + 2*y[0])/(4*y[2]) + 2 + x*y[5] + y[4] - y[7]/y[6]);	// This is \partial_{\rho} e^{-2\hat{B}}
	deriv_ret[0] = y[1];
	deriv_ret[1] = -2*y[1] + (y[1] + 2*y[0])*(y[7]/y[6] - x*y[5] - y[4] + y[3]/y[2] - 0.5*B/A) + sqrt(y[2])*(y[2] + 4*g2gamma2*y[6]*y[6])*y[0]/A;
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

double large_zeta_from_g2gamma2(double x){
	// Based on the constant g^2 \gamma^2, this function returns the larger allowed value of \zeta for an AdS spacetime
	// This function is no currently used by anything
	double C = pow(36*x*sqrt((27+99*x-8*x*x-80*x*x*x)/x) - 64*sqrt(3)*x*x*x + 8*pow(3,3.5)*x + pow(3,3.5),1/3.0); // Helper constant

	double large_zeta = (sqrt(-((12*x*C*C + (32*pow(3,7.0/6.0)*x*x + 16*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)*sqrt((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C) - 32*pow(3,13.0/4.0)*x*x*C)/C) + pow(((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C),0.75) + 4*pow(3,13.0/12.0)*x*pow(((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C),0.25))/(8*pow(3,13.0/12.0)*x*pow(((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C),0.25));
	return large_zeta;
}

double small_zeta_from_g2gamma2(double x){
	// Based on the constant g^2 \gamma^2, this function returns the smaller allowed value of \zeta for an AdS spacetime
	double C = pow(36*x*sqrt((27+99*x-8*x*x-80*x*x*x)/x) - 64*sqrt(3)*x*x*x + 8*pow(3,3.5)*x + pow(3,3.5),1/3.0); // Helper constant
	double small_zeta = -(sqrt(-((12*x*C*C + (32*pow(3,7.0/6.0)*x*x + 16*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)*sqrt((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C) - 32*pow(3,13.0/4.0)*x*x*C)/C) - pow(((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C),0.75) - 4*pow(3,13.0/12.0)*x*pow(((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C),0.25))/(8*pow(3,13.0/12.0)*x*pow(((12*x*C*C - (16*pow(3,7.0/6.0)*x*x + 8*pow(3,13.0/6.0)*x)*C + 64*pow(3,4.0/3.0)*x*x*x + 4*pow(3,10.0/3.0)*x)/C),0.25));
	return small_zeta;
}

void init_AdS_large_zeta56(double *x, double *y, double init_param){
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

void init_AdS_large_zeta56B(double *x, double *y, double init_param){
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
	double yp_5[8],yp_6[8],yp_3[8];
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
	// Also need some perturbation in the \hat{\beta} direction
	yp_3[0] = 1;	// \hat{\beta}
	yp_3[1] = -1.5 + 0.5*sqrt(1 - 24*(2*zeta_init*zeta_init - 4*zeta_init + 1)/(zeta_init*(2-zeta_init)));	// \partial_{\rho} \hat{\beta}
	yp_3[2] = yp_3[3] = yp_3[4] = yp_3[5] = yp_3[6] = yp_3[7] = 0;
	// Normalise these
	double yp_5_norm = 0, yp_6_norm = 0, yp_3_norm = 0;
	short i;
	for(i=0;i<8;i++) {yp_5_norm += yp_5[i]*yp_5[i]; yp_6_norm += yp_6[i]*yp_6[i]; yp_3_norm += yp_3[i]*yp_3[i];}
	yp_5_norm = sqrt(yp_5_norm); yp_6_norm = sqrt(yp_6_norm); yp_3_norm = sqrt(yp_3_norm);
	for(i=0;i<8;i++) {yp_5[i] = yp_5[i]/yp_5_norm; yp_6[i] = yp_6[i]/yp_6_norm; yp_3[i] = yp_3[i]/yp_3_norm;}
	// Now perturb the vector
	double pi = M_PI;
	for(i=0;i<8;i++) {y[i] += SUGRA_init_e*yp_5[i]*cos(init_param*pi) + SUGRA_init_e*yp_6[i]*sin(init_param*pi) + SUGRA_init_e2*yp_3[i];}
}

void init_AdS_large_zeta35(double *x, double *y, double init_param){
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
	double yp_5[8],yp_3[8];
	// The following is the eigenvector which is always unstable
	yp_5[0] = yp_5[1] = yp_5[4] = yp_5[5] = 0;
	yp_5[6] = 1;	// e^{-2\hat{h}}
	yp_5[2] = yp_5[6]*(4*zeta_init*sqrt(zeta_init)*(4*zeta_init - sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1)))/(1-2*zeta_init+2*zeta_init*zeta_init)*(1-3*zeta_init);	// \zeta
	yp_5[3]	= (-1.5 + 0.5*sqrt(3)*sqrt((-7*zeta_init*zeta_init + 22*zeta_init - 4 + 4*(1-zeta_init)*sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1))/(zeta_init*(2-zeta_init))))*yp_5[2];	// \partial_{\rho} \zeta
	yp_5[7]	= (-1.5 + 0.5*sqrt(3)*sqrt((-7*zeta_init*zeta_init + 22*zeta_init - 4 + 4*(1-zeta_init)*sqrt(25*zeta_init*zeta_init - 6*zeta_init + 1))/(zeta_init*(2-zeta_init))))*yp_5[6];	// \partial_{\rho} e^{-2\hat{h}}
	// The following eigenvector is unstable for zeta_init > (1 - sqrt(2/5)) ~ 0.3675
	yp_3[0] = 1;	// \hat{\beta}
	yp_3[1] = -1.5 + 0.5*sqrt(1 - 24*(2*zeta_init*zeta_init - 4*zeta_init + 1)/(zeta_init*(2-zeta_init)));	// \partial_{\rho} \hat{\beta}
	yp_3[2] = yp_3[3] = yp_3[4] = yp_3[5] = yp_3[6] = yp_3[7] = 0;
	// Normalise these
	double yp_5_norm = 0, yp_3_norm = 0;
	short i;
	for(i=0;i<8;i++) {yp_5_norm += yp_5[i]*yp_5[i]; yp_3_norm += yp_3[i]*yp_3[i];}
	yp_5_norm = sqrt(yp_5_norm); yp_3_norm = sqrt(yp_3_norm);
	for(i=0;i<8;i++) {yp_5[i] = yp_5[i]/yp_5_norm; yp_3[i] = yp_3[i]/yp_3_norm;}
	// Now perturb the vector
	double pi = M_PI;
	for(i=0;i<8;i++) {y[i] += SUGRA_init_e*yp_5[i]*cos(init_param*pi) + SUGRA_init_e*yp_3[i]*sin(init_param*pi);}
}

void init_ls_Lifshitz(double *x, double *y, double init_param){
	*x = 1.0;
	// Set y[] to be the Li fixed point for the "lower" s sign choice
	y[0] = pow((6 + z_init + 2*sqrt(2*(z_init + 4)))/(z_init*z_init*(z_init + 4)),0.25)*sqrt(z_init - 1);	// \hat{\beta}
	y[1] = 0;	// \partial_{\rho} \hat{\beta}
	y[2] = sqrt((6 + z_init + 2*sqrt(2*(z_init + 4)))/(z_init*z_init*(z_init + 4)));	// \zeta
	y[3] = 0;	// \partial_{\rho} \zeta
	y[4] = z_init;	// \hat{F}
	y[5] = 0;	// \partial_{\rho} \hat{F}
	y[6] = (6 + 3*z_init + 2*sqrt(2*(z_init + 4)))/(2*sqrt(z_init)*pow((z_init + 4),0.75)*pow(6 + z_init + 2*sqrt(2*(z_init + 4)),0.25));	// e^{-2\hat{h}}
	y[7] = 0;	// \partial_{\rho} e^{-2\hat{h}}

	/* PERTURBATION VECTORS */
	// I do not currently have any code in here to generate the vectors along the unstable directions. They must be entered by hand
	// They are currently set to the values appropriate for z =  
	double yp_3[] = {0.0029352971374196,0.068625551700454,0.0017479543923451,0.040866164107448,0.04254255230109,0.99462032842194,-0.0015627302704824,-0.036535731120268};
	double yp_5[] = {0.019468982821528,0.066306325488752,0.0012027012380506,0.0040960896872197,0.28052408380601,0.95539252064666,0.0055840312707358,0.019017767819249};
	/* ------------------- */
	// Switch from 'CAS' variables to 'new' variables
	yp_3[2] = -2*sqrt(2)*yp_3[2]; yp_5[2] = -2*sqrt(2)*yp_5[2];
	yp_3[3] = -2*sqrt(2)*yp_3[3]; yp_5[3] = -2*sqrt(2)*yp_5[3];
	yp_3[6] = -2*yp_3[6]; yp_5[6] = -2*yp_5[6];
	yp_3[7] = -2*yp_3[7]; yp_5[7] = -2*yp_5[7];

	// Normalise these
	double yp_3_norm = 0, yp_5_norm = 0;
	short i;
	for(i=0;i<8;i++) {yp_3_norm += yp_3[i]*yp_3[i]; yp_5_norm += yp_5[i]*yp_5[i];}
	yp_3_norm = sqrt(yp_3_norm); yp_5_norm = sqrt(yp_5_norm);
	for(i=0;i<8;i++) {yp_3[i] = yp_3[i]/yp_3_norm; yp_5[i] = yp_5[i]/yp_5_norm;}
	// Now perturb the vector
	double pi = M_PI;
	for(i=0;i<8;i++) {y[i] += SUGRA_init_e*yp_3[i]*cos(init_param*pi) + SUGRA_init_e*yp_5[i]*sin(init_param*pi);}
}

void init_us_Lifshitz3(double *x, double *y, double init_param){
	*x = 1.0;
	// Set y[] to be the Li fixed point for the "upper" s sign choice
	y[0] = pow((6 + z_init - 2*sqrt(2*(z_init + 4)))/(z_init*z_init*(z_init + 4)),0.25)*sqrt(z_init - 1);	// \hat{\beta}
	y[1] = 0;	// \partial_{\rho} \hat{\beta}
	y[2] = sqrt((6 + z_init - 2*sqrt(2*(z_init + 4)))/(z_init*z_init*(z_init + 4)));	// \zeta
	y[3] = 0;	// \partial_{\rho} \zeta
	y[4] = z_init;	// \hat{F}
	y[5] = 0;	// \partial_{\rho} \hat{F}
	y[6] = (6 + 3*z_init - 2*sqrt(2*(z_init + 4)))/(2*sqrt(z_init)*pow((z_init + 4),0.75)*pow(6 + z_init - 2*sqrt(2*(z_init + 4)),0.25));	// e^{-2\hat{h}}
	y[7] = 0;	// \partial_{\rho} e^{-2\hat{h}}

	/* PERTURBATION VECTORS */
	// I do not currently have any code in here to generate the vectors along the unstable directions. They must be entered by hand
	// They are currently set to the values appropriate for z =  
	double yp_3[] = {-0.0096841534977808,-0.084248385406061,-0.0094678909365473,-0.082366984867348,-0.11272679075647,-0.98067941379207,0,0};
	double yp_5[] = {0,0,0,0,0,0,0,0};
	/* ------------------- */
	// Switch from 'CAS' variables to 'new' variables
	yp_3[2] = -2*sqrt(2)*yp_3[2]; yp_5[2] = -2*sqrt(2)*yp_5[2];
	yp_3[3] = -2*sqrt(2)*yp_3[3]; yp_5[3] = -2*sqrt(2)*yp_5[3];
	yp_3[6] = -2*yp_3[6]; yp_5[6] = -2*yp_5[6];
	yp_3[7] = -2*yp_3[7]; yp_5[7] = -2*yp_5[7];

	// Normalise these
	double yp_3_norm = 0, yp_5_norm = 0;
	short i;
	for(i=0;i<8;i++) {yp_3_norm += yp_3[i]*yp_3[i]; yp_5_norm += yp_5[i]*yp_5[i];}
	yp_3_norm = sqrt(yp_3_norm); yp_5_norm = sqrt(yp_5_norm);
	if(yp_3_norm==0) {yp_3_norm = 1;}
	if(yp_5_norm==0) {yp_5_norm = 1;}
	for(i=0;i<8;i++) {yp_3[i] = yp_3[i]/yp_3_norm; yp_5[i] = yp_5[i]/yp_5_norm;}
	// Now perturb the vector
	double pi = M_PI;
	for(i=0;i<8;i++) {y[i] += SUGRA_init_e*yp_3[i]*cos(init_param*pi) + SUGRA_init_e*yp_5[i]*sin(init_param*pi);}
}

int hit_AdS_small_zeta(double x, double *y){
	// double zeta = other_zeta(zeta_init); // The (small) \zeta value of the spacetime that we're aiming for
	double zeta = small_zeta_from_g2gamma2(g2gamma2);	// The \zeta value of the spacetime that we're aiming for
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

int hit_PFP_or_singularity(double x, double *y){
	// Decides whether we've hit the "pseudo-fixed point" at \zeta=1/3, e^{-2\hat{h}}=\hat{\beta} = 0
	// Returns -1 if we have, otherwise +1
	if((y[2]-1.0/3.0 < SUGRA_hit_e)	// Is \zeta -1/3 small enough?
	& (y[0] < SUGRA_hit_e)		// Is \hat{\beta} small enough?
	& (y[6] < SUGRA_hit_e)		// Is e^{-2\hat{h}} small enough?
	){return -1;}
	else {return +1;}
}

int hit_us_Li(double x, double *y){
	double zeta = sqrt((6+z_target-2*sqrt(2*z_target+8))/(4+z_target))/z_target;	// The \zeta value of the spacetime that we're aiming for
	if((y[0] - 0)*(y[0] - 0)	// \hat{\beta}
	+ (y[1] - 0)*(y[1] - 0)		// \partial_{\rho} \hat{\beta}
	+ (y[2] - zeta)*(y[2] - zeta)	// \zeta
	+ (y[3] - 0)*(y[3] - 0)		// \partial_{\rho} \zeta
	+ (y[4] - 1)*(y[4] - 1)		// \hat{F}
	+ (y[5] - 0)*(y[5] - 0)		// \partial_{\rho} \hat{F}
	+ (y[6] - (1 - 2*zeta + 2*zeta*zeta)/(2*sqrt(zeta)))*(y[6] - (1 - 2*zeta + 2*zeta*zeta)/(2*sqrt(zeta)))	// e^{-2\hat{h}}
	+ (y[7] - 0)*(y[7] - 0)		// \partial_{\rho} e^{-2\hat{h}}
	< SUGRA_hit_e*SUGRA_hit_e){return 0;}
	// else if(y[2] - zeta >= 0){return +1;}
	// else {return -1;} // Our decision is purely based on whether \zeta is too large or too small
	else if(y[0] > pow(((6+z_target-2*sqrt(2*z_target+8))/(z_target*z_target*(4+z_target))),0.25)*sqrt(z_target-1)){return +1;}
	else {return -1;}	// Our decision is based purely on whether \beta is too large or too small
}

/* End of function definitions */
