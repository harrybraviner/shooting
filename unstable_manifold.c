#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const unsigned short dim = 8;
const unsigned long step_limit = 10000;
const double delta_x = 0.001;
char *outfile_name = "unstable.dat";

#include "./shoot.c"
#include "./SUGRA.c"

int main(){
	FILE *outfile = fopen(outfile_name,"w");
	if(outfile==NULL) {fprintf(stderr,"Can't open %s for writing. Quitting.\n",outfile_name); return -1;}

	double x_0,y_0[dim];
	struct hit_data HIT;
	
	init_AdS_large_zeta(&x_0,y_0,0.0);
	HIT = shoot(x_0,y_0,delta_x,step_limit,outfile,f_SUGRA,hit_AdS_small_zeta);
	printf("HIT.miss_direction = %d\n",HIT.miss_direction);

	return 0;
}
