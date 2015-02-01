#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int dim = 8;

#include "shoot.c"
#include "SUGRA.c"

int main(){
	int i = 0;
	double step = 0.001;
	double x_0,y_0[dim];
	for (i=0; i <= 2/step; i++){
		init_Lifshitz(&x_0,y_0,step*i);
		if(shoot(x_0,y_0,0.01,1000,NULL,&f_SUGRA,&hit_AdS_small_zeta).hit==1){fprintf(stderr,"Missed by +1 with zeta = %lf\n",step*i);}
	}
	return 0;
}
