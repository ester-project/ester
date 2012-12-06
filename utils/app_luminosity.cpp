#include"star.h"
#include<stdlib.h>

int main(int argc,char *argv[]) {

	if(argc!=3) {
		printf("Usage: app_luminosity model_file i\n\t i: Inclination angle (deg.)\n");
		return 1;
	}

	star2d A;
	A.read(argv[1]);
	
	double i=atof(argv[2]);
	
	printf("L= %.3f Lsun\n", A.luminosity()/L_SUN);
	printf("L_app (%.2fº)= %.3f Lsun\n",i, A.apparent_luminosity(i/180.*PI)/L_SUN);
	return 0;
	

}
