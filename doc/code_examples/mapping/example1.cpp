/* The following example illustrates the use of the mapping
library in spheroidal coordinates */

#include<stdio.h>
#include"mapping.h"
#include"constants.h" //For the defintion of PI

//Function prototypes
double laplacian(matrix y,double r0,double th0);
double integral(matrix y);

// Define mapping objects as global variable
mapping map;
// Create references for spherical coordinates
matrix &r=map.r,&th=map.th;

int main() {

	//Initialize map. In the example we will use 2 domains
	map.set_ndomains(2);
	map.set_npts(100);
	map.set_nt(50);
	//map.set_nex(20); // We won't use the external domain
	map.init();
	
	map.R.setrow(0,1e-3*ones(1,map.nt)); 
							// Use 1e-3 as the interior limit (instead of 0)
							// to avoid a division by zero in the 
							// calculation of the laplacian
	map.R.setrow(1,0.5+0.1*sin(th)*sin(th));
	map.R.setrow(2,ones(1,map.nt));
	map.remap();
	
	matrix y;
	
	//Define the function y
	y=r*r*r*(1+sin(th)*sin(th));
	
	double lap_y,int_y;
	double r0=0.3,th0=PI/3;
	
	lap_y=laplacian(y,r0,th0);
	int_y=integral(y);
	
	printf("The value of the laplacian at (%f,%f) is %e\n",r0,th0,lap_y);
	printf("The volume integral is %e\n",int_y);
	
	return 0;	
	
}

// Function for calculating the laplacian of y at (r0,th0)
double laplacian(matrix y,double r0,double th0) {

	matrix lap_y;
	
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt;
	matrix &rz=map.rz,&rzz=map.rzz,&rzt=map.rzt,&rt=map.rt,&rtt=map.rtt;
	matrix &Dt=map.Dt,&Dt2=map.Dt2;
	matrix_block_diag &D=map.D;
	
	lap_y=gzz*(D,D,y)+2*gzt*(D,y,Dt)+(y,Dt2)/r/r
		+(2./r/rz-rtt/r/r/rz-gzz*rzz/rz-gzt*(2*rzt/rz-cos(th)/sin(th)))*(D,y)
		+cos(th)/r/r/sin(th)*(y,Dt);
	
	lap_y=map.eval(lap_y,r0*ones(1,1),th0*ones(1,1));
	
	return lap_y(0);

}
// Function for calculating the volume integral of y
double integral(matrix y) {

	return 2*PI*(map.I,r*r*map.rz*y,map.It)(0);

}

