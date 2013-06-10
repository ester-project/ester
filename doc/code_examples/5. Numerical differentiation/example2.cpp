/* The following example illustrates the use of the numerical differentiation
library in 2D in spherical coordinates */

#include<stdio.h>
#include"numdiff.h"
#include"constants.h" //For the defintion of PI

//Function prototypes
double laplacian(matrix y,double r0,double th0);
double integral(matrix y);

// Define diff_gl and diff_leg objects as global variables
diff_gl gl; 
diff_leg leg;
// Create references for spherical coordinates
matrix &r=gl.x,&th=leg.th;

int main() {

	//Initialize gl. In the example we will use 2 domains
	gl.set_ndomains(2);
	gl.set_xif(1e-12,0.2,1.); // Use 1e-12 as the interior limit (instead of 0)
							  //to avoid the central singularity
	gl.set_npts(100,100);
	gl.init();
	//Initialize leg
	leg.npts=50;
	leg.init();
	
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
	
	lap_y=(gl.D,r*r*gl.D,y)/r/r+(y,leg.lap_00)/r/r;
	
	//Interpolate in the direction of r
	lap_y=gl.eval(lap_y,r0); // Now lap_y is (1,nth)
	//Interpolate in the direction of theta
	lap_y=leg.eval_00(lap_y,th0); // Now lap_y is (1,1)
	
	return lap_y(0);
	// lap_y has only 1 element, but we must include (0) 
	// at the end to return a double
}
// Function for calculating the volume integral of y
double integral(matrix y) {

	return 2*PI*(gl.I,r*r*y,leg.I_00)(0);

}

