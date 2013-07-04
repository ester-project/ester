// The following example solves the differential equation
//		x^2*y'' + x*y' - y = x^2
// with boundary conditions y(0)=0 and y(1)=0
// whose exact solution is
//      y = x*(x-1)/3		

#include<stdio.h>
#include"numdiff.h"

int main() {
	
	//Initialize a diff_gl object with 1 domain 

	int n=20; // Number of points. 
		//In fact, this example can be solved using only 3 points.
	diff_gl gl(1);
	
	gl.set_npts(n);
	gl.set_xif(0.,1.);
	gl.init();
	
	// We will work with only 1 domain, so we create a reference to the 
	// first (and only) block of gl.D
	
	matrix &D=gl.D.block(0);
	matrix &x=gl.x;
	
	// Set up the operator matrix and the right hand side
	
	matrix op,rhs;
	
	op=x*x*(D,D)+x*D-eye(n);
	rhs=x*x;
	
	// Introduce boundary conditions
	
	op.setrow(0,zeros(1,n));op(0,0)=1;
	rhs(0)=0;
	op.setrow(-1,zeros(1,n));op(-1,-1)=1;
	rhs(-1)=0;
	
	// Solve the system
	
	matrix y;
	
	y=op.solve(rhs);
	
	// Interpolate the solution into a finer grid
	
	matrix x_fine,y_fine;
	
	x_fine=vector_t(0,1,100);
	y_fine=gl.eval(y,x_fine);
	
	// Compare with the exact solution
	
	matrix y_exact;
	
	y_exact=x_fine*(x_fine-1)/3;
	
	printf("Solved using %d points\n",gl.N);
	printf("Max. error=%e\n",max(abs(y_fine-y_exact)));
	
	return 0;	
	
}
