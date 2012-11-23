Utilities to calculate a 1d/2d stellar model using as input parameters the radius 
and effective temperature instead of mass and Xc.

Compile with

 $ make
 
This creates 2 executables and place them in $(ESTER_DIR)/bin. The executables are

- star1dR. 
	For calculating a non-rotating stellar model. The input parameters M and Xc have been replaced by 
		- R : Radius
		- Teff : Effective temperature.
	If these parameters are not set in the command-line or in a parameter file, the new model will 
	have the same radius and Teff that the input model.
	
- star2dR.
	For calculating a rotating stellar model. The input parameters M, Xc and Omega_bk have been replaced by 
		- Rp or R : Polar radius
		- Re : Equatorial radius
		- Teff : Effective temperature.
	If these parameters are not set in the command-line or in a parameter file, the new model will 
	have the same radius and Teff that the input model.
	
Known limitations:

	These routines are designed to refine an existing model. If the input model is not specified or
	it is too far from the solution, convergence is not guaranteed.
	
	The 2d version star2dR has problems to converge when the rotation velocity of the input model
	is zero, or it is much lower than the final velocity.
	
Recommended workflow:

	1. Calculate an initial model using star1d.
	2. Refine this model with star1dR to have the correct values of R and Teff.
	3. Use star2d to obtain the corresponding rotating model using an initial estimation of Omega_bk.
	4. Refine the model using star2dR.
	

	
