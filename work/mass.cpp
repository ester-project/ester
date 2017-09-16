
matrix calcMass(star1d &a) {

	solver massSolver;
	massSolver.init(a.map.ndomains, 1, "full");
	massSolver.regvar("mass");
	massSolver.set_nr(a.map.npts);

	massSolver.add_l("mass", "mass", ones(a.nr, 1), a.D);
	matrix rhs = 4 * PI * a.r * a.r * a.rho;

	massSolver.bc_bot2_add_d(0, "mass", "mass", ones(1, 1));
	rhs(0) = 0;
	int jfirst = 0;
	for(int i = 1; i < a.ndomains; i++) {
		jfirst += a.map.gl.npts[i-1];
		massSolver.bc_bot2_add_d(i, "mass", "mass", ones(1, 1));
		massSolver.bc_bot1_add_d(i, "mass", "mass", -ones(1, 1));
		rhs(jfirst) = 0;
	}
	massSolver.set_rhs("mass", rhs);
	massSolver.solve();

	matrix mass = massSolver.get_var("mass");

	return mass;

}
