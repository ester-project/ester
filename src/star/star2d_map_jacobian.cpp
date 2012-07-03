#include"star.h"

void star2d::add_dr(solver *op,const char *eqn,const matrix &q) {

	op->add_d(eqn,"eta",q*map.J[0]);
	op->add_d(eqn,"deta",q*map.J[1]);
	op->add_d(eqn,"Ri",q*map.J[2]);
	op->add_d(eqn,"dRi",q*map.J[3]);

}

void star2d::add_drz(solver *op,const char *eqn,const matrix &q) {

	op->add_d(eqn,"eta",q*(D,map.J[0]));
	op->add_d(eqn,"deta",q*(D,map.J[1]));
	op->add_d(eqn,"Ri",q*(D,map.J[2]));
	op->add_d(eqn,"dRi",q*(D,map.J[3]));

}

void star2d::add_drzz(solver *op,const char *eqn,const matrix &q) {

	op->add_d(eqn,"eta",q*(D,D,map.J[0]));
	op->add_d(eqn,"deta",q*(D,D,map.J[1]));
	op->add_d(eqn,"Ri",q*(D,D,map.J[2]));
	op->add_d(eqn,"dRi",q*(D,D,map.J[3]));

}

void star2d::add_drt(solver *op,const char *eqn,const matrix &q) {

	op->add_d(eqn,"eta",q*(map.J[0],Dt));
	op->add_d(eqn,"deta",q*(map.J[1],Dt));
	op->add_d(eqn,"Ri",q*(map.J[2],Dt));
	op->add_d(eqn,"dRi",q*(map.J[3],Dt));
	op->add_r(eqn,"Ri",q*map.J[2],Dt);
	op->add_r(eqn,"dRi",q*map.J[3],Dt);

}

void star2d::add_drtt(solver *op,const char *eqn,const matrix &q) {

	op->add_d(eqn,"eta",q*(map.J[0],Dt2));
	op->add_d(eqn,"deta",q*(map.J[1],Dt2));
	op->add_d(eqn,"Ri",q*(map.J[2],Dt2));
	op->add_d(eqn,"dRi",q*(map.J[3],Dt2));
	op->add_r(eqn,"Ri",q*map.J[2],Dt2);
	op->add_r(eqn,"dRi",q*map.J[3],Dt2);
	op->add_r(eqn,"Ri",2*q*(map.J[2],Dt),Dt);
	op->add_r(eqn,"dRi",2*q*(map.J[3],Dt),Dt);

}

void star2d::add_drzt(solver *op,const char *eqn,const matrix &q) {

	op->add_d(eqn,"eta",q*(D,map.J[0],Dt));
	op->add_d(eqn,"deta",q*(D,map.J[1],Dt));
	op->add_d(eqn,"Ri",q*(D,map.J[2],Dt));
	op->add_d(eqn,"dRi",q*(D,map.J[3],Dt));
	op->add_r(eqn,"Ri",q*(D,map.J[2]),Dt);
	op->add_r(eqn,"dRi",q*(D,map.J[3]),Dt);

}

void star2d::add_drex(solver *op,const char *eqn,const matrix &q) {

	op->add_d(ndomains(),eqn,"eta",q*map.ex.J[0]);
	op->add_d(ndomains(),eqn,"Ri",q*map.ex.J[2]);

}

void star2d::add_drzex(solver *op,const char *eqn,const matrix &q) {

	op->add_d(ndomains(),eqn,"eta",q*(Dex,map.ex.J[0]));
	op->add_d(ndomains(),eqn,"Ri",q*(Dex,map.ex.J[2]));

}

void star2d::add_drzzex(solver *op,const char *eqn,const matrix &q) {

	op->add_d(ndomains(),eqn,"eta",q*(Dex,Dex,map.ex.J[0]));
	op->add_d(ndomains(),eqn,"Ri",q*(Dex,Dex,map.ex.J[2]));

}

void star2d::add_drtex(solver *op,const char *eqn,const matrix &q) {

	op->add_d(ndomains(),eqn,"eta",q*(map.ex.J[0],Dt));
	op->add_d(ndomains(),eqn,"Ri",q*(map.ex.J[2],Dt));
	op->add_r(ndomains(),eqn,"Ri",q*map.ex.J[2],Dt);

}

void star2d::add_drttex(solver *op,const char *eqn,const matrix &q) {

	op->add_d(ndomains(),eqn,"eta",q*(map.ex.J[0],Dt2));
	op->add_d(ndomains(),eqn,"Ri",q*(map.ex.J[2],Dt2));
	op->add_r(ndomains(),eqn,"Ri",q*map.ex.J[2],Dt2);
	op->add_r(ndomains(),eqn,"Ri",2*q*(map.ex.J[2],Dt),Dt);

}

void star2d::add_drztex(solver *op,const char *eqn,const matrix &q) {

	op->add_d(ndomains(),eqn,"eta",q*(Dex,map.ex.J[0],Dt));
	op->add_d(ndomains(),eqn,"Ri",q*(Dex,map.ex.J[2],Dt));
	op->add_r(ndomains(),eqn,"Ri",q*(Dex,map.ex.J[2]),Dt);

}

void star2d::add_dr(solver *op,int iblock,const char *eqn,const matrix &q) {

	if(iblock==ndomains()) {
		add_drex(op,eqn,q);
		return;
	}
	
	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];	
	n=iblock;
		
	op->add_d(n,eqn,"eta",q*map.J[0].block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
	op->add_d(n,eqn,"deta",q*map.J[1].block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
	op->add_d(n,eqn,"Ri",q*map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
	op->add_d(n,eqn,"dRi",q*map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1));

}

void star2d::add_drz(solver *op,int iblock,const char *eqn,const matrix &q) {

	if(iblock==ndomains()) {
		add_drzex(op,eqn,q);
		return;
	}

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;

	op->add_d(n,eqn,"eta",q*(D.block(n),map.J[0].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));
	op->add_d(n,eqn,"deta",q*(D.block(n),map.J[1].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));
	op->add_d(n,eqn,"Ri",q*(D.block(n),map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));
	op->add_d(n,eqn,"dRi",q*(D.block(n),map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));

}

void star2d::add_drzz(solver *op,int iblock,const char *eqn,const matrix &q) {

	if(iblock==ndomains()) {
		add_drzzex(op,eqn,q);
		return;
	}

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	
	op->add_d(n,eqn,"eta",q*(D.block(n),D.block(n),map.J[0].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));
	op->add_d(n,eqn,"deta",q*(D.block(n),D.block(n),map.J[1].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));
	op->add_d(n,eqn,"Ri",q*(D.block(n),D.block(n),map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));
	op->add_d(n,eqn,"dRi",q*(D.block(n),D.block(n),map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)));

}

void star2d::add_drt(solver *op,int iblock,const char *eqn,const matrix &q) {

	if(iblock==ndomains()) {
		add_drtex(op,eqn,q);
		return;
	}

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	
	op->add_d(n,eqn,"eta",q*(map.J[0].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_d(n,eqn,"deta",q*(map.J[1].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_d(n,eqn,"Ri",q*(map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_d(n,eqn,"dRi",q*(map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_r(n,eqn,"Ri",q*map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt);
	op->add_r(n,eqn,"dRi",q*map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt);

}

void star2d::add_drtt(solver *op,int iblock,const char *eqn,const matrix &q) {

	if(iblock==ndomains()) {
		add_drttex(op,eqn,q);
		return;
	}

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;

	op->add_d(n,eqn,"eta",q*(map.J[0].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt2));
	op->add_d(n,eqn,"deta",q*(map.J[1].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt2));
	op->add_d(n,eqn,"Ri",q*(map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt2));
	op->add_d(n,eqn,"dRi",q*(map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt2));
	op->add_r(n,eqn,"Ri",q*map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt2);
	op->add_r(n,eqn,"dRi",q*map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt2);
	op->add_r(n,eqn,"Ri",2*q*(map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt),Dt);
	op->add_r(n,eqn,"dRi",2*q*(map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt),Dt);

}

void star2d::add_drzt(solver *op,int iblock,const char *eqn,const matrix &q) {

	if(iblock==ndomains()) {
		add_drztex(op,eqn,q);
		return;
	}

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	
	op->add_d(n,eqn,"eta",q*(D.block(n),map.J[0].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_d(n,eqn,"deta",q*(D.block(n),map.J[1].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_d(n,eqn,"Ri",q*(D.block(n),map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_d(n,eqn,"dRi",q*(D.block(n),map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1),Dt));
	op->add_r(n,eqn,"Ri",q*(D.block(n),map.J[2].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)),Dt);
	op->add_r(n,eqn,"dRi",q*(D.block(n),map.J[3].block(j0,j0+map.gl.npts[n]-1,0,nth()-1)),Dt);

}

void star2d::add_bc_bot2_dr(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
		
	if(iblock<ndomains()) {
		op->bc_bot2_add_d(n,eqn,"eta",q*map.J[0].row(j0));
		op->bc_bot2_add_d(n,eqn,"deta",q*map.J[1].row(j0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*map.J[2].row(j0));
		op->bc_bot2_add_d(n,eqn,"dRi",q*map.J[3].row(j0));
	} else {
		op->bc_bot2_add_d(n,eqn,"eta",q*map.ex.J[0].row(0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*map.ex.J[2].row(0));
	}
	
}

void star2d::add_bc_bot2_drz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;

	if(iblock<ndomains()) {
		op->bc_bot2_add_d(n,eqn,"eta",q*(D,map.J[0]).row(j0));
		op->bc_bot2_add_d(n,eqn,"deta",q*(D,map.J[1]).row(j0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(D,map.J[2]).row(j0));
		op->bc_bot2_add_d(n,eqn,"dRi",q*(D,map.J[3]).row(j0));
	} else {
		op->bc_bot2_add_d(n,eqn,"eta",q*(Dex,map.ex.J[0]).row(0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(Dex,map.ex.J[2]).row(0));
	}
	
}

void star2d::add_bc_bot2_drzz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;	
		
	if(iblock<ndomains()) {
		op->bc_bot2_add_d(n,eqn,"eta",q*(D,D,map.J[0]).row(j0));
		op->bc_bot2_add_d(n,eqn,"deta",q*(D,D,map.J[1]).row(j0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(D,D,map.J[2]).row(j0));
		op->bc_bot2_add_d(n,eqn,"dRi",q*(D,D,map.J[3]).row(j0));
	} else {
		op->bc_bot2_add_d(n,eqn,"eta",q*(Dex,Dex,map.ex.J[0]).row(0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(Dex,Dex,map.ex.J[2]).row(0));
	}
	
}

void star2d::add_bc_bot2_drt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	
	if(iblock<ndomains()) {
		op->bc_bot2_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt));
		op->bc_bot2_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt));
		op->bc_bot2_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt));
		op->bc_bot2_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt);
		op->bc_bot2_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt);
	} else {
		op->bc_bot2_add_d(n,eqn,"eta",q*(map.ex.J[0].row(0),Dt));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(map.ex.J[2].row(0),Dt));
		op->bc_bot2_add_r(n,eqn,"Ri",q*map.ex.J[2].row(0),Dt);
	}
	
}

void star2d::add_bc_bot2_drtt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	
	if(iblock<ndomains()) {
		op->bc_bot2_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt2));
		op->bc_bot2_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt2));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt2));
		op->bc_bot2_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt2));
		op->bc_bot2_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt2);
		op->bc_bot2_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt2);
		op->bc_bot2_add_r(n,eqn,"Ri",2*q*(map.J[2].row(j0),Dt),Dt);
		op->bc_bot2_add_r(n,eqn,"dRi",2*q*(map.J[3].row(j0),Dt),Dt);
	} else {
		op->bc_bot2_add_d(n,eqn,"eta",q*(map.ex.J[0].row(0),Dt2));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(map.ex.J[2].row(0),Dt2));
		op->bc_bot2_add_r(n,eqn,"Ri",q*map.ex.J[2].row(0),Dt2);
		op->bc_bot2_add_r(n,eqn,"Ri",2*q*(map.ex.J[2].row(0),Dt),Dt);
	}
	
}

void star2d::add_bc_bot2_drzt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	
	if(iblock<ndomains()) {
		op->bc_bot2_add_d(n,eqn,"eta",q*(D,map.J[0],Dt).row(j0));
		op->bc_bot2_add_d(n,eqn,"deta",q*(D,map.J[1],Dt).row(j0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(D,map.J[2],Dt).row(j0));
		op->bc_bot2_add_d(n,eqn,"dRi",q*(D,map.J[3],Dt).row(j0));
		op->bc_bot2_add_r(n,eqn,"Ri",q*(D,map.J[2]).row(j0),Dt);
		op->bc_bot2_add_r(n,eqn,"dRi",q*(D,map.J[3]).row(j0),Dt);
	} else {
		op->bc_bot2_add_d(n,eqn,"eta",q*(Dex,map.ex.J[0],Dt).row(0));
		op->bc_bot2_add_d(n,eqn,"Ri",q*(Dex,map.ex.J[2],Dt).row(0));
		op->bc_bot2_add_r(n,eqn,"Ri",q*(Dex,map.ex.J[2]).row(0),Dt);
	}
	
}

void star2d::add_bc_top1_dr(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	if(n<ndomains()) j0+=map.gl.npts[n]-1;
	
	if(iblock<ndomains()) {
		op->bc_top1_add_d(n,eqn,"eta",q*map.J[0].row(j0));
		op->bc_top1_add_d(n,eqn,"deta",q*map.J[1].row(j0));
		op->bc_top1_add_d(n,eqn,"Ri",q*map.J[2].row(j0));
		op->bc_top1_add_d(n,eqn,"dRi",q*map.J[3].row(j0));
	} else {
		op->bc_top1_add_d(n,eqn,"eta",q*map.ex.J[0].row(nex()-1));
		op->bc_top1_add_d(n,eqn,"Ri",q*map.ex.J[2].row(nex()-1));
	}
	
}

void star2d::add_bc_top1_drz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	if(n<ndomains()) j0+=map.gl.npts[n]-1;
	
	if(iblock<ndomains()) {
		op->bc_top1_add_d(n,eqn,"eta",q*(D,map.J[0]).row(j0));
		op->bc_top1_add_d(n,eqn,"deta",q*(D,map.J[1]).row(j0));
		op->bc_top1_add_d(n,eqn,"Ri",q*(D,map.J[2]).row(j0));
		op->bc_top1_add_d(n,eqn,"dRi",q*(D,map.J[3]).row(j0));
	} else {
		op->bc_top1_add_d(n,eqn,"eta",q*(Dex,map.ex.J[0]).row(nex()-1));
		op->bc_top1_add_d(n,eqn,"Ri",q*(Dex,map.ex.J[2]).row(nex()-1));
	}
	
}

void star2d::add_bc_top1_drzz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;	
	if(n<ndomains()) j0+=map.gl.npts[n]-1;
		
	if(iblock<ndomains()) {
		op->bc_top1_add_d(n,eqn,"eta",q*(D,D,map.J[0]).row(j0));
		op->bc_top1_add_d(n,eqn,"deta",q*(D,D,map.J[1]).row(j0));
		op->bc_top1_add_d(n,eqn,"Ri",q*(D,D,map.J[2]).row(j0));
		op->bc_top1_add_d(n,eqn,"dRi",q*(D,D,map.J[3]).row(j0));
	} else {
		op->bc_top1_add_d(n,eqn,"eta",q*(Dex,Dex,map.ex.J[0]).row(nex()-1));
		op->bc_top1_add_d(n,eqn,"Ri",q*(Dex,Dex,map.ex.J[2]).row(nex()-1));
	}
	
}

void star2d::add_bc_top1_drt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	if(n<ndomains()) j0+=map.gl.npts[n]-1;
	
	if(iblock<ndomains()) {
		op->bc_top1_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt));
		op->bc_top1_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt));
		op->bc_top1_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt));
		op->bc_top1_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt));
		op->bc_top1_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt);
		op->bc_top1_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt);
	} else {
		op->bc_top1_add_d(n,eqn,"eta",q*(map.ex.J[0].row(nex()-1),Dt));
		op->bc_top1_add_d(n,eqn,"Ri",q*(map.ex.J[2].row(nex()-1),Dt));
		op->bc_top1_add_r(n,eqn,"Ri",q*map.ex.J[2].row(nex()-1),Dt);
	}
	
}

void star2d::add_bc_top1_drtt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	if(n<ndomains()) j0+=map.gl.npts[n]-1;
	
	if(iblock<ndomains()) {
		op->bc_top1_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt2));
		op->bc_top1_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt2));
		op->bc_top1_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt2));
		op->bc_top1_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt2));
		op->bc_top1_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt2);
		op->bc_top1_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt2);
		op->bc_top1_add_r(n,eqn,"Ri",2*q*(map.J[2].row(j0),Dt),Dt);
		op->bc_top1_add_r(n,eqn,"dRi",2*q*(map.J[3].row(j0),Dt),Dt);
	} else {
		op->bc_top1_add_d(n,eqn,"eta",q*(map.ex.J[0].row(nex()-1),Dt2));
		op->bc_top1_add_d(n,eqn,"Ri",q*(map.ex.J[2].row(nex()-1),Dt2));
		op->bc_top1_add_r(n,eqn,"Ri",q*map.ex.J[2].row(nex()-1),Dt2);
		op->bc_top1_add_r(n,eqn,"Ri",2*q*(map.ex.J[2].row(nex()-1),Dt),Dt);
	}
	
}

void star2d::add_bc_top1_drzt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	if(n<ndomains()) j0+=map.gl.npts[n]-1;
	
	if(iblock<ndomains()) {
		op->bc_top1_add_d(n,eqn,"eta",q*(D,map.J[0],Dt).row(j0));
		op->bc_top1_add_d(n,eqn,"deta",q*(D,map.J[1],Dt).row(j0));
		op->bc_top1_add_d(n,eqn,"Ri",q*(D,map.J[2],Dt).row(j0));
		op->bc_top1_add_d(n,eqn,"dRi",q*(D,map.J[3],Dt).row(j0));
		op->bc_top1_add_r(n,eqn,"Ri",q*(D,map.J[2]).row(j0),Dt);
		op->bc_top1_add_r(n,eqn,"dRi",q*(D,map.J[3]).row(j0),Dt);
	} else {
		op->bc_top1_add_d(n,eqn,"eta",q*(Dex,map.ex.J[0],Dt).row(nex()-1));
		op->bc_top1_add_d(n,eqn,"Ri",q*(Dex,map.ex.J[2],Dt).row(nex()-1));
		op->bc_top1_add_r(n,eqn,"Ri",q*(Dex,map.ex.J[2]).row(nex()-1),Dt);
	}
	
}

void star2d::add_bc_bot1_dr(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0--;
		
	op->bc_bot1_add_d(n,eqn,"eta",q*map.J[0].row(j0));
	op->bc_bot1_add_d(n,eqn,"deta",q*map.J[1].row(j0));
	op->bc_bot1_add_d(n,eqn,"Ri",q*map.J[2].row(j0));
	op->bc_bot1_add_d(n,eqn,"dRi",q*map.J[3].row(j0));
	
}

void star2d::add_bc_bot1_drz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0--;

	op->bc_bot1_add_d(n,eqn,"eta",q*(D,map.J[0]).row(j0));
	op->bc_bot1_add_d(n,eqn,"deta",q*(D,map.J[1]).row(j0));
	op->bc_bot1_add_d(n,eqn,"Ri",q*(D,map.J[2]).row(j0));
	op->bc_bot1_add_d(n,eqn,"dRi",q*(D,map.J[3]).row(j0));
		
}

void star2d::add_bc_bot1_drzz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0--;	
		
	op->bc_bot1_add_d(n,eqn,"eta",q*(D,D,map.J[0]).row(j0));
	op->bc_bot1_add_d(n,eqn,"deta",q*(D,D,map.J[1]).row(j0));
	op->bc_bot1_add_d(n,eqn,"Ri",q*(D,D,map.J[2]).row(j0));
	op->bc_bot1_add_d(n,eqn,"dRi",q*(D,D,map.J[3]).row(j0));
	
}

void star2d::add_bc_bot1_drt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0--;
	
	op->bc_bot1_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt));
	op->bc_bot1_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt));
	op->bc_bot1_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt));
	op->bc_bot1_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt));
	op->bc_bot1_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt);
	op->bc_bot1_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt);
	
}

void star2d::add_bc_bot1_drtt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0--;
	
	op->bc_bot1_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt2));
	op->bc_bot1_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt2));
	op->bc_bot1_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt2));
	op->bc_bot1_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt2));
	op->bc_bot1_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt2);
	op->bc_bot1_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt2);
	op->bc_bot1_add_r(n,eqn,"Ri",2*q*(map.J[2].row(j0),Dt),Dt);
	op->bc_bot1_add_r(n,eqn,"dRi",2*q*(map.J[3].row(j0),Dt),Dt);
	
}

void star2d::add_bc_bot1_drzt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0--;
	
	op->bc_bot1_add_d(n,eqn,"eta",q*(D,map.J[0],Dt).row(j0));
	op->bc_bot1_add_d(n,eqn,"deta",q*(D,map.J[1],Dt).row(j0));
	op->bc_bot1_add_d(n,eqn,"Ri",q*(D,map.J[2],Dt).row(j0));
	op->bc_bot1_add_d(n,eqn,"dRi",q*(D,map.J[3],Dt).row(j0));
	op->bc_bot1_add_r(n,eqn,"Ri",q*(D,map.J[2]).row(j0),Dt);
	op->bc_bot1_add_r(n,eqn,"dRi",q*(D,map.J[3]).row(j0),Dt);
	
}

void star2d::add_bc_top2_dr(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0+=map.gl.npts[n];
	
	if(iblock<ndomains()-1) {
		op->bc_top2_add_d(n,eqn,"eta",q*map.J[0].row(j0));
		op->bc_top2_add_d(n,eqn,"deta",q*map.J[1].row(j0));
		op->bc_top2_add_d(n,eqn,"Ri",q*map.J[2].row(j0));
		op->bc_top2_add_d(n,eqn,"dRi",q*map.J[3].row(j0));
	} else {
		op->bc_top2_add_d(n,eqn,"eta",q*map.ex.J[0].row(0));
		op->bc_top2_add_d(n,eqn,"Ri",q*map.ex.J[2].row(0));
	}
	
}

void star2d::add_bc_top2_drz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0+=map.gl.npts[n];
	
	if(iblock<ndomains()-1) {
		op->bc_top2_add_d(n,eqn,"eta",q*(D,map.J[0]).row(j0));
		op->bc_top2_add_d(n,eqn,"deta",q*(D,map.J[1]).row(j0));
		op->bc_top2_add_d(n,eqn,"Ri",q*(D,map.J[2]).row(j0));
		op->bc_top2_add_d(n,eqn,"dRi",q*(D,map.J[3]).row(j0));
	} else {
		op->bc_top2_add_d(n,eqn,"eta",q*(Dex,map.ex.J[0]).row(0));
		op->bc_top2_add_d(n,eqn,"Ri",q*(Dex,map.ex.J[2]).row(0));
	}
	
}

void star2d::add_bc_top2_drzz(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;	
	j0+=map.gl.npts[n];
		
	if(iblock<ndomains()-1) {
		op->bc_top2_add_d(n,eqn,"eta",q*(D,D,map.J[0]).row(j0));
		op->bc_top2_add_d(n,eqn,"deta",q*(D,D,map.J[1]).row(j0));
		op->bc_top2_add_d(n,eqn,"Ri",q*(D,D,map.J[2]).row(j0));
		op->bc_top2_add_d(n,eqn,"dRi",q*(D,D,map.J[3]).row(j0));
	} else {
		op->bc_top2_add_d(n,eqn,"eta",q*(Dex,Dex,map.ex.J[0]).row(0));
		op->bc_top2_add_d(n,eqn,"Ri",q*(Dex,Dex,map.ex.J[2]).row(0));
	}
	
}

void star2d::add_bc_top2_drt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0+=map.gl.npts[n];
	
	if(iblock<ndomains()-1) {
		op->bc_top2_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt));
		op->bc_top2_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt));
		op->bc_top2_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt));
		op->bc_top2_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt));
		op->bc_top2_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt);
		op->bc_top2_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt);
	} else {
		op->bc_top2_add_d(n,eqn,"eta",q*(map.ex.J[0].row(0),Dt));
		op->bc_top2_add_d(n,eqn,"Ri",q*(map.ex.J[2].row(0),Dt));
		op->bc_top2_add_r(n,eqn,"Ri",q*map.ex.J[2].row(0),Dt);
	}
	
}

void star2d::add_bc_top2_drtt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0+=map.gl.npts[n];
	
	if(iblock<ndomains()-1) {
		op->bc_top2_add_d(n,eqn,"eta",q*(map.J[0].row(j0),Dt2));
		op->bc_top2_add_d(n,eqn,"deta",q*(map.J[1].row(j0),Dt2));
		op->bc_top2_add_d(n,eqn,"Ri",q*(map.J[2].row(j0),Dt2));
		op->bc_top2_add_d(n,eqn,"dRi",q*(map.J[3].row(j0),Dt2));
		op->bc_top2_add_r(n,eqn,"Ri",q*map.J[2].row(j0),Dt2);
		op->bc_top2_add_r(n,eqn,"dRi",q*map.J[3].row(j0),Dt2);
		op->bc_top2_add_r(n,eqn,"Ri",2*q*(map.J[2].row(j0),Dt),Dt);
		op->bc_top2_add_r(n,eqn,"dRi",2*q*(map.J[3].row(j0),Dt),Dt);
	} else {
		op->bc_top2_add_d(n,eqn,"eta",q*(map.ex.J[0].row(0),Dt2));
		op->bc_top2_add_d(n,eqn,"Ri",q*(map.ex.J[2].row(0),Dt2));
		op->bc_top2_add_r(n,eqn,"Ri",q*map.ex.J[2].row(0),Dt2);
		op->bc_top2_add_r(n,eqn,"Ri",2*q*(map.ex.J[2].row(0),Dt),Dt);
	}
	
}

void star2d::add_bc_top2_drzt(solver *op,int iblock,const char *eqn,const matrix &q) {

	int n,j0=0;
	for(n=0;n<iblock;n++) j0+=map.gl.npts[n];
	n=iblock;
	j0+=map.gl.npts[n];
	
	if(iblock<ndomains()-1) {
		op->bc_top2_add_d(n,eqn,"eta",q*(D,map.J[0],Dt).row(j0));
		op->bc_top2_add_d(n,eqn,"deta",q*(D,map.J[1],Dt).row(j0));
		op->bc_top2_add_d(n,eqn,"Ri",q*(D,map.J[2],Dt).row(j0));
		op->bc_top2_add_d(n,eqn,"dRi",q*(D,map.J[3],Dt).row(j0));
		op->bc_top2_add_r(n,eqn,"Ri",q*(D,map.J[2]).row(j0),Dt);
		op->bc_top2_add_r(n,eqn,"dRi",q*(D,map.J[3]).row(j0),Dt);
	} else {
		op->bc_top2_add_d(n,eqn,"eta",q*(Dex,map.ex.J[0],Dt).row(0));
		op->bc_top2_add_d(n,eqn,"Ri",q*(Dex,map.ex.J[2],Dt).row(0));
		op->bc_top2_add_r(n,eqn,"Ri",q*(Dex,map.ex.J[2]).row(0),Dt);
	}
	
}



