#include"solver.h"
#include<string.h>
#include<stdlib.h>


void solver::add_d(const char *eqn, const char *varn,const matrix &d) {
	
	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		block[i].add_d(ieq,ivar,d.block(j0,j1,0,d.ncols()-1));
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l) {

	int ieq,ivar;
	int i,j0;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<l.nblocks();i++) {
		block[i].add_l(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i));
		j0+=l.block(i).nrows();
	}
}
void solver::add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		block[i].add_r(ieq,ivar,d.block(j0,j1,0,d.ncols()-1),r);
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r) {

	int ieq,ivar;
	int i,j0;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<l.nblocks();i++) {
		block[i].add_lr(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),r);
		j0+=l.block(i).nrows();
	}
}
void solver::add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;j1=0;
	for(i=0;i<l.nblocks();i++) {
		block[i].add_li(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),ii.block(j1,j1+l.block(i).ncols()-1,0,ii.ncols()-1));
		j0+=l.block(i).nrows();
		j1+=l.block(i).ncols();
	}
}
void solver::add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		block[i].add_ri(ieq,ivar,d.block(j0,j1,0,d.ncols()-1),r,ii.block(j0,j1,0,ii.ncols()-1));
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	j1=0;
	for(i=0;i<l.nblocks();i++) {
		block[i].add_lri(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),r,ii.block(j1,j1+l.block(i).ncols()-1,0,ii.ncols()-1));
		j0+=l.block(i).nrows();
		j1+=l.block(i).ncols();		
	}
}


void solver::bc_pol_add_d(const char *eqn, const char *varn,const matrix &d) {
	
	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		bc_pol[i].add_d(ieq,ivar,d.block(j0,j1,0,d.ncols()-1));
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::bc_pol_add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l) {

	int ieq,ivar;
	int i,j0;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<l.nblocks();i++) {
		bc_pol[i].add_l(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i));
		j0+=l.block(i).nrows();
	}
}
void solver::bc_pol_add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		bc_pol[i].add_r(ieq,ivar,d.block(j0,j1,0,d.ncols()-1),r);
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::bc_pol_add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r) {

	int ieq,ivar;
	int i,j0;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<l.nblocks();i++) {
		bc_pol[i].add_lr(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),r);
		j0+=l.block(i).nrows();
	}
}
void solver::bc_pol_add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;j1=0;
	for(i=0;i<l.nblocks();i++) {
		bc_pol[i].add_li(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),ii.block(j1,j1+l.block(i).ncols()-1,0,ii.ncols()-1));
		j0+=l.block(i).nrows();
		j1+=l.block(i).ncols();
	}
}
void solver::bc_pol_add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		bc_pol[i].add_ri(ieq,ivar,d.block(j0,j1,0,d.ncols()-1),r,ii.block(j0,j1,0,ii.ncols()-1));
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::bc_pol_add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	j1=0;
	for(i=0;i<l.nblocks();i++) {
		bc_pol[i].add_lri(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),r,ii.block(j1,j1+l.block(i).ncols()-1,0,ii.ncols()-1));
		j0+=l.block(i).nrows();
		j1+=l.block(i).ncols();		
	}
}

void solver::bc_eq_add_d(const char *eqn, const char *varn,const matrix &d) {
	
	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		bc_eq[i].add_d(ieq,ivar,d.block(j0,j1,0,d.ncols()-1));
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::bc_eq_add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l) {

	int ieq,ivar;
	int i,j0;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<l.nblocks();i++) {
		bc_eq[i].add_l(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i));
		j0+=l.block(i).nrows();
	}
}
void solver::bc_eq_add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		bc_eq[i].add_r(ieq,ivar,d.block(j0,j1,0,d.ncols()-1),r);
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::bc_eq_add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r) {

	int ieq,ivar;
	int i,j0;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<l.nblocks();i++) {
		bc_eq[i].add_lr(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),r);
		j0+=l.block(i).nrows();
	}
}
void solver::bc_eq_add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;j1=0;
	for(i=0;i<l.nblocks();i++) {
		bc_eq[i].add_li(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),ii.block(j1,j1+l.block(i).ncols()-1,0,ii.ncols()-1));
		j0+=l.block(i).nrows();
		j1+=l.block(i).ncols();
	}
}
void solver::bc_eq_add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	for(i=0;i<nb;i++) {
		j1=j0+def_nr[i]-1;
		bc_eq[i].add_ri(ieq,ivar,d.block(j0,j1,0,d.ncols()-1),r,ii.block(j0,j1,0,ii.ncols()-1));
		j0=j1+1;
		if(j0>=d.nrows()) break;
	}
}
void solver::bc_eq_add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &ii) {

	int ieq,ivar;
	int i,j0,j1;
	
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	j0=0;
	j1=0;
	for(i=0;i<l.nblocks();i++) {
		bc_eq[i].add_lri(ieq,ivar,d.block(j0,j0+l.block(i).nrows()-1,0,d.ncols()-1),l.block(i),r,ii.block(j1,j1+l.block(i).ncols()-1,0,ii.ncols()-1));
		j0+=l.block(i).nrows();
		j1+=l.block(i).ncols();		
	}
}

void solver::add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_d(ieq,ivar,d);
}
void solver::add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_l(ieq,ivar,d,l);
}
void solver::add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_r(ieq,ivar,d,r);
}
void solver::add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	block[iblock].add_lri(ieq,ivar,d,l,r,i);
}


void solver::bc_pol_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_d(ieq,ivar,d);
}
void solver::bc_pol_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_l(ieq,ivar,d,l);
}
void solver::bc_pol_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_r(ieq,ivar,d,r);
}
void solver::bc_pol_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::bc_pol_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::bc_pol_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::bc_pol_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_pol[iblock].add_lri(ieq,ivar,d,l,r,i);
}


void solver::bc_eq_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_d(ieq,ivar,d);
}
void solver::bc_eq_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_l(ieq,ivar,d,l);
}
void solver::bc_eq_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_r(ieq,ivar,d,r);
}
void solver::bc_eq_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::bc_eq_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::bc_eq_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::bc_eq_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_eq[iblock].add_lri(ieq,ivar,d,l,r,i);
}


void solver::bc_bot1_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_d(ieq,ivar,d);
}
void solver::bc_bot1_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_l(ieq,ivar,d,l);
}
void solver::bc_bot1_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_r(ieq,ivar,d,r);
}
void solver::bc_bot1_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::bc_bot1_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::bc_bot1_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::bc_bot1_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot1[iblock].add_lri(ieq,ivar,d,l,r,i);
}

void solver::bc_bot2_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_d(ieq,ivar,d);
}
void solver::bc_bot2_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_l(ieq,ivar,d,l);
}
void solver::bc_bot2_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_r(ieq,ivar,d,r);
}
void solver::bc_bot2_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::bc_bot2_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::bc_bot2_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::bc_bot2_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_bot2[iblock].add_lri(ieq,ivar,d,l,r,i);
}

void solver::bc_top1_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_d(ieq,ivar,d);
}
void solver::bc_top1_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_l(ieq,ivar,d,l);
}
void solver::bc_top1_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_r(ieq,ivar,d,r);
}
void solver::bc_top1_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::bc_top1_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::bc_top1_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::bc_top1_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top1[iblock].add_lri(ieq,ivar,d,l,r,i);
}

void solver::bc_top2_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_d(ieq,ivar,d);
}
void solver::bc_top2_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_l(ieq,ivar,d,l);
}
void solver::bc_top2_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_r(ieq,ivar,d,r);
}
void solver::bc_top2_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_lr(ieq,ivar,d,l,r);
}
void solver::bc_top2_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_li(ieq,ivar,d,l,i);
}
void solver::bc_top2_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_ri(ieq,ivar,d,r,i);
}
void solver::bc_top2_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);
	sync=0;
	bc_top2[iblock].add_lri(ieq,ivar,d,l,r,i);
}

