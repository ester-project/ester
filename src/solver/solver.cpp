#include"solver.h"
#include<string.h>
#include<stdlib.h>
extern "C" {
#include CBLAS
}

void solver_block::init(int nvar) {
	
	int i,j;
	
	eq=new solver_elem**[nvar];
	for(i=0;i<nvar;i++) {
		eq[i]=new solver_elem*[nvar];
		for(j=0;j<nvar;j++) eq[i][j]=NULL;
	}
	eq_set=zeros(nvar,1);
	nv=nvar;
	nr=new int[nv];
	nth=new int[nv];
	for(i=0;i<nv;i++) {nr[i]=0;nth[i]=0;}

}

void solver_block::destroy() {

	int i;
	
	reset();
	
	for(i=0;i<nv;i++) delete [] eq[i];
	delete [] eq;
	delete [] nr;
	delete [] nth;
	
}

void solver_block::reset() {

	int i;
	
	for(i=0;i<nv;i++) {
		reset(i);
	}
}
	
void solver_block::reset(int ieq) {

	solver_elem *p;
	solver_elem *p0;
	int j;
	
	for(j=0;j<nv;j++) {
		while(eq[ieq][j]!=NULL) {
			p=eq[ieq][j];
			if(p->next==NULL) {
				delete p;
				eq[ieq][j]=NULL;
			} else {
				while(p->next!=NULL) {p0=p;p=p->next;}
				delete p;
				p0->next=NULL;
			}
		}
	}
	eq_set(ieq)=0;
	nr[ieq]=0;
	nth[ieq]=0;
}
	 
void solver_block::add_d(int ieq,int ivar,const matrix &d) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='d';
	p->D=d;
	p->next=NULL;

	eq_set(ieq)=1;

}

void solver_block::add_l(int ieq,int ivar,const matrix &d,const matrix &l) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='l';
	p->D=d;
	p->L=l;
	p->next=NULL;

	eq_set(ieq)=1;

}

void solver_block::add_r(int ieq,int ivar,const matrix &d,const matrix &r) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='r';
	p->D=d;
	p->R=r;
	p->next=NULL;

	eq_set(ieq)=1;

}

void solver_block::add_lr(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='f';
	p->D=d;
	p->L=l;
	p->R=r;
	p->next=NULL;

	eq_set(ieq)=1;

}

void solver_block::add_li(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &i) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='m';
	p->D=d;
	p->L=l;
	p->I=i;
	p->next=NULL;

	eq_set(ieq)=1;

}

void solver_block::add_ri(int ieq,int ivar,const matrix &d,const matrix &r,const matrix &i) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='s';
	p->D=d;
	p->R=r;
	p->I=i;
	p->next=NULL;

	eq_set(ieq)=1;

}

void solver_block::add_lri(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {

	solver_elem *p,*pnew;
	
	nr[ieq]=d.nrows();
	nth[ieq]=d.ncols();
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
		p=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
		p=pnew;
	}
	p->type='g';
	p->D=d;
	p->L=l;
	p->R=r;
	p->I=i;
	p->next=NULL;

	eq_set(ieq)=1;

}

solver::solver() {
	initd=0;
	verbose=0;
	use_cgs=0;
	maxit_ref=0;
	maxit_cgs=10;
	rel_tol=1e-12;
	abs_tol=1e-20;
	def_nr=NULL;
	mult_fcn=NULL;
	mult_context=NULL;
};

void solver::init(int nblock,int nvar,const char *solver_type) {

	int i;

	if(initd) destroy();
	sync=0;
	initd=1;
	nv=nvar;
	nb=nblock;
	solver_N=0;
	strncpy(type,solver_type,20);
	type[20]='\0';
	op=NULL;
	
	block=new solver_block[nb];
	bc_eq=new solver_block[nb];
	bc_pol=new solver_block[nb];
	bc_bot1=new solver_block[nb];
	bc_bot2=new solver_block[nb];
	bc_top1=new solver_block[nb];
	bc_top2=new solver_block[nb];
	
	for(i=0;i<nb;i++) {
		block[i].init(nv);
		bc_eq[i].init(nv);
		bc_pol[i].init(nv);
		if(i) bc_bot1[i].init(nv);
		bc_bot2[i].init(nv);
		bc_top1[i].init(nv);
		if(i<nb-1) bc_top2[i].init(nv);
	}
	var_nr=new int*[nb];
	var_nbot=new int*[nb];
	var_ntop=new int*[nb];
	var_nth=new int*[nb];
	def_nr=new int[nb];
	for(i=0;i<nb;i++) {
		var_nr[i]=new int[nv];
		var_nbot[i]=new int[nv];
		var_ntop[i]=new int[nv];
		var_nth[i]=new int[nv];
	}
	var=new char*[nv];
	for(i=0;i<nv;i++) {
		var[i]=new char[32];var[i][0]='\0';
	}
	
	x=new matrix[nv];
	rhs=new matrix[nv];

}



void solver::destroy() {

	int i;
	
	initd=0;
	delete op;
	
	for(i=0;i<nb;i++) {
		block[i].destroy();
		bc_eq[i].destroy();
		bc_pol[i].destroy();
		if(i) bc_bot1[i].destroy();
		bc_bot2[i].destroy();
		bc_top1[i].destroy();
		if(i<nb-1) bc_top2[i].destroy();
		delete [] var_nr[i];
		delete [] var_nbot[i];
		delete [] var_ntop[i];
		delete [] var_nth[i];
	}

	delete [] block;
	delete [] bc_eq;
	delete [] bc_pol;
	delete [] bc_bot1;
	delete [] bc_bot2;
	delete [] bc_top1;
	delete [] bc_top2;
	
	delete [] var_nr;
	delete [] var_nbot;
	delete [] var_ntop;
	delete [] var_nth;
	delete [] def_nr;
	
	for(i=0;i<nv;i++) {
		delete [] var[i];
	}
	delete [] var;
	delete [] x;
	delete [] rhs;
	
}

void solver::reset() {

	int i;

	sync=0;
	for(i=0;i<nb;i++) {
		block[i].reset();
		bc_eq[i].reset();
		bc_pol[i].reset();
		if(i) bc_bot1[i].reset();
		bc_bot2[i].reset();
		bc_top1[i].reset();
		if(i<nb-1) bc_top2[i].reset();
	}
	
}

void solver::reset(int iblock) {

	sync=0;
	block[iblock].reset();
	bc_eq[iblock].reset();
	bc_pol[iblock].reset();
	if(iblock) bc_bot1[iblock].reset();
	bc_bot2[iblock].reset();
	bc_top1[iblock].reset();
	if(iblock<nb-1) bc_top2[iblock].reset();
	
}

void solver::reset(int iblock,int ieq) {

	sync=0;
	block[iblock].reset(ieq);
	bc_eq[iblock].reset(ieq);
	bc_pol[iblock].reset(ieq);
	if(iblock) bc_bot1[iblock].reset(ieq);
	bc_bot2[iblock].reset(ieq);
	bc_top1[iblock].reset(ieq);
	if(iblock<nb-1) bc_top2[iblock].reset(ieq);
	
}

void solver::set_nr(int *nr) {
	int i;

	for(i=0;i<nb;i++) {
		def_nr[i]=nr[i];
	}
}

void solver::regvar(const char *var_name) {

	int i,j;
	j=0;
	while (strlen(var[j])) {
		j++;
		if(j==nv) {
			printf("Can't register variable (increase nvar)\n");
			exit(1);
		}
	}	

	strncpy(var[j],var_name,31);	
	
}


int solver::get_nvar() {

	return nv;
	
}

int solver::get_nblocks() {

	return nb;
	
}

int solver::get_id(const char *varn) {

	int i=0;

	while(strcmp(varn,var[i])) {
		i++;
		if(i==nv) {
			printf("Unknown variable %s\n",varn);
			exit(1);
		}
	}
	return i;

}

void solver::check() {

	int i,j,n,nt;
	
	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(!block[i].eq_set(j)) {
				n=0;
				nt=0;
				if(bc_bot2[i].eq_set(j)) {n+=bc_bot2[i].nr[j];nt=bc_bot2[i].nth[j];}
				if(bc_top1[i].eq_set(j)) {n+=bc_top1[i].nr[j];nt=bc_top1[i].nth[j];}
				if(n) block[i].add_d(j,j,zeros(n,nt));
			}
		}
	}

}

void solver::set_rhs(const char *eqn,const matrix &b) {

	int ieq;
	ieq=get_id(eqn);
	rhs[ieq]=b;
}

matrix solver::get_rhs(const char *eqn) {
	
	int ieq;
	ieq=get_id(eqn);
	return rhs[ieq];
}

matrix solver::get_var(const char *varn) {
	
	int ivar;
	ivar=get_id(varn);
	return x[ivar];
}

matrix solver::get_rhs(int ieq) {
	
	return rhs[ieq];
}

matrix solver::get_var(int ivar) {
	
	return x[ivar];
}

//#define PERF_LOG

void solver::solve(int *info) {

	matrix y;
	matrix rhsw;
	int err;
	int i;
	
#ifdef PERF_LOG
	static tiempo tref,tcgs,tlu,tcreate,ttot;
	static int nref=0,ncgs=0,nlu=0,ncreate=0;
	ttot.start();
#endif
	if(info!=NULL) for(i=0;i<5;i++) info[i]=0;
	if(!sync) check();
	calc_struct();
	wrap(x,&y);
	wrap(rhs,&rhsw);
	if(op!=NULL&&sync==0&&use_cgs) {
		if(verbose)
			printf("CGS iteration:\n");
#ifdef PERF_LOG
		tcgs.start();
#endif
		y*=0;
		err=cgs(rhsw,y,maxit_cgs);
#ifdef PERF_LOG
		tcgs.stop();
		if(err>=0) ncgs+=err; else ncgs+=maxit_cgs;
#endif
		if(info!=NULL) {info[2]=1;info[4]=err;}
		if(err>=0) {
			unwrap(x,&y);
#ifdef PERF_LOG
			ttot.stop();
			printf("Total: %2.3fs \nLU(%d): %2.3fs (%2.3fs) cgs(%d): %2.3fs (%2.3fs)  ref(%d): %2.3fs (%2.3fs)  create(%d): %2.3fs (%2.3fs)\n",ttot.value(),nlu,tlu.value(),tlu.value()/nlu,ncgs,tcgs.value(),tcgs.value()/ncgs,nref,tref.value(),tref.value()/nref,ncreate,tcreate.value(),tcreate.value()/ncreate);
#endif
			return;
		}
		if(verbose)
			printf("Not converged after %d iterations\n",maxit_cgs);
	}
	if(!sync) {
		if(verbose)
			printf("Creating matrix...\n");
#ifdef PERF_LOG
		tcreate.start();
#endif
		create();
#ifdef PERF_LOG
		tcreate.stop();
		ncreate++;
#endif
		if(info!=NULL) info[0]=1;
	}
#ifdef PERF_LOG
	tlu.start();
#endif
	if(verbose)
		printf("LU factorization:\n");
	y=op->solve(rhsw);
#ifdef PERF_LOG
	tlu.stop();
	nlu++;
#endif
	if(maxit_ref) {
#ifdef PERF_LOG
		tref.start();
#endif
		if(verbose)
			printf("CGS refinement:\n");
		err=cgs(rhsw,y,maxit_ref);
#ifdef PERF_LOG
		tref.stop();
		if(err>=0) nref+=err; else nref+=maxit_cgs;
#endif
		if(info!=NULL) {info[1]=1;info[3]=err;}
		if(verbose)
			if(err<0)
				printf("CGS refinement not converged after %d iterations. Singular matrix?\n",maxit_ref);
	}
	unwrap(x,&y);
#ifdef PERF_LOG
	ttot.stop();
	printf("Total: %2.3fs \nLU(%d): %2.3fs (%2.3fs) cgs(%d): %2.3fs (%2.3fs)  ref(%d): %2.3fs (%2.3fs)  create(%d): %2.3fs (%2.3fs)\n",ttot.value(),nlu,tlu.value(),tlu.value()/nlu,ncgs,tcgs.value(),tcgs.value()/ncgs,nref,tref.value(),tref.value()/nref,ncreate,tcreate.value(),tcreate.value()/ncreate);
#endif	

}


void solver::mult(matrix *y) {

	if(mult_fcn!=NULL) (*mult_fcn)(y,mult_context);
	else mult_op(y);

}

void solver::mult_op(matrix *y) {

	matrix z[nv],yy,zz,zz0;
	int i,j,n,j0[nv];
	solver_elem *p;
	int **nr=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth,&N=solver_N;
	
	for(i=0;i<nv;i++) {
		z[i]=zeros(y[i].nrows(),y[i].ncols());
		for(j=0;j<nv;j++) j0[j]=0;
		for(n=0;n<nb;n++) {
			zz=zeros(nr[n][i],nth[n][i]);
			for(j=0;j<nv;j++) {
				yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
				p=block[n].eq[i][j];
				while(p!=NULL) {
					switch(p->type) {
						case 'd':
							zz+=p->D*yy;
							break;
						case 'l':
							zz+=p->D*(p->L,yy);
							break;
						case 'r':
							zz+=p->D*(yy,p->R);
							break;
						case 'f':
							zz+=p->D*(p->L,yy,p->R);
							break;
						case 'm':
							zz+=p->D*(p->L,(p->I*yy));
							break;
						case 's':
							zz+=p->D*((p->I*yy),p->R);
							break;
						case 'g':
							zz+=p->D*(p->L,(p->I*yy),p->R);
					}
					p=p->next;
				}
			}
			if(bc_pol[n].eq_set(i)) {
				zz0=zeros(nr[n][i],1);
				for(j=0;j<nv;j++) {
					yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
					p=bc_pol[n].eq[i][j];
					while(p!=NULL) {
						switch(p->type) {
							case 'd':
								zz0+=p->D*yy.col(nth[n][j]-1);
								break;
							case 'l':
								zz0+=p->D*(p->L,yy.col(nth[n][j]-1));
								break;
							case 'r':
								zz0+=p->D*(yy,p->R);
								break;
							case 'f':
								zz0+=p->D*(p->L,yy,p->R);
								break;
							case 'm':
								zz0+=p->D*(p->L,(p->I*yy.col(nth[n][j]-1)));
								break;
							case 's':
								zz0+=p->D*((p->I*yy),p->R);
								break;
							case 'g':
								zz0+=p->D*(p->L,(p->I*yy),p->R);
						}
						p=p->next;
					}
				}
				zz.setblock(0,zz0.nrows()-1,zz.ncols()-1,zz.ncols()-1,zz0);
			}
			if(bc_eq[n].eq_set(i)) {
				zz0=zeros(nr[n][i],1);
				for(j=0;j<nv;j++) {
					yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
					p=bc_eq[n].eq[i][j];
					while(p!=NULL) {
						switch(p->type) {
							case 'd':
								zz0+=p->D*yy.col(0);
								break;
							case 'l':
								zz0+=p->D*(p->L,yy.col(0));
								break;
							case 'r':
								zz0+=p->D*(yy,p->R);
								break;
							case 'f':
								zz0+=p->D*(p->L,yy,p->R);
								break;
							case 'm':
								zz0+=p->D*(p->L,(p->I*yy.col(0)));
								break;
							case 's':
								zz0+=p->D*((p->I*yy),p->R);
								break;
							case 'g':
								zz0+=p->D*(p->L,(p->I*yy),p->R);
						}
						p=p->next;
					}
				}
				zz.setblock(0,zz0.nrows()-1,0,0,zz0);
			}
			if(bc_bot2[n].eq_set(i)) {
				zz0=zeros(nbot[n][i],nth[n][i]);
				for(j=0;j<nv;j++) {
					yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
					p=bc_bot2[n].eq[i][j];
					while(p!=NULL) {
						switch(p->type) {
							case 'd':
								zz0+=p->D*yy.row(0);
								break;
							case 'l':
								zz0+=p->D*(p->L,yy);
								break;
							case 'r':
								zz0+=p->D*(yy.row(0),p->R);
								break;
							case 'f':
								zz0+=p->D*(p->L,yy,p->R);
								break;
							case 'm':
								zz0+=p->D*(p->L,(p->I*yy));
								break;
							case 's':
								zz0+=p->D*((p->I*yy.row(0)),p->R);
								break;
							case 'g':
								zz0+=p->D*(p->L,(p->I*yy),p->R);
						}
						p=p->next;
					}
				}
				zz.setblock(0,zz0.nrows()-1,0,zz0.ncols()-1,zz0);
			}
			if(n) if(bc_bot1[n].eq_set(i)) {
				zz0=zeros(nbot[n][i],nth[n][i]);
				for(j=0;j<nv;j++) {
					yy=y[j].block(j0[j]-nr[n-1][j],j0[j]-1,0,nth[n-1][j]-1);
					p=bc_bot1[n].eq[i][j];
					while(p!=NULL) {
						switch(p->type) {
							case 'd':
								zz0+=p->D*yy.row(yy.nrows()-1);
								break;
							case 'l':
								zz0+=p->D*(p->L,yy);
								break;
							case 'r':
								zz0+=p->D*(yy.row(yy.nrows()-1),p->R);
								break;
							case 'f':
								zz0+=p->D*(p->L,yy,p->R);
								break;
							case 'm':
								zz0+=p->D*(p->L,(p->I*yy));
								break;
							case 's':
								zz0+=p->D*((p->I*yy.row(yy.nrows()-1)),p->R);
								break;
							case 'g':
								zz0+=p->D*(p->L,(p->I*yy),p->R);
						}
						p=p->next;
					}
				}
				zz.setblock(0,zz0.nrows()-1,0,zz0.ncols()-1,
					zz.block(0,zz0.nrows()-1,0,zz0.ncols()-1)+zz0);
			}
			if(bc_top1[n].eq_set(i)) {
				zz0=zeros(ntop[n][i],nth[n][i]);
				for(j=0;j<nv;j++) {
					yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
					p=bc_top1[n].eq[i][j];
					while(p!=NULL) {
						switch(p->type) {
							case 'd':
								zz0+=p->D*yy.row(yy.nrows()-1);
								break;
							case 'l':
								zz0+=p->D*(p->L,yy);
								break;
							case 'r':
								zz0+=p->D*(yy.row(yy.nrows()-1),p->R);
								break;
							case 'f':
								zz0+=p->D*(p->L,yy,p->R);
								break;
							case 'm':
								zz0+=p->D*(p->L,(p->I*yy));
								break;
							case 's':
								zz0+=p->D*((p->I*yy.row(yy.nrows()-1)),p->R);
								break;
							case 'g':
								zz0+=p->D*(p->L,(p->I*yy),p->R);
						}
						p=p->next;
					}
				}
				zz.setblock(zz.nrows()-zz0.nrows(),zz.nrows()-1,0,zz0.ncols()-1,zz0);
			}
			if(n<nb-1) if(bc_top2[n].eq_set(i)) {
				zz0=zeros(ntop[n][i],nth[n][i]);
				for(j=0;j<nv;j++) {
					yy=y[j].block(j0[j]+nr[n][j],j0[j]+nr[n][j]+nr[n+1][j]-1,0,nth[n+1][j]-1);
					p=bc_top2[n].eq[i][j];
					while(p!=NULL) {
						switch(p->type) {
							case 'd':
								zz0+=p->D*yy.row(0);
								break;
							case 'l':
								zz0+=p->D*(p->L,yy);
								break;
							case 'r':
								zz0+=p->D*(yy.row(0),p->R);
								break;
							case 'f':
								zz0+=p->D*(p->L,yy,p->R);
								break;
							case 'm':
								zz0+=p->D*(p->L,(p->I*yy));
								break;
							case 's':
								zz0+=p->D*((p->I*yy.row(0)),p->R);
								break;
							case 'g':
								zz0+=p->D*(p->L,(p->I*yy),p->R);
						}
						p=p->next;
					}
				}
				zz.setblock(zz.nrows()-zz0.nrows(),zz.nrows()-1,0,zz0.ncols()-1,
					zz.block(zz.nrows()-zz0.nrows(),zz.nrows()-1,0,zz0.ncols()-1)+zz0);
			}
			z[i].setblock(j0[i],j0[i]+nr[n][i]-1,0,nth[n][i]-1,zz);
			for(j=0;j<nv;j++) j0[j]+=nr[n][j];
		}
	}
	for(i=0;i<nv;i++) y[i]=z[i];

}

int solver::cgs(const matrix &rhs,matrix &x,int maxit) {
	
	matrix r,r_,u,p,q,v,err_rel,err_abs,y[nv];
	int k,error=0,n,fin,i;
	double s,a,s_1;
	
	n=rhs.nrows();
	
	r=x;
	unwrap(y,&r);
	mult(y);
	wrap(y,&r);
	r=rhs-r;
	op->fwd_subs(r);
	r_=r;
	k=0;
	fin=0;
	while(!fin) {
		if(k>=maxit) {
			error=1;
			break;
		}
		s=cblas_ddot(n,r_.data(),1,r.data(),1);
		if(!s) {
			error=2;
			break;
		}
		if(!k) {
			u=r;
			p=r;
		} else {
			a=s/s_1;
			u=r+a*q;
			p=u+a*(q+a*p);
		}
		v=p;
		op->back_subs(v);
		unwrap(y,&v);
		mult(y);
		wrap(y,&v);
		op->fwd_subs(v);
		a=s/cblas_ddot(n,r_.data(),1,v.data(),1);
		q=u-a*v;
		v=u+q;
		op->back_subs(v);
		err_abs=abs(a*v);
		err_rel=abs(err_abs/(x+(x==0)));
		x=x+a*v;
		unwrap(y,&v);
		mult(y);
		wrap(y,&v);
		op->fwd_subs(v);
		r=r-a*v;
		s_1=s;
		fin=!exist((err_rel>rel_tol)&&(err_abs>abs_tol));
		if (verbose)
			printf("%d max. rel. error:%e max. abs. error:%e\n",k,max(err_rel*(err_abs>=abs_tol)),max(err_abs*(err_rel>=rel_tol)));
		k++;
	}
	if(error) {
		return -error;
	}
	
	return k;
}

void solver::create() {

	sync=1;
	if(!strcmp(type,"full")) {
		if(op==NULL) op=new solver_full(nb);
		create_full();
	}else if(!strcmp(type,"full-oc")) {
		if(op==NULL) op=new solver_full(nb,1);
		create_full();
	} else if(!strcmp(type,"iter")) {
		if(op==NULL) op=new solver_iter();
	} else {
		strcpy(type,"full");
		if(op==NULL) op=new solver_full(nb);
		create_full();
	}
	op->verbose=verbose;

}

void solver::create_full() {

	matrix opi;
	solver_elem *p;
	int **nr=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth;
	int i,j,k,l,ivar,ieq,n,nn,i0,j0,set,nn2,kk,kk0,kk1,jj;
	
	for(n=0;n<nb;n++) {
		nn=0;
		for(ivar=0;ivar<nv;ivar++) nn+=nr[n][ivar]*nth[n][ivar];
		opi=zeros(nn,nn);
		j0=0;
		
		for(kk=0;kk<3;kk++) {
		for(ivar=0;ivar<nv;ivar++) {
			if(kk==0) {kk0=0;kk1=nbot[n][ivar];}
			if(kk==1) {kk0=nbot[n][ivar];kk1=nr[n][ivar]-ntop[n][ivar];}
			if(kk==2) {kk0=nr[n][ivar]-ntop[n][ivar];kk1=nr[n][ivar];}
			for(k=kk0;k<kk1;k++) {
				for(l=0;l<nth[n][ivar];l++) {
					i0=0;
					for(ieq=0;ieq<nv;ieq++) {
						for(i=0;i<nbot[n][ieq];i++) {
							for(j=0;j<nth[n][ieq];j++) {
								p=bc_bot2[n].eq[ieq][ivar];
								while(p!=NULL) {
									switch(p->type) {
									case 'd':
										if(k==0&&(j==l||nth[n][ivar]==1))
											opi(i0,j0)+=p->D(i,j);
										break;
									case 'l':
										if(j==l||nth[n][ivar]==1)
											opi(i0,j0)+=p->D(i,j)*p->L(i,k);
										break;
									case 'r':
										if(k==0)
											opi(i0,j0)+=p->D(i,j)*p->R(l,j);
										break;
									case 'f':
										opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j);
										break;
									case 'm':
										if(j==l||nth[n][ivar]==1)
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->I(k,l);
										break;
									case 's':
										if(k==0)
											opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(0,l);
										break;
									case 'g':
										opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j)*p->I(k,l);
									}
									p=p->next;	
								}
								i0++;
							}
						}
					}
					for(ieq=0;ieq<nv;ieq++) {
						for(i=nbot[n][ieq];i<nr[n][ieq]-ntop[n][ieq];i++) {
							for(j=0;j<nth[n][ieq];j++) {
								if(j==0&&bc_eq[n].eq_set(ieq)) {
									p=bc_eq[n].eq[ieq][ivar];
									jj=0;
								} else if(j==nth[n][ieq]-1&&bc_pol[n].eq_set(ieq)) {
									p=bc_pol[n].eq[ieq][ivar];
									jj=0;
								} else {
									p=block[n].eq[ieq][ivar];
									jj=j;
								}
								while(p!=NULL) {
									switch(p->type) {
									case 'd':
										if((i==k||nr[n][ivar]==1)&&(j==l||nth[n][ivar]==1))
											opi(i0,j0)+=p->D(i,jj);
										break;
									case 'l':
										if((j==l)||(nth[n][ivar]==1))
											opi(i0,j0)+=p->D(i,jj)*p->L(i,k);
										break;
									case 'r':
										if((i==k)||(nr[n][ivar]==1))
											opi(i0,j0)+=p->D(i,jj)*p->R(l,jj);
										break;
									case 'f':
										opi(i0,j0)+=p->D(i,jj)*p->L(i,k)*p->R(l,jj);
										break;
									case 'm':
										if((j==l)||(nth[n][ivar]==1))
											opi(i0,j0)+=p->D(i,jj)*p->L(i,k)*p->I(k,jj);
										break;
									case 's':
										if((i==k)||(nr[n][ivar]==1))
											opi(i0,j0)+=p->D(i,jj)*p->R(l,jj)*p->I(k,l);
										break;
									case 'g':
										opi(i0,j0)+=p->D(i,jj)*p->L(i,k)*p->R(l,jj)*p->I(k,l);
									}
									p=p->next;	
								}
								i0++;
							}
						}
					}
					for(ieq=0;ieq<nv;ieq++) {
						for(i=0;i<ntop[n][ieq];i++) {
							for(j=0;j<nth[n][ieq];j++) {
								p=bc_top1[n].eq[ieq][ivar];
								while(p!=NULL) {
									switch(p->type) {
									case 'd':
										if(k==(nr[n][ivar]-1)&&(j==l||nth[n][ivar]==1))
											opi(i0,j0)+=p->D(i,j);
										break;
									case 'l':
										if(j==l||nth[n][ivar]==1)
											opi(i0,j0)+=p->D(i,j)*p->L(i,k);
										break;
									case 'r':
										if(k==(nr[n][ivar]-1))
											opi(i0,j0)+=p->D(i,j)*p->R(l,j);
										break;
									case 'f':
										opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j);
										break;
									case 'm':
										if(j==l||nth[n][ivar]==1)
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->I(k,l);
										break;
									case 's':
										if(k==(nr[n][ivar]-1))
											opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(0,l);
										break;
									case 'g':
										opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j)*p->I(k,l);
									}
									p=p->next;	
								}
								i0++;
							}
						}
					}
					j0++;
				}
			}
		}
		}
		op->set_block(n,opi);
		
		if(n>0) {
			nn=0;
			for(ivar=0;ivar<nv;ivar++) nn+=nbot[n][ivar]*nth[n][ivar];
			nn2=0;
			for(ivar=0;ivar<nv;ivar++) nn2+=nr[n-1][ivar]*nth[n-1][ivar];
			opi=zeros(nn,nn2);
			j0=0;
			set=0;
			for(kk=0;kk<3;kk++) {
			for(ivar=0;ivar<nv;ivar++) {
				if(kk==0) {kk0=0;kk1=nbot[n-1][ivar];}
				if(kk==1) {kk0=nbot[n-1][ivar];kk1=nr[n-1][ivar]-ntop[n-1][ivar];}
				if(kk==2) {kk0=nr[n-1][ivar]-ntop[n-1][ivar];kk1=nr[n-1][ivar];}
				for(k=kk0;k<kk1;k++) {
					for(l=0;l<nth[n-1][ivar];l++) {
						i0=0;
						for(ieq=0;ieq<nv;ieq++) {
							for(i=0;i<nbot[n][ieq];i++) {
								for(j=0;j<nth[n][ieq];j++) {
									p=bc_bot1[n].eq[ieq][ivar];
									while(p!=NULL) {
										set=1;
										switch(p->type) {
										case 'd':
											if(k==(nr[n-1][ivar]-1)&&(j==l||nth[n-1][ivar]==1)) 
												opi(i0,j0)+=p->D(i,j);
											break;
										case 'l':
											if(j==l||nth[n-1][ivar]==1)
												opi(i0,j0)+=p->D(i,l)*p->L(i,k);
											break;
										case 'r':
											if(k==(nr[n-1][ivar]-1))
												opi(i0,j0)+=p->D(i,j)*p->R(l,j);
											break;
										case 'f':
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j);
											break;
										case 'm':
											if(j==l||nth[n-1][ivar]==1)
												opi(i0,j0)+=p->D(i,l)*p->L(i,k)*p->I(k,l);
											break;
										case 's':
											if(k==(nr[n-1][ivar]-1))
												opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(0,l);
											break;
										case 'g':
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j)*p->I(k,l);
										}
										p=p->next;	
									}
									i0++;
								}
							}
						}
						j0++;
					}
				}
			}
			}
			if(set) op->set_blockinf(n-1,opi);

		}
		if(n<nb-1) {
			nn=0;
			for(ivar=0;ivar<nv;ivar++) nn+=ntop[n][ivar]*nth[n][ivar];
			nn2=0;
			for(ivar=0;ivar<nv;ivar++) nn2+=nr[n+1][ivar]*nth[n+1][ivar];
			opi=zeros(nn,nn2);
			j0=0;
			set=0;
			for(kk=0;kk<3;kk++) {
			for(ivar=0;ivar<nv;ivar++) {
				if(kk==0) {kk0=0;kk1=nbot[n+1][ivar];}
				if(kk==1) {kk0=nbot[n+1][ivar];kk1=nr[n+1][ivar]-ntop[n+1][ivar];}
				if(kk==2) {kk0=nr[n+1][ivar]-ntop[n+1][ivar];kk1=nr[n+1][ivar];}
				for(k=kk0;k<kk1;k++) {
					for(l=0;l<nth[n+1][ivar];l++) {
						i0=0;
						for(ieq=0;ieq<nv;ieq++) {
							for(i=0;i<ntop[n][ieq];i++) {
								for(j=0;j<nth[n][ieq];j++) {
									p=bc_top2[n].eq[ieq][ivar];
									while(p!=NULL) {
										set=1;
										switch(p->type) {
										case 'd':
											if(k==0&&(j==l||nth[n+1][ivar]==1))
												opi(i0,j0)+=p->D(i,j);
											break;
										case 'l':
											if(j==l||nth[n+1][ivar]==1)
												opi(i0,j0)+=p->D(i,l)*p->L(i,k);
											break;
										case 'r':
											if(k==0)
												opi(i0,j0)+=p->D(i,j)*p->R(l,j);
											break;
										case 'f':
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j);
											break;
										case 'm':
											if(j==l||nth[n+1][ivar]==1)
												opi(i0,j0)+=p->D(i,l)*p->L(i,k)*p->I(k,l);
											break;
										case 's':
											if(k==0)
												opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(0,l);
											break;
										case 'g':
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j)*p->I(k,l);
										}
										p=p->next;	
									}
									i0++;
								}
							}
						}
						j0++;
					}
				}
			}
			}
			if(set) op->set_blocksup(n,opi);

		}
		
	}

}

void solver::wrap(const matrix *y,matrix *x) {

	int i,j,j0[nv],i0,**n=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth,&N=solver_N,nn,nnt;
	matrix z,q;

	for(j=0;j<nv;j++) {
		nn=0;
		nnt=0;
		for(i=0;i<nb;i++) {nn+=n[i][j];nnt=(nnt>nth[i][j])?nnt:nth[i][j];}
		if(nn) if(y[j].nrows()!=nn||y[j].ncols()!=nnt) {*x=zeros(N,1);return;}
	}
	z.dim(N,1);
	i0=0;
	for(j=0;j<nv;j++) j0[j]=0;
	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(nbot[i][j]) {
				q=y[j].block(j0[j],j0[j]+nbot[i][j]-1,0,nth[i][j]-1).transpose();
				q.redim(nbot[i][j]*nth[i][j],1);
				z.setblock(i0,i0+q.nrows()-1,0,0,q);
				i0+=q.nrows();
			}
		}
		for(j=0;j<nv;j++) {
			nn=n[i][j]-nbot[i][j]-ntop[i][j];
			if(nn) {
				q=y[j].block(j0[j]+nbot[i][j],j0[j]+nbot[i][j]+nn-1,0,nth[i][j]-1).transpose();
				q.redim(nn*nth[i][j],1);
				z.setblock(i0,i0+q.nrows()-1,0,0,q);
				i0+=q.nrows();
			}
		}
		for(j=0;j<nv;j++) {
			if(ntop[i][j]) {
				q=y[j].block(j0[j]+n[i][j]-ntop[i][j],j0[j]+n[i][j]-1,0,nth[i][j]-1).transpose();
				q.redim(ntop[i][j]*nth[i][j],1);
				z.setblock(i0,i0+q.nrows()-1,0,0,q);
				i0+=q.nrows();
			}
		}
		for(j=0;j<nv;j++) j0[j]+=n[i][j];
	}
	
	*x=z;
	
}


void solver::unwrap(matrix *y,const matrix *x) {

	int i,j,j0[nv],i0,**n=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth,&N=solver_N,nn,nnt;
	matrix q;

	for(j=0;j<nv;j++) {
		nn=0;
		nnt=0;
		for(i=0;i<nb;i++) {nn+=n[i][j];nnt=(nnt>nth[i][j])?nnt:nth[i][j];}
		y[j]=zeros(nn,nnt);
	}
	i0=0;
	for(j=0;j<nv;j++) j0[j]=0;
	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(nbot[i][j]) {
				q=x->block(i0,i0+nbot[i][j]*nth[i][j]-1,0,0);
				i0+=q.nrows();
				q.redim(nth[i][j],nbot[i][j]);
				y[j].setblock(j0[j],j0[j]+nbot[i][j]-1,0,nth[i][j]-1,q.transpose());
			}
		}
		for(j=0;j<nv;j++) {
			nn=n[i][j]-nbot[i][j]-ntop[i][j];
			if(nn) {
				q=x->block(i0,i0+nn*nth[i][j]-1,0,0);
				i0+=q.nrows();
				q.redim(nth[i][j],nn);
				y[j].setblock(j0[j]+nbot[i][j],j0[j]+nbot[i][j]+nn-1,0,nth[i][j]-1,q.transpose());
			}
		}
		for(j=0;j<nv;j++) {
			if(ntop[i][j]) {
				q=x->block(i0,i0+ntop[i][j]*nth[i][j]-1,0,0);
				i0+=q.nrows();
				q.redim(nth[i][j],ntop[i][j]);
				y[j].setblock(j0[j]+n[i][j]-ntop[i][j],j0[j]+n[i][j]-1,0,nth[i][j]-1,q.transpose());
			}
		}
		for(j=0;j<nv;j++) j0[j]+=n[i][j];
	}
	
}


void solver::calc_struct() {

	int &N=solver_N,i,j,**nr=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth;

	N=0;
	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(block[i].eq_set(j)) {
				nr[i][j]=block[i].nr[j];
				nth[i][j]=block[i].nth[j];
				N+=nr[i][j]*nth[i][j];
			} else {
				nr[i][j]=0;
				nth[i][j]=0;
				nbot[i][j]=0;
				ntop[i][j]=0;
				continue;
			}
			if(bc_bot2[i].eq_set(j)) nbot[i][j]=bc_bot2[i].nr[j];
			else nbot[i][j]=0;
			if(bc_top1[i].eq_set(j)) ntop[i][j]=bc_top1[i].nr[j];
			else ntop[i][j]=0;
		}
	}
	/*for(j=0;j<nv;j++) {
		printf("%d %s:\n",j,var[j]);
		for(i=0;i<nb;i++) {
			printf("\t%d: %dx%d (nbot=%d,ntop=%d)\n",i,nr[i][j],nth[i][j],nbot[i][j],ntop[i][j]);
		}
	}*/
}

