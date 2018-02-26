#include "ester-config.h"
#include "utils.h"
#include "solver.h"
#include <string.h>
#include <stdlib.h>
#include <cmath>

//#define PERF_LOG

void solver::solver_block::init(int nvar) {

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

void solver::solver_block::destroy() {

	int i;

	reset();

	for(i=0;i<nv;i++) delete [] eq[i];
	delete [] eq;
	delete [] nr;
	delete [] nth;

}

void solver::solver_block::reset() {

	int i;

	for(i=0;i<nv;i++)
		reset(i);
}

void solver::solver_block::reset(int ieq) {

	solver_elem *p;
	solver_elem *p0;
	int j;

	for(j=0;j<nv;j++) {
		p=eq[ieq][j];
		while(p!=NULL) {
			p0=p;
			p=p0->next;
			delete p0;
		}
		eq[ieq][j]=NULL;
	}
	eq_set(ieq)=0;
	nr[ieq]=0;
	nth[ieq]=0;
}

void solver::solver_block::add(int ieq,int ivar,char type,const matrix *d,const matrix *l,const matrix *r,const matrix *i) {

	solver_elem *p,*pnew;

	if(!nr[ieq]&&!nth[ieq]) {
		nr[ieq]=d->nrows();
		nth[ieq]=d->ncols();
	}
	pnew=new solver_elem;
	if(eq[ieq][ivar]==NULL) {
		eq[ieq][ivar]=pnew;
	} else {
		p=eq[ieq][ivar];
		while(p->next!=NULL) p=p->next;
		p->next=pnew;
	}
	p=pnew;
	p->type=type;
	p->D=*d;
	if(type=='l'||type=='m'||type=='f'||type=='g') p->L=*l;
	if(type=='r'||type=='s'||type=='f'||type=='g') p->R=*r;
	if(type=='m'||type=='s'||type=='g') p->I=*i;

	p->next=NULL;

	eq_set(ieq)=1;

}

solver::solver() {
	initd=0;
	verbose=0;
	use_cgs=1;
	maxit_ref=10;
	maxit_cgs=10;
	rel_tol=1e-12;
	abs_tol=1e-20;
	def_nr=NULL;
	mult_fcn=NULL;
	mult_context=NULL;
	debug=0;
};


/// \brief Initialize the solver object for \p nblock blocks, \p nvar variables,
/// using a solver operator of type \p solver_type (either "full", "full-oc" or
/// "iter").
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
		for(int j=0;j<nv;j++) {
			var_nr[i][j]=0;
			var_nbot[i][j]=0;
			var_ntop[i][j]=0;
			var_nth[i][j]=0;
		}
	}
	var=new char*[nv];
	for(i=0;i<nv;i++) {
		var[i]=new char[32];var[i][0]='\0';
	}

	reg=zeros(1,nv);
	dep=zeros(1,nv);

	sol=new matrix[nv];
	rhs=new matrix[nv];

}



void solver::destroy() {

	int i;

	initd=0;
	if(op!=NULL) delete op;

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
	delete [] sol;
	delete [] rhs;

}

void solver::reset() {

	int i;

	sync=0;
	for(i=0;i<nb;i++)
		reset(i);

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

void solver::reset(int iblock,const char *eq_name) {
	reset(iblock,get_id(eq_name));
}

void solver::reset(const char *eq_name) {

	for(int i=0;i<nb;i++)
		reset(i,eq_name);
}

void solver::set_nr(int *nr) {
	int i;

	for(i=0;i<nb;i++) {
		def_nr[i]=nr[i];
	}
}

void solver::regvar(const char *var_name,int dependent) {

	int j;
	j=0;
	while (strlen(var[j])) {
		if(!strcmp(var[j],var_name)) {
			ester_err("ERROR: Can't register variable (already registered)");
			exit(1);
		}
		j++;
		if(j==nv) {
			ester_err("ERROR: Can't register variable (increase nvar)");
			exit(1);
		}
	}

	strncpy(var[j],var_name,31);
	dep(j)=dependent;
	reg(j)=1;

}


int solver::get_nvar() {

	return nv;

}

int solver::get_nblocks() {

	return nb;

}

int solver::get_id(const char *varn) {

	int i=0;

	while(strcmp(varn,var[i])||!reg(i)) {
		i++;
		if(i==nv) {
			ester_err("ERROR: Unknown variable %s", varn);
			exit(1);
		}
	}
	return i;

}

void solver::fill_void_blocks() {

	int i,j,n,nt;

	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(!reg(j)) continue;
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
	if(dep(ieq)) {
		ester_err("ERROR (solver):\n\t");
		ester_err("RHS not used in the definition of dependent variable \"%s\"",
                eqn);
        std::exit(EXIT_FAILURE);
	}

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
	return sol[ivar];
}

matrix_map solver::get_vars() {

	matrix_map map;
	for(int i=0;i<nv;i++)
		if(reg(i)) map[std::string(var[i])]=sol[i];

	return map;
}

matrix solver::get_rhs(int ieq) {

	return rhs[ieq];
}

matrix solver::get_var(int ivar) {

	return sol[ivar];
}

/// \brief Solves the set of equations stored in the operator object.
void solver::solve(int *info) {

    matrix y;
    matrix rhsw;
    int err,struct_changed=0;
    int i;

#ifdef PERF_LOG
    static tiempo tref,tcgs,tlu,tcreate,ttot;
    static int nref=0,ncgs=0,nlu=0,ncreate=0;
    ttot.start();
#endif
    if(info!=NULL) for(i=0;i<5;i++) info[i]=0;
    if(!sync) fill_void_blocks();
    if(verbose)
        printf("Checking structure...\n");
    err=check_struct();
    if(err&2) exit(1);
    if(err) struct_changed=1;
    if(verbose)
        printf("Substitution of dependent variables...\n");
    subst_dep();
    wrap(sol,&y);
    wrap(rhs,&rhsw);
    if(op!=NULL&&sync==0&&use_cgs&&struct_changed&&verbose)
        ester_warn("Operator structure has changed, skipping CGS iteration");
    if(op!=NULL&&sync==0&&use_cgs&&!struct_changed) {
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
            unwrap(&y,sol);
#ifdef PERF_LOG
            ttot.stop();
            printf("Total: %2.3fs \nLU(%d): %2.3fs (%2.3fs) cgs(%d): %2.3fs (%2.3fs)  ref(%d): %2.3fs (%2.3fs)  create(%d): %2.3fs (%2.3fs)\n",ttot.value(),nlu,tlu.value(),tlu.value()/nlu,ncgs,tcgs.value(),tcgs.value()/ncgs,nref,tref.value(),tref.value()/nref,ncreate,tcreate.value(),tcreate.value()/ncreate);
#endif
            if(verbose)
                printf("Solving dependent variables...\n");
            solve_dep();
            return;
        }
        if(verbose)
            printf("Not converged after %d iterations\n",maxit_cgs);
    }
    if(!sync) {
        if(verbose)
            printf("Creating matrix...\n");
        if(debug) printf("SOLVER: Debug mode is ON\n");
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
    if(!sync) {
        if(verbose)
            printf("LU Factorization...\n");
#ifdef PERF_LOG
        tlu.start();
#endif
    } else {
        if(verbose)
            printf("Solving...\n");
    }
    y=op->solve(rhsw);
#ifdef PERF_LOG
    if(!sync) {
        tlu.stop();
        nlu++;
    }
#endif
    sync=1;
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
    unwrap(&y,sol);
    if(verbose)
        printf("Solving dependent variables...\n");
    solve_dep();
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

	matrix *z,yy,zz,zz0;
    z = new matrix[nv];
	int i,j,n,j0[nv];
	solver_elem *p;
	int **nr=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth;

	for(i=0;i<nv;i++) {
		if(!reg(i)) continue;
		z[i]=zeros(y[i].nrows(),y[i].ncols());
		for(j=0;j<nv;j++) j0[j]=0;
		for(n=0;n<nb;n++) {
			if(block[n].eq_set(i)) {
				zz=zeros(nr[n][i],nth[n][i]);
				for(j=0;j<nv;j++) {
					if(!reg(j)) continue;
					if((p=block[n].eq[i][j])!=NULL)
						yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
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
						if(!reg(j)) continue;
						if((p=bc_pol[n].eq[i][j])!=NULL)
							yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
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
						if(!reg(j)) continue;
						if((p=bc_eq[n].eq[i][j])!=NULL)
							yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
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
						if(!reg(j)) continue;
						if((p=bc_bot2[n].eq[i][j])!=NULL)
							yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
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
						if(!reg(j)) continue;
						if((p=bc_bot1[n].eq[i][j])!=NULL)
							yy=y[j].block(j0[j]-nr[n-1][j],j0[j]-1,0,nth[n-1][j]-1);
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
						if(!reg(j)) continue;
						if((p=bc_top1[n].eq[i][j])!=NULL)
							yy=y[j].block(j0[j],j0[j]+nr[n][j]-1,0,nth[n][j]-1);
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
						if(!reg(j)) continue;
						if((p=bc_top2[n].eq[i][j])!=NULL)
							yy=y[j].block(j0[j]+nr[n][j],j0[j]+nr[n][j]+nr[n+1][j]-1,0,nth[n+1][j]-1);
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
			}
			for(j=0;j<nv;j++) j0[j]+=nr[n][j];
		}
	}
	for(i=0;i<nv;i++) y[i]=z[i];
    delete [] z;
}

int solver::cgs(const matrix &rhs,matrix &x,int maxit) {

	matrix r,r_,u,p,q,v,err_rel,err_abs;
    matrix *y = new matrix[nv];
	int k,error=0,fin;
	double s,a,s_1=.0;

	r=x;
	unwrap(&r,y);
	mult(y);
	wrap(y,&r);
	r=rhs-r;
	op->left_precond(r);
	r_=r;
	k=0;
	fin=0;
	while(!fin) {
		if(k>=maxit) {
			error=1;
			break;
		}
		s=sum(r_*r);
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
		op->right_precond(v);
		unwrap(&v,y);
		mult(y);
		wrap(y,&v);
		op->left_precond(v);
		a=s/sum(r_*v);
		q=u-a*v;
		v=u+q;
		op->right_precond(v);
		err_abs=abs(a*v);
		err_rel=abs(err_abs/(x+(x==0)));
		x=x+a*v;
		unwrap(&v,y);
		mult(y);
		wrap(y,&v);
		op->left_precond(v);
		r=r-a*v;
		s_1=s;
		fin=!exist((err_rel>rel_tol)&&(err_abs>abs_tol));
		if (verbose)
			printf("%d max. rel. error:%e max. abs. error:%e\n",k,max(err_rel*(err_abs>=abs_tol)),max(err_abs*(err_rel>=rel_tol)));
		k++;
	}
    delete [] y;
	if(error) {
		return -error;
	}

	return k;
}

/// \brief Creates the solver solver operator depending on the \p type of solver
/// configured (defined when calling solver::init).
void solver::create() {

	if(op!=NULL) {
		delete op;
		op=NULL;
	}

	if(!strcmp(type,"full")) {
		op=new solver_full(nb);
		create_full();
	}else if(!strcmp(type,"full-oc")) {
		op=new solver_full(nb,1);
		create_full();
	} else if(!strcmp(type,"iter")) {
		op=new solver_iter();
	} else {
		ester_err("Unknown solver type \"%s\"", type);
		exit(1);
	}
	op->verbose=verbose;

}

/// \brief Creates a solver operator of type "full".
///
/// Stores all equation's terms defined in the fields \p block, \p bc_bot1,
/// \p bc_bot2, \p bc_top1, \p bc_top2, \p bc_pol and \p bc_eq into
/// the solver operator \p op.
void solver::create_full() {

	matrix opi;
	solver_elem *p;
	int **nr=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth;
	int i,j,k,l,ivar,ieq,n,nn,i0,j0,set,nn2,kk,kk0=0,kk1=0,jj;

	if(verbose) printf("\t");
	for(n=0;n<nb;n++) {
		nn=0;
		for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) nn+=nr[n][ivar]*nth[n][ivar];
		opi.zero(nn,nn);
		j0=0;

		for(kk=0;kk<3;kk++) {
		for(ivar=0;ivar<nv;ivar++) {
			if(!reg(ivar)||dep(ivar)) continue;
			if(kk==0) {kk0=0;kk1=nbot[n][ivar];}
			if(kk==1) {kk0=nbot[n][ivar];kk1=nr[n][ivar]-ntop[n][ivar];}
			if(kk==2) {kk0=nr[n][ivar]-ntop[n][ivar];kk1=nr[n][ivar];}
			for(k=kk0;k<kk1;k++) {
                for(l=0;l<nth[n][ivar];l++) {
					i0=0;
					for(ieq=0;ieq<nv;ieq++) {
						if(!reg(ieq)||dep(ieq)) continue;
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
										if(j==l||nth[n][ivar]==1) {
											int lI;
											if(p->I.ncols()==1) lI=0;
											else lI=j;
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->I(k,lI);
										}
										break;
									case 's':
										if(k==0) {
											int kI;
											if(p->I.nrows()==1) kI=0;
											else kI=i;
											opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(kI,l);
										}
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
						if(!reg(ieq)||dep(ieq)) continue;
						for(i=nbot[n][ieq];i<nr[n][ieq]-ntop[n][ieq];i++) {
							for(j=0;j<nth[n][ieq];j++) {
								int j_check;
								if(j==0&&bc_eq[n].eq_set(ieq)) {
									p=bc_eq[n].eq[ieq][ivar];
									jj=0;
									j_check=0;
								} else if(j==nth[n][ieq]-1&&bc_pol[n].eq_set(ieq)) {
									p=bc_pol[n].eq[ieq][ivar];
									jj=0;
									j_check=nth[n][ivar]-1;
								} else {
									p=block[n].eq[ieq][ivar];
									jj=j;
									j_check=j;
								}
								while(p!=NULL) {
									switch(p->type) {
									case 'd':
										if((i==k||nr[n][ivar]==1)&&(j_check==l||nth[n][ivar]==1))
											opi(i0,j0)+=p->D(i,jj);
										break;
									case 'l':
										if((j_check==l)||(nth[n][ivar]==1))
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
										if((j_check==l)||(nth[n][ivar]==1)) {
											int lI;
											if(p->I.ncols()==1) lI=0;
											else lI=jj;
											opi(i0,j0)+=p->D(i,jj)*p->L(i,k)*p->I(k,lI);
										}
										break;
									case 's':
										if((i==k)||(nr[n][ivar]==1)) {
											int kI;
											if(p->I.nrows()==1) kI=0;
											else kI=i;
											opi(i0,j0)+=p->D(i,jj)*p->R(l,jj)*p->I(kI,l);
										}
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
						if(!reg(ieq)||dep(ieq)) continue;
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
										if(j==l||nth[n][ivar]==1) {
											int lI;
											if(p->I.ncols()==1) lI=0;
											else lI=j;
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->I(k,lI);
										}
										break;
									case 's':
										if(k==(nr[n][ivar]-1)) {
											int kI;
											if(p->I.nrows()==1) kI=0;
											else kI=i;
											opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(kI,l);
										}
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
		if (debug) check_full(n,opi,0);
		((solver_full *) op)->set_block(n,opi);

		if(n>0) {
			nn=0;
			for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) nn+=nbot[n][ivar]*nth[n][ivar];
			nn2=0;
			for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) nn2+=nr[n-1][ivar]*nth[n-1][ivar];
			if(nn) opi.zero(nn,nn2);
			j0=0;
			set=0;
			for(kk=0;kk<3;kk++) {
			for(ivar=0;ivar<nv;ivar++) {
				if(!reg(ivar)||dep(ivar)) continue;
				if(kk==0) {kk0=0;kk1=nbot[n-1][ivar];}
				if(kk==1) {kk0=nbot[n-1][ivar];kk1=nr[n-1][ivar]-ntop[n-1][ivar];}
				if(kk==2) {kk0=nr[n-1][ivar]-ntop[n-1][ivar];kk1=nr[n-1][ivar];}
				for(k=kk0;k<kk1;k++) {
					for(l=0;l<nth[n-1][ivar];l++) {
						i0=0;
						for(ieq=0;ieq<nv;ieq++) {
							if(!reg(ieq)||dep(ieq)) continue;
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
												opi(i0,j0)+=p->D(i,j)*p->L(i,k);
											break;
										case 'r':
											if(k==(nr[n-1][ivar]-1))
												opi(i0,j0)+=p->D(i,j)*p->R(l,j);
											break;
										case 'f':
											opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->R(l,j);
											break;
										case 'm':
											if(j==l||nth[n-1][ivar]==1) {
												int lI;
												if(p->I.ncols()==1) lI=0;
												else lI=j;
												opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->I(k,lI);
											}
											break;
										case 's':
											if(k==(nr[n-1][ivar]-1)) {
												int kI;
												if(p->I.nrows()==1) kI=0;
												else kI=i;
												opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(kI,l);
											}
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
			if (debug&set) check_full(n,opi,-1);
			if(set) ((solver_full *) op)->set_blockinf(n-1,opi);

		}
		if(n<nb-1) {
			nn=0;
			for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) nn+=ntop[n][ivar]*nth[n][ivar];
			nn2=0;
			for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) nn2+=nr[n+1][ivar]*nth[n+1][ivar];
			if(nn) opi.zero(nn,nn2);
			j0=0;
			set=0;
			for(kk=0;kk<3;kk++) {
			for(ivar=0;ivar<nv;ivar++) {
				if(!reg(ivar)||dep(ivar)) continue;
				if(kk==0) {kk0=0;kk1=nbot[n+1][ivar];}
				if(kk==1) {kk0=nbot[n+1][ivar];kk1=nr[n+1][ivar]-ntop[n+1][ivar];}
				if(kk==2) {kk0=nr[n+1][ivar]-ntop[n+1][ivar];kk1=nr[n+1][ivar];}
				for(k=kk0;k<kk1;k++) {
					for(l=0;l<nth[n+1][ivar];l++) {
						i0=0;
						for(ieq=0;ieq<nv;ieq++) {
							if(!reg(ieq)||dep(ieq)) continue;
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
											if(j==l||nth[n+1][ivar]==1) {
												int lI;
												if(p->I.ncols()==1) lI=0;
												else lI=j;
												opi(i0,j0)+=p->D(i,j)*p->L(i,k)*p->I(k,lI);
											}
											break;
										case 's':
											if(k==0) {
												int kI;
												if(p->I.nrows()==1) kI=0;
												else kI=i;
												opi(i0,j0)+=p->D(i,j)*p->R(l,j)*p->I(kI,l);
											}
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
			if(debug&set) check_full(n,opi,1);
			if(set) ((solver_full *) op)->set_blocksup(n,opi);
		}
		if(verbose) {printf("#");fflush(stdout);}
	}
	if(verbose) printf("\n");

}

void solver::check_full(int n, const matrix &opi, int pos) {


	int ivar,iblock,j0=0,nj=0,i,j,i0=0,ni=0;
	matrix x;
    matrix *y = new matrix[nv];

	for(iblock=0;iblock<n+pos;iblock++)
		for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) j0+=var_nr[iblock][ivar]*var_nth[iblock][ivar];
	for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) nj+=var_nr[n+pos][ivar]*var_nth[n+pos][ivar];
	if(!pos) {
		ni=nj;
		i0=j0;
	} else if(pos==-1) {
		for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) ni+=var_nbot[n][ivar]*var_nth[n][ivar];
		i0=j0+nj;
	} else if(pos==1) {
		for(ivar=0;ivar<nv;ivar++) if(reg(ivar)&&!dep(ivar)) ni+=var_ntop[n][ivar]*var_nth[n][ivar];
		i0=j0-ni;
	}
	for(j=0;j<nj;j++) {
		x=zeros(solver_N,1);
		x(j0+j)=1;
		unwrap(&x,y);
		mult_op(y);
		wrap(y,&x);
		for(i=0;i<ni;i++) {
			if(fabs((x(i+i0)-opi(i,j))/opi(i,j))>1e-8)
				if(fabs(x(i+i0)-opi(i,j))>1e-12) {
					printf("\nError in operator");
					if(pos==-1) printf(" (inf) ");
					if(pos==1) printf(" (sup) ");
					printf(":\n\tblock=%d,row=%d,col=%d (%e,%e)\n",
						n,i,j,opi(i,j),x(i+i0));
					exit(1);
				}
		}
	}
    delete [] y;
}

/// \brief Stores distributed (ie. split into several blocks) matrix \p y into
/// the single vector \p x
void solver::wrap(const matrix *y,matrix *x) {

	int i,j,j0[nv],i0,**n=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth,&N=solver_N,nn,nnt;
	matrix z,q;

	for(j=0;j<nv;j++) {
		if(!reg(j)||dep(j)) continue;
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
			if(!reg(j)||dep(j)) continue;
			if(nbot[i][j]) {
				q=y[j].block(j0[j],j0[j]+nbot[i][j]-1,0,nth[i][j]-1).transpose();
				q.redim(nbot[i][j]*nth[i][j],1);
				z.setblock(i0,i0+q.nrows()-1,0,0,q);
				i0+=q.nrows();
			}
		}
		for(j=0;j<nv;j++) {
			if(!reg(j)||dep(j)) continue;
			nn=n[i][j]-nbot[i][j]-ntop[i][j];
			if(nn) {
				q=y[j].block(j0[j]+nbot[i][j],j0[j]+nbot[i][j]+nn-1,0,nth[i][j]-1).transpose();
				q.redim(nn*nth[i][j],1);
				z.setblock(i0,i0+q.nrows()-1,0,0,q);
				i0+=q.nrows();
			}
		}
		for(j=0;j<nv;j++) {
			if(!reg(j)||dep(j)) continue;
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


/// \brief Performs the inverse operation of the \p wrap function: restore the
/// single vector \p x into the blocked vector y.
void solver::unwrap(const matrix *x,matrix *y) {

	int i,j,j0[nv],i0,**n=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth,nn,nnt;
	matrix q;
    bool nan = false;

	for(j=0;j<nv;j++) {
		if(!reg(j)) continue;
		nn=0;
		nnt=0;
		for(i=0;i<nb;i++) {nn+=n[i][j];nnt=(nnt>nth[i][j])?nnt:nth[i][j];}
		y[j]=zeros(nn,nnt);
	}
	i0=0;
	for(j=0;j<nv;j++) j0[j]=0;
	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(!reg(j)||dep(j)) continue;
			if(nbot[i][j]) {
				q=x->block(i0,i0+nbot[i][j]*nth[i][j]-1,0,0);
				i0+=q.nrows();
				q.redim(nth[i][j],nbot[i][j]);
				y[j].setblock(j0[j],j0[j]+nbot[i][j]-1,0,nth[i][j]-1,q.transpose());
                if (std::isnan(max(abs(q)))) {
                    LOGE("var %8s (%2dx%2d), block %d (bc_bot): NaN\n", var[j],
                            nth[i][j], nbot[i][j], i);
                    nan = true;
                }
			}
		}
		for(j=0;j<nv;j++) {
			if(!reg(j)||dep(j)) continue;
			nn=n[i][j]-nbot[i][j]-ntop[i][j];
			if(nn) {
				q=x->block(i0,i0+nn*nth[i][j]-1,0,0);
				i0+=q.nrows();
				q.redim(nth[i][j],nn);
				y[j].setblock(j0[j]+nbot[i][j],j0[j]+nbot[i][j]+nn-1,0,nth[i][j]-1,q.transpose());
                if (std::isnan(max(abs(q)))) {
                    LOGE("var %8s (%2dx%2d), block %d (eq): NaN\n", var[j],
                            nth[i][j], nn, i);
                    nan = true;
                }
			}
		}
		for(j=0;j<nv;j++) {
			if(!reg(j)||dep(j)) continue;
			if(ntop[i][j]) {
				q=x->block(i0,i0+ntop[i][j]*nth[i][j]-1,0,0);
				i0+=q.nrows();
				q.redim(nth[i][j],ntop[i][j]);
				y[j].setblock(j0[j]+n[i][j]-ntop[i][j],j0[j]+n[i][j]-1,0,nth[i][j]-1,q.transpose());
                if (std::isnan(max(abs(q)))) {
                    LOGE("var %8s (%2dx%2d), block %d (bc_top): NaN\n", var[j],
                            nth[i][j], ntop[i][j], i);
                    nan = true;
                }
			}
		}
		for(j=0;j<nv;j++) j0[j]+=n[i][j];
	}
    if (nan) ester_err("NaN values found in unwrap");
}



void solver::add(const char *eqn, const char *varn,const char *block_type,char type,const matrix *d,const matrix_block_diag *l,const matrix *r,const matrix *i) {

	int k,j0,j1;
	const matrix *ll;
	matrix *dd,*ii;
	char err_msg[256];
	int error=0;

	sync=0;
	dd=new matrix;
	ll=NULL;
	if(type=='m'||type=='s'||type=='g') ii=new matrix;
	else ii=NULL;
	j0=0;
	if(strcmp(block_type,"block")&&strcmp(block_type,"bc_eq")&&strcmp(block_type,"bc_pol")) {
		ester_err("(solver::add): Invalid block_type %s", block_type);
		exit(1);
	}

	for(k=0;k<nb;k++) {
		j1=j0+def_nr[k]-1;
		if(j1>=d->nrows()) {
			sprintf(err_msg,"Matrix D has incorrect size");
			error=1;
			break;
		}
		*dd=d->block(j0,j1,0,d->ncols()-1);
		if(type=='l'||type=='m'||type=='f'||type=='g') {
			if(k>=l->nblocks()) {
				sprintf(err_msg,"Matrix L has incorrect number of blocks");
				error=1;
				break;
			}
			ll=&(l->block(k));
		}
		if(type=='m'||type=='s'||type=='g') {
			if(j1>=i->nrows()) {
				sprintf(err_msg,"Matrix I has incorrect size");
				error=1;
				break;
			}
			*ii=i->block(j0,j1,0,i->ncols()-1);
		} else ii=NULL;
		add(k,eqn,varn,block_type,type,dd,ll,r,ii);
		j0=j1+1;
		if(j0==d->nrows()) {
			if(type=='l'||type=='m'||type=='f'||type=='g') {
				if(l->nblocks()>k+1) {
					sprintf(err_msg,"Matrix L has incorrect number of blocks");
					error=1;
					break;
				}
			}
			if(type=='m'||type=='s'||type=='g') {
				if(i->nrows()>j0) {
					sprintf(err_msg,"Matrix I has incorrect size");
					error=1;
					break;
				}
			}
			break;
		} else if(k==nb-1) {
			sprintf(err_msg,"Matrix D has incorrect size");
			error=1;
		}
	}

	if(error) {
		ester_err("ERROR (solver):\n\t%s\n\tin eq \"%s\", var \"%s\"",
                err_msg,eqn,varn);
		switch(type) {
			case 'd':
				ester_err(" (type: d)");
				break;
			case 'l':
				ester_err(" (type: l)");
				break;
			case 'r':
				ester_err(" (type: r)");
				break;
			case 'f':
				ester_err(" (type: lr)");
				break;
			case 'm':
				ester_err(" (type: li)");
				break;
			case 's':
				ester_err(" (type: ri)");
				break;
			case 'g':
				ester_err(" (type: lri)");
		}
		exit(1);
	}
	delete dd;
	if(type=='m'||type=='s'||type=='g') delete ii;

}

void solver::add(int iblock,const char *eqn, const char *varn,const char *block_type,char type,const matrix *d,const matrix *l,const matrix *r,const matrix *i) {

	solver_block *bb;

	int ieq,ivar;
	ieq=get_id(eqn);
	ivar=get_id(varn);

	if(dep(ieq)&&type!='d') {
		fprintf(stderr,"ERROR (solver):\n\t");
		fprintf(stderr,"Only type D terms are allowed in the definition of dependent variables\n");
		fprintf(stderr,"\tin block %d, eq \"%s\", var \"%s\"\n",iblock,eqn,varn);
		exit(1);
	}
	if(dep(ieq)&&strcmp(block_type,"block")) {
		fprintf(stderr,"ERROR (solver):\n\t");
		fprintf(stderr,"No boundary conditions are allowed in the definition of dependent variables\n");
		fprintf(stderr,"\tin block %d, eq \"%s\", var \"%s\"\n",iblock,eqn,varn);
		exit(1);
	}
	if(iblock==0&&!strcmp(block_type,"bc_bot1")) {
		fprintf(stderr,"ERROR (solver):\n\t");
		fprintf(stderr,"\"bc_bot1\" terms are not allowed in first domain\n");
		fprintf(stderr,"\tin block %d, eq \"%s\", var \"%s\"\n",iblock,eqn,varn);
		exit(1);
	}
	if(iblock==nb-1&&!strcmp(block_type,"bc_top2")) {
		fprintf(stderr,"ERROR (solver):\n\t");
		fprintf(stderr,"\"bc_top2\" terms are not allowed in last domain\n");
		fprintf(stderr,"\tin block %d, eq \"%s\", var \"%s\"\n",iblock,eqn,varn);
		exit(1);
	}

	sync=0;
	if(!strcmp(block_type,"block"))
		bb=block;
	else if(!strcmp(block_type,"bc_eq"))
		bb=bc_eq;
	else if(!strcmp(block_type,"bc_pol"))
		bb=bc_pol;
	else if(!strcmp(block_type,"bc_bot1"))
		bb=bc_bot1;
	else if(!strcmp(block_type,"bc_bot2"))
		bb=bc_bot2;
	else if(!strcmp(block_type,"bc_top1"))
		bb=bc_top1;
	else if(!strcmp(block_type,"bc_top2"))
		bb=bc_top2;
	else {
		fprintf(stderr,"ERROR (solver::add) : Unknown block_type %s\n",block_type);
		exit(1);
	}

	matrix *ll,L;
	ll=NULL;
	if(type=='l'||type=='f'||type=='m'||type=='g') {
		L=*l;
		ll=&L;
		if((!strcmp(block_type,"bc_bot1"))||(!strcmp(block_type,"bc_bot2"))
			||(!strcmp(block_type,"bc_top1"))||(!strcmp(block_type,"bc_top2"))) {
			if(l->nrows()==1||d->nrows()>1) L=L*ones(d->nrows(),1);
		}
	}

	bb[iblock].add(ieq,ivar,type,d,ll,r,i);

}


int solver::check_struct() {

	int &N=solver_N,i,j,**nr=var_nr,**nbot=var_nbot,**ntop=var_ntop,**nth=var_nth;
	int changed=0;

	N=0;
	for(i=0;i<nb;i++) {
		for(j=0;j<nv;j++) {
			if(block[i].eq_set(j)) {
				if(nr[i][j]!=block[i].nr[j]) changed=1;
				nr[i][j]=block[i].nr[j];
				if(nth[i][j]!=block[i].nth[j]) changed=1;
				nth[i][j]=block[i].nth[j];
				if(!dep(j)) N+=nr[i][j]*nth[i][j];
			} else {
				if(nr[i][j]||nth[i][j]||nbot[i][j]||ntop[i][j]) changed=1;
				nr[i][j]=0;
				nth[i][j]=0;
				nbot[i][j]=0;
				ntop[i][j]=0;
				continue;
			}
			if(bc_bot2[i].eq_set(j)) {
				if(nbot[i][j]!=bc_bot2[i].nr[j]) changed=1;
				nbot[i][j]=bc_bot2[i].nr[j];
			} else {
				if(nbot[i][j]!=0) changed=1;
				nbot[i][j]=0;
			}
			if(bc_top1[i].eq_set(j)) {
				if(ntop[i][j]!=bc_top1[i].nr[j]) changed=1;
				ntop[i][j]=bc_top1[i].nr[j];
			} else {
				if(ntop[i][j]!=0) changed=1;
				ntop[i][j]=0;
			}
		}
	}

	if(debug) {
		for(j=0;j<nv;j++) {
			if(!reg(j)) continue;
			printf("%d %s",j,var[j]);
			if(dep(j)) printf(" (dep) ");
			printf(":\n");
			for(i=0;i<nb;i++) {
				printf("\t%2d: %dx%d (nbot=%d,ntop=%d)\n",i,nr[i][j],nth[i][j],nbot[i][j],ntop[i][j]);
			}
		}
	}

	int error=0;

	for(int n=0;n<nb;n++) {
		for(int i=0;i<nv;i++) {
			if(!reg(i)) continue;
			if(block[n].eq_set(i))
				for(int j=0;j<nv;j++)
					error=error||check_struct_block(n,i,j);
			if(bc_pol[n].eq_set(i))
				for(int j=0;j<nv;j++)
					error=error||check_struct_bc_th(n,i,j,"pol");
			if(bc_eq[n].eq_set(i))
				for(int j=0;j<nv;j++)
					error=error||check_struct_bc_th(n,i,j,"eq");
			if(n)
				if(bc_bot1[n].eq_set(i))
					for(int j=0;j<nv;j++)
						error=error||check_struct_bc(n,i,j,"bot1");
			if(bc_bot2[n].eq_set(i))
				for(int j=0;j<nv;j++)
					error=error||check_struct_bc(n,i,j,"bot2");
			if(bc_top1[n].eq_set(i))
				for(int j=0;j<nv;j++)
					error=error||check_struct_bc(n,i,j,"top1");
			if(n<nb-1)
				if(bc_top2[n].eq_set(i))
					for(int j=0;j<nv;j++)
						error=error||check_struct_bc(n,i,j,"top2");
		}
	}

	int Nr,Nt;
	for(int i=0;i<nv;i++) {
		if (dep(i)||!reg(i)) continue;
		if(strlen(var[i])) {
			Nr=0;Nt=0;
			for(int n=0;n<nb;n++) {
				Nr+=nr[n][i];
				Nt=nth[n][i]>Nt?nth[n][i]:Nt;
			}
			if(rhs[i].nrows()!=Nr||rhs[i].ncols()!=Nt) {
				fprintf(stderr,"ERROR (solver):\n\tRHS for var \"%s\" has not the correct size\n",var[i]);
				fprintf(stderr,"\tIt is (%d,%d) and should be (%d,%d)\n",rhs[i].nrows(),rhs[i].ncols(),Nr,Nt);
				error=1;
			}
		}
	}

	return changed+2*error;

}

int solver::check_struct_block(int n,int i,int j) {

	int error=0,n1,m1,n2,m2;
	char err_msg[256];
	solver_elem *p;

	p=block[n].eq[i][j];
	n1=var_nr[n][i];m1=var_nth[n][i];
	n2=var_nr[n][j];m2=var_nth[n][j];


	while(p!=NULL) {

		if(!n2||!m2) {
			sprintf(err_msg,"Variable does not exist in this block");
			check_struct_error(err_msg,n,i,j,p);
			error=1;
		}

		if(p->D.nrows()!=n1||p->D.ncols()!=m1) {
			sprintf(err_msg,"Matrix D has incorrect size");
			check_struct_error(err_msg,n,i,j,p);
			error=1;
		}

		if(p->type=='l'||p->type=='f'||p->type=='m'||p->type=='g') {
			if(p->L.nrows()!=n1||p->L.ncols()!=n2) {
				sprintf(err_msg,"Matrix L has incorrect size");
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='r'||p->type=='f'||p->type=='s'||p->type=='g') {
			if(p->R.ncols()!=m1||p->R.nrows()!=m2) {
				sprintf(err_msg,"Matrix R has incorrect size");
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='m'||p->type=='s'||p->type=='g') {
			bool cond = false;
			if(p->type=='g') cond=p->I.nrows()==n2&&p->I.ncols()==m2;
			if(p->type=='m') cond=p->I.nrows()==n2&&(p->I.ncols()==m1||p->I.ncols()==m2);
			if(p->type=='s') cond=(p->I.nrows()==n1||p->I.nrows()==n2)&&p->I.ncols()==m2;
			if(!cond) {
				sprintf(err_msg,"Matrix I has incorrect size");
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='d') {
			bool cond;
			cond= (n2==n1)||(n2==1);
			cond= cond && ((m2==m1)||(m2==1));
			if(!cond) {
				sprintf(err_msg,"Term incompatible with size of variables");
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}
		if(p->type=='l'||p->type=='m') {
			bool cond;
			cond= (m2==m1)||(m2==1);
			if(!cond) {
				sprintf(err_msg,"Term incompatible with size of variables");
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}
		if(p->type=='r'||p->type=='s') {
			bool cond;
			cond= (n2==n1)||(n2==1);
			if(!cond) {
				sprintf(err_msg,"Term incompatible with size of variables");
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}
		p=p->next;
	}
	return error;

}

int solver::check_struct_bc_th(int n,int i,int j,const char *bctype) {

	int error=0,n1,n2,m2;
	char err_msg[256];
	solver_elem *p;

	if(!strcmp(bctype,"pol")) p=bc_pol[n].eq[i][j];
	else if(!strcmp(bctype,"eq")) p=bc_eq[n].eq[i][j];
	else return 1;

	n1=var_nr[n][i];
	n2=var_nr[n][j];m2=var_nth[n][j];

	while(p!=NULL) {

		if(!n2||!m2) {
			sprintf(err_msg,"Variable does not exist in this block (bc_%s)",bctype);
			check_struct_error(err_msg,n,i,j,p);
			error=1;
		}

		if(p->D.nrows()!=n1||p->D.ncols()!=1) {
			sprintf(err_msg,"Matrix D has incorrect size (bc_%s)",bctype);
			check_struct_error(err_msg,n,i,j,p);
			error=1;
		}

		if(p->type=='l'||p->type=='f'||p->type=='m'||p->type=='g') {
			if(p->L.nrows()!=n1||p->L.ncols()!=n2) {
				sprintf(err_msg,"Matrix L has incorrect size (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='r'||p->type=='f'||p->type=='s'||p->type=='g') {
			if(p->R.ncols()!=1||p->R.nrows()!=m2) {
				sprintf(err_msg,"Matrix R has incorrect size (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='m'||p->type=='s'||p->type=='g') {
			bool cond = false;
			if(p->type=='g') cond=p->I.nrows()==n2&&p->I.ncols()==m2;
			if(p->type=='m') cond=p->I.nrows()==n2&&p->I.ncols()==1;
			if(p->type=='s') cond=(p->I.nrows()==n1||p->I.nrows()==n2)&&p->I.ncols()==m2;
			if(!cond) {
				sprintf(err_msg,"Matrix I has incorrect size (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='d'||p->type=='r'||p->type=='s') {
			bool cond;
			cond= (n2==n1)||(n2==1);
			if(!cond) {
				sprintf(err_msg,"Term incompatible with size of variables (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		p=p->next;
	}
	return error;

}


int solver::check_struct_bc(int n,int i,int j,const char *bctype) {

	int error=0,nc,m1,n2,m2;
	char err_msg[256];
	solver_elem *p;

	if(!strcmp(bctype,"bot1")) {
		p=bc_bot1[n].eq[i][j];
		nc=var_nbot[n][i];m1=var_nth[n][i];
		n2=var_nr[n-1][j];m2=var_nth[n-1][j];
	} else if(!strcmp(bctype,"bot2")) {
		p=bc_bot2[n].eq[i][j];
		nc=var_nbot[n][i];m1=var_nth[n][i];
		n2=var_nr[n][j];m2=var_nth[n][j];
	} else if(!strcmp(bctype,"top1")) {
		p=bc_top1[n].eq[i][j];
		nc=var_ntop[n][i];m1=var_nth[n][i];
		n2=var_nr[n][j];m2=var_nth[n][j];
	} else if(!strcmp(bctype,"top2")) {
		p=bc_top2[n].eq[i][j];
		nc=var_ntop[n][i];m1=var_nth[n][i];
		n2=var_nr[n+1][j];m2=var_nth[n+1][j];
	} else return 1;

	while(p!=NULL) {

		if(!n2||!m2) {
			sprintf(err_msg,"Variable does not exist in this block (bc_%s)",bctype);
			check_struct_error(err_msg,n,i,j,p);
			error=1;
		}

		if(p->D.nrows()!=nc||p->D.ncols()!=m1) {
			sprintf(err_msg,"Matrix D has incorrect size (bc_%s)",bctype);
			check_struct_error(err_msg,n,i,j,p);
			error=1;
		}

		if(p->type=='l'||p->type=='f'||p->type=='m'||p->type=='g') {
			if(p->L.nrows()!=nc||p->L.ncols()!=n2) {
				sprintf(err_msg,"Matrix L has incorrect size (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='r'||p->type=='f'||p->type=='s'||p->type=='g') {
			if(p->R.ncols()!=m1||p->R.nrows()!=m2) {
				sprintf(err_msg,"Matrix R has incorrect size (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='m'||p->type=='s'||p->type=='g') {
			bool cond = false;
			if(p->type=='g') cond=p->I.nrows()==n2&&p->I.ncols()==m2;
			if(p->type=='m') cond=p->I.nrows()==n2&&(p->I.ncols()==m1||p->I.ncols()==m2);
			if(p->type=='s') cond=(p->I.nrows()==nc||p->I.nrows()==1)&&p->I.ncols()==m2;
			if(!cond) {
				sprintf(err_msg,"Matrix I has incorrect size (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		if(p->type=='d'||p->type=='l'||p->type=='m') {
			bool cond;
			cond= (m2==m1)||(m2==1);
			if(!cond) {
				sprintf(err_msg,"Term incompatible with size of variables (bc_%s)",bctype);
				check_struct_error(err_msg,n,i,j,p);
				error=1;
			}
		}

		p=p->next;
	}
	return error;

}


void solver::check_struct_error(const char *err_msg,int n,int i,int j,solver_elem *p) {

	fprintf(stderr,"ERROR (solver):\n\t%s\n\tin block %d, eq \"%s\", var \"%s\"",err_msg,n,var[i],var[j]);
	switch(p->type) {
		case 'd':
			fprintf(stderr," (type: d)\n");
			break;
		case 'l':
			fprintf(stderr," (type: l)\n");
			break;
		case 'r':
			fprintf(stderr," (type: r)\n");
			break;
		case 'f':
			fprintf(stderr," (type: lr)\n");
			break;
		case 'm':
			fprintf(stderr," (type: li)\n");
			break;
		case 's':
			fprintf(stderr," (type: ri)\n");
			break;
		case 'g':
			fprintf(stderr," (type: lri)\n");
	}

}

/// \brief Substitutes dependent variables defined in the solver.
///
/// Dependent variables are defined as an expression of other variables of the
/// problem. This function is responsible for replacing these dependent
/// variables with their expression.
void solver::subst_dep() {

	solver_elem *p;

	for(int n=0;n<nb;n++) {
		int substd=1;
		while(substd==1) {
			substd=0;
			for(int i=0;i<nv;i++) {
				if(dep(i)&&block[n].eq[i][i]!=NULL) {
					fprintf(stderr,"ERROR (solver):\n\tFound loop in definition of dependent variable \"%s\"\n",var[i]);
					exit(1);
				}
			}
			for(int i=0;i<nv;i++) {
				if(block[n].eq_set(i)&&dep(i)) {
					for(int j=0;j<nv;j++) {
						p=block[n].eq[i][j];
						if(dep(j)&&p!=NULL) {
							substd=1;
							subst_dep_eq("block",block,n,i,j);
						}
					}
				}
			}
		}
	}

	solver_block *bb = NULL;
	char block_type[8];

	for(int k=0;k<7;k++) {
		switch(k) {
			case 0:
				bb=block;strcpy(block_type,"block");break;
			case 1:
				bb=bc_eq;strcpy(block_type,"bc_eq");break;
			case 2:
				bb=bc_pol;strcpy(block_type,"bc_pol");break;
			case 3:
				bb=bc_bot1;strcpy(block_type,"bc_bot1");break;
			case 4:
				bb=bc_bot2;strcpy(block_type,"bc_bot2");break;
			case 5:
				bb=bc_top1;strcpy(block_type,"bc_top1");break;
			case 6:
				bb=bc_top2;strcpy(block_type,"bc_top2");
		}
		for(int n=0;n<nb;n++) {
			if(n==0&&!strcmp(block_type,"bc_bot1")) continue;
			if(n==nb-1&&!strcmp(block_type,"bc_top2")) continue;
			for(int i=0;i<nv;i++) {
				if(dep(i)) continue;
				if(bb[n].eq_set(i)) {
					for(int j=0;j<nv;j++) {
						p=bb[n].eq[i][j];
						if(dep(j)&&p!=NULL)
							subst_dep_eq(block_type,bb,n,i,j);
					}
				}
			}
		}
	}
}

void solver::subst_dep_eq(const char *block_type,solver_block *bb,int n,int i,int j) {

	solver_elem *p,*pdep,*p0;
	matrix d;
	int n2,m2,ndep;

	p=bb[n].eq[i][j];
	ndep=n;
	if(!strcmp(block_type,"bc_bot1")) ndep--;
	if(!strcmp(block_type,"bc_top2")) ndep++;

	while(p!=NULL) {
		for(int k=0;k<nv;k++) {
			pdep=block[ndep].eq[j][k];
			while(pdep!=NULL) {
				n2=var_nr[ndep][k];m2=var_nth[ndep][k];
				d=pdep->D;
				if(!strcmp(block_type,"bc_eq"))
					if(p->type=='d'||p->type=='l'||p->type=='m') {
						m2=1;
						d=d.col(0);
					}
				if(!strcmp(block_type,"bc_pol"))
					if(p->type=='d'||p->type=='l'||p->type=='m') {
						m2=1;
						d=d.col(d.ncols()-1);
					}
				if(!strcmp(block_type,"bc_bot2")||!strcmp(block_type,"bc_top2"))
					if(p->type=='d'||p->type=='r'||p->type=='s') {
						n2=1;
						d=d.row(0);
					}
				if(!strcmp(block_type,"bc_top1")||!strcmp(block_type,"bc_bot1"))
					if(p->type=='d'||p->type=='r'||p->type=='s') {
						n2=1;
						d=d.row(d.nrows()-1);
					}
				subst_dep_elem(i,k,&bb[n],p,d,n2,m2);
				pdep=pdep->next;
			}
		}
		p0=p;
		p=p->next;
		delete p0;
	}

	bb[n].eq[i][j]=NULL;

}

void solver::subst_dep_elem(int i,int k,solver_block *bb,solver_elem *p,const matrix &d,int n2,int m2) {

	int n1,m1;
	matrix D,L,R,I;
	char type_new = 'x';

	n1=d.nrows();
	m1=d.ncols();

	D=p->D;L=p->L;R=p->R;I=p->I;

    switch(p->type) {

        case 'd':
            type_new='d';
            D=D*d;
            break;
        case 'l':
            if(n1==n2) {
                type_new='m';
                I=d;
            } else {
                type_new='d';
                D=D*(L,d);
            }
            break;
        case 'r':
            if(m1==m2) {
                type_new='s';
                I=d;
            } else {
                type_new='d';
                D=D*(d,R);
            }
            break;
        case 'f':
            if(n1==n2&&m1==m2) {
                type_new='g';
                I=d;
            } else if(n1==n2) {
                type_new='l';
                D=D*(d,R);
            } else if(m1==m2) {
                type_new='r';
                D=D*(L,d);
            } else {
                type_new='d';
                D=D*(L,d,R);
            }
            break;
        case 'm':
            if(n1==n2) {
                type_new='m';
                I=I*d;
            } else {
                type_new='d';
                D=D*(L,I*d);
            }
            break;
        case 's':
            if(m1==m2) {
                type_new='s';
                I=I*d;
            } else {
                type_new='d';
                D=D*(I*d,R);
            }
            break;
        case 'g':
            if(n1==n2&&m1==m2) {
                type_new='g';
                I=I*d;
            } else if(n1==n2) {
                type_new='m';
                I=(I*d,R);
            } else if(m1==m2) {
                type_new='s';
                I=(L,I*d);
            } else {
                type_new='d';
                D=D*(L,I*d,R);
            }
            break;
        default:
            fprintf(stderr,
                    "Error in solver::subst_dep_elem: unknown type (%c)\n",
                    p->type);
            exit(EXIT_FAILURE);
    }

	bb->add(i,k,type_new,&D,&L,&R,&I);

}

/// \brief Solves dependent variables: compute in the \p sol field the values of
/// dependent variables (have to be called after solver::solve).
void solver::solve_dep() {

	matrix *x = new matrix[nv];

	for(int i=0;i<nv;i++) {
		x[i]=sol[i];
	}
	mult(x);
	for(int i=0;i<nv;i++) {
		if(dep(i))
			sol[i]=x[i];
	}
    delete [] x;
}
