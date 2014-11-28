#include "config.h"
#include "solver.h"
#include <stdlib.h>
#include <cmath>
#include <time.h>
extern "C" {
#ifdef USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#else
#include <cblas.h>
#include <lapack.h>
#endif
}

solver_full::solver_full(int nblocks,int offcore) {
	
	nb=nblocks;
	N=new int[nb];
	m=new matrix[nb];
	c=new matrix[nb];
	r=new matrix[nb];
	if(nb>1) {
		msup=new matrix[nb-1];
		minf=new matrix[nb-1];
		blk_index.zero(nb-1,2);
	} else {
		msup=NULL;
		minf=NULL;
	}
	ipiv=new int*[nb];
	ipiv_flag.zero(nb,1);
	lu_flag=0;
	
	oc=offcore;
	
	if(oc) {
		sprintf(tempdir,"solver_full_oc_tmp_XXXXXX");
		mkdtemp(tempdir);
	}

}

solver_full::~solver_full() {
	
	int i;
	
	delete [] N;
	delete [] m;
	delete [] r;
	delete [] c;
	 
	for(i=0;i<nb;i++) if(ipiv_flag(i)) delete [] ipiv[i];
	delete [] ipiv;
	if(nb>1) {
		delete [] msup;
		delete [] minf;
	}
	if(oc) {
		char tmp[100];
		sprintf(tmp,"rm -r %s",tempdir);
		system(tmp);
	}
}


void solver_full::set_block(int iblock,matrix &a) {

	
	N[iblock]=a.nrows();
	if(!oc) m[iblock].swap(a);
	else write_block(iblock,a);
	lu_flag=0;
	if(iblock>0) blk_index(iblock-1,0)=0;
	if(iblock<nb-1) blk_index(iblock,1)=0;
	
}

void solver_full::set_blocksup(int iblock,matrix &a) {

	if(!oc) msup[iblock].swap(a);
	else write_blocksup(iblock,a);
	blk_index(iblock,1)=1;
	lu_flag=0;
}

void solver_full::set_blockinf(int iblock,matrix &a) {

	if(!oc) minf[iblock].swap(a);
	else write_blockinf(iblock,a);
	blk_index(iblock,0)=1;
	lu_flag=0;
}

/// \brief Solves the set of equation stored in the object given the right hand
/// side \p rhs.
///
/// Solving the equation is done with Crout's method: LU factorization followed
/// forward and backward substitution.
matrix solver_full::solve(const matrix &rhs) {

	matrix x;
	
	if(!lu_flag) lu_calc();

	x=rhs;
	fwd_subs(x);
	back_subs(x);

	return x;

}

void solver_full::lu_calc() {

	int i;
	matrix UU;

	if(oc) read_block(0);
	for(i=0;i<nb;i++) {
		if(i) {
			if(oc) {
				if(blk_index(i-1,0)) read_blockinf(i-1);
				if(blk_index(i-1,1)) read_blocksup(i-1);
			}
			if(blk_index(i-1,0)) {
				minf[i-1]=minf[i-1].transpose();
				solve_block(i-1,'T',minf[i-1]);
				minf[i-1]=minf[i-1].transpose();
			}
			if(oc) {
				write_block(i-1,m[i-1]);
				read_block(i);
			}
			if(blk_index(i-1,0)&&blk_index(i-1,1)) {
				UU=zeros(minf[i-1].nrows(),m[i].ncols());
				cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,minf[i-1].nrows(),m[i].ncols(),msup[i-1].nrows(),1,
					minf[i-1].data()+(minf[i-1].ncols()-msup[i-1].nrows())*minf[i-1].nrows(),minf[i-1].nrows(),
					msup[i-1].data(),msup[i-1].nrows(),0,
					UU.data(),minf[i-1].nrows());
				m[i].setblock(0,UU.nrows()-1,0,UU.ncols()-1,
					m[i].block(0,UU.nrows()-1,0,UU.ncols()-1)-UU);	
			}
			if(oc) {
				if(blk_index(i-1,0)) write_blockinf(i-1,minf[i-1]);
				if(blk_index(i-1,1)) flush_blocksup(i-1);
			}
		}

		lu_block(i);
	}
	if(oc) write_block(nb-1,m[nb-1]);
	lu_flag=1;	
}

void solver_full::lu_block(int i) {
	
	int n,info;
	double rowcnd,colcnd,amax;
	
	n=m[i].nrows();

	r[i]=ones(n,1);
	c[i]=ones(1,n);
	dgeequ_(&n,&n,m[i].data(),&n,r[i].data(),c[i].data(),&rowcnd,&colcnd,&amax,&info);

	m[i]*=r[i];
	m[i]*=c[i];

	c[i].redim(n,1);
	//char str[25];
	//sprintf(str,"matriz%d",i);
	//FILE *fp;int l;fp=fopen(str,"wb");l=m[i].nrows();fwrite(&l,sizeof(int),1,fp);l=m[i].ncols();fwrite(&l,sizeof(int),1,fp);m[i].write(fp,'b');fclose(fp);	

	char norm='1';
	double work[4*n],rcond,anorm;
	int iwork[n];
	
	if (verbose)
		anorm=dlange_(&norm,&n,&n,m[i].data(),&n,work);//printf("%e\n",anorm);
	
	if(ipiv_flag(i)) delete [] ipiv[i];
	ipiv[i]=new int[n];
	dgetrf_(&n,&n,m[i].data(),&n,ipiv[i],&info);
	ipiv_flag(i)=1;
	
	if(verbose) {
		dgecon_(&norm,&n,m[i].data(),&n,&anorm,&rcond,work,iwork,&info);
		printf("Block %d rcond: %e\n",i,rcond);
	}
	
}

void solver_full::solve_block(int i,char trans,matrix &x) {

	int n,nrhs,info=0;
	
	
	if(trans=='T') x=x*c[i];
	else x=x*r[i];
	
	n=m[i].nrows();nrhs=x.ncols();
	dgetrs_(&trans,&n,&nrhs,m[i].data(),&n,ipiv[i],x.data(),&n,&info);
	
	if(trans=='T') x=x*r[i];
	else x=x*c[i];

}



void solver_full::fwd_subs(matrix &x) {

	int i,j0,j1;
	
	j1=N[0];
	j0=0;
	for(i=1;i<nb;i++) {
		if(blk_index(i-1,0)) {
			if(oc) read_blockinf(i-1); 
			x.setblock(j1,j1+minf[i-1].nrows()-1,0,x.ncols()-1,
				x.block(j1,j1+minf[i-1].nrows()-1,0,x.ncols()-1)-
				( minf[i-1] , x.block(j0,j1-1,0,x.ncols()-1) )
				);
			if(oc) flush_blockinf(i-1);
			}
		j0=j1;
		j1+=N[i];
	}
	
}

void solver_full::back_subs(matrix &x) {
	
	int i,j0,j1;
	matrix xx;
		
	j1=0;
	for(i=0;i<nb;i++) j1+=N[i];
	for(i=nb-1;i>=0;i--) {
		j0=j1-N[i];
		if(i<nb-1)
			if(blk_index(i,1)) {
				if(oc) read_blocksup(i);
				x.setblock(j1-msup[i].nrows(),j1-1,0,x.ncols()-1,
					x.block(j1-msup[i].nrows(),j1-1,0,x.ncols()-1)-
					( msup[i] , x.block(j1,j1+N[i+1]-1,0,x.ncols()-1) )
					);
				if(oc) flush_blocksup(i);
			}
		xx=x.block(j0,j1-1,0,x.ncols()-1);
		if(oc) read_block(i);
		solve_block(i,'N',xx);
		if(oc) flush_block(i);
		x.setblock(j0,j1-1,0,x.ncols()-1,xx);
		j1=j0; 
	}
	
}


void solver_full::read_block(int i) {
	
	FILE *fp;
	char tmp[100];
	int nr,nc;
	
	sprintf(tmp,"%s/block_%d",tempdir,i);
	fp=fopen(tmp,"rb");
	fread(&nr,sizeof(int),1,fp);
	fread(&nc,sizeof(int),1,fp);
	m[i].read(nr,nc,fp,'b');
	fclose(fp);
	
}

void solver_full::write_block(int i,const matrix &a) {

	FILE *fp;
	char tmp[100];
	int nr,nc;
	
	sprintf(tmp,"%s/block_%d",tempdir,i);
	fp=fopen(tmp,"wb");
	nr=a.nrows();
	nc=a.ncols();
	fwrite(&nr,sizeof(int),1,fp);
	fwrite(&nc,sizeof(int),1,fp);
	a.write(fp,'b');
	fclose(fp);	
	m[i]=zeros(1,1);
	
}

void solver_full::flush_block(int i) {

	m[i]=zeros(1,1);
	
}

void solver_full::read_blockinf(int i) {
	
	FILE *fp;
	char tmp[100];
	int nr,nc;
	
	sprintf(tmp,"%s/blockinf_%d",tempdir,i);
	fp=fopen(tmp,"rb");
	fread(&nr,sizeof(int),1,fp);
	fread(&nc,sizeof(int),1,fp);
	minf[i].read(nr,nc,fp,'b');
	fclose(fp);
	
}

void solver_full::write_blockinf(int i,const matrix &a) {

	FILE *fp;
	char tmp[100];
	int nr,nc;
	
	sprintf(tmp,"%s/blockinf_%d",tempdir,i);
	fp=fopen(tmp,"wb");
	nr=a.nrows();
	nc=a.ncols();
	fwrite(&nr,sizeof(int),1,fp);
	fwrite(&nc,sizeof(int),1,fp);
	a.write(fp,'b');
	fclose(fp);	
	minf[i]=zeros(1,1);
	
}

void solver_full::flush_blockinf(int i) {

	minf[i]=zeros(1,1);
	
}

void solver_full::read_blocksup(int i) {
	
	FILE *fp;
	char tmp[100];
	int nr,nc;
	
	sprintf(tmp,"%s/blocksup_%d",tempdir,i);
	fp=fopen(tmp,"rb");
	fread(&nr,sizeof(int),1,fp);
	fread(&nc,sizeof(int),1,fp);
	msup[i].read(nr,nc,fp,'b');
	fclose(fp);
	
}

void solver_full::write_blocksup(int i,const matrix &a) {

	FILE *fp;
	char tmp[100];
	int nr,nc;
	
	sprintf(tmp,"%s/blocksup_%d",tempdir,i);
	fp=fopen(tmp,"wb");
	nr=a.nrows();
	nc=a.ncols();
	fwrite(&nr,sizeof(int),1,fp);
	fwrite(&nc,sizeof(int),1,fp);
	a.write(fp,'b');
	fclose(fp);	
	msup[i]=zeros(1,1);
	
}

void solver_full::flush_blocksup(int i) {

	msup[i]=zeros(1,1);
	
}



