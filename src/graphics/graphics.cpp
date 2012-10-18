#include"matrix.h"
#include"graphics.h"
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<math.h>

#if USE_PGPLOT==1
#include<cpgplot.h>

figure::figure(const char *device) {

	int i,cmax;
	float f;
	
	strcpy(dev,device);
	id=cpgopen(dev);
	cpgask(0);
	xlog=0;ylog=0;
	axis(0,0,0,0,0);
	hold_state=0;
	draw_state=1;
	draw_colorbar=0;
	axis_set=0;
	caxis_set=0;
	cpgslct(id);
	cpgqcol(&i,&cmax);
	for(i=16;i<=cmax;i++) {
		f=float(i-16)/(cmax-16);
		if(f<0.125)
			cpgscr(i,0,0,(f+0.125)/0.25);
		else if(f<0.375)
			cpgscr(i,0,(f-0.125)/0.25,1);
		else if(f<0.625)
			cpgscr(i,(f-0.375)/0.25,1,1-(f-0.375)/0.25);
		else if(f<0.875)
			cpgscr(i,1,1-(f-0.625)/0.25,0);
		else
			cpgscr(i,1-(f-0.875)/0.25,0,0);
	}

}
	
figure::~figure() {

	cpgslct(id);
	cpgclos();
}	


void figure::plot(const matrix &x,const matrix &y,const char *line) {

	int i,j,jx,jy,n,mx,my,marker=0,line_style=1,color=-1,c=0;
	matrix xx,yy;
	
	if(x.nrows()==1&&x.ncols()>1) {
		xx=x.transpose();
		yy=y.transpose();
	} else {
		xx=x;
		yy=y;
	}
	
	float xf[xx.nrows()],yf[yy.nrows()];
	
	n=xx.nrows();
	mx=xx.ncols();
	my=yy.ncols();
	if(my>1&&mx!=my) mx=1;
	if(!hold_state) {
		if(!axis_set) {
			x0=min(xx);x1=max(xx);y0=min(yy);y1=max(yy);just=0;
			if((float)x0==0&&(float)x1==0) {x0=x0-1;x1=x1+1;}
			else if((float)x0==(float)x1) {
				x0=x0-0.1*x0;
				x1=x1+0.1*x1;
			}
			else {
				x0=x0-(x1-x0)*0.1;
				x1=x1+(x1-x0)*0.1;
			}
			if((float)y0==0&&(float)y1==0) {y0=y0-1;y1=y1+1;}
			else if((float)y0==(float)y1) {
				y0=y0-0.1*y0;
				y1=y1+0.1*y1;
			}
			else {
				y0=y0-(y1-y0)*0.1;
				y1=y1+(y1-y0)*0.1;
			}
		}
	}
	axis_set=0;
	cpgslct(id);
	cpgsch(1);
	if(!hold_state) 
		cpgenv(x0,x1,y0,y1,just,10*xlog+20*ylog);
	for(i=0;i<strlen(line);i++) {
		switch(line[i]) {
			case '.': marker=1;break;
			case '+': marker=2;break;
			case '*': marker=3;break;
			case 'o': marker=4;break;
			case 'x': marker=5;
		}
	}
	if(marker) line_style=0;
	for(i=0;i<strlen(line);i++) {
		switch(line[i]) {
			case '-': line_style=1;break;
			case '=': line_style=2;break;
			case ';': line_style=3;break;
			case ':': line_style=4;
		}
	}
	for(i=0;i<strlen(line);i++) {
		switch(line[i]) {
			case 'w': color=0;break;
			case 'k': color=1;break;
			case 'r': color=2;break;
			case 'g': color=3;break;
			case 'b': color=4;break;
			case 'c': color=5;break;
			case 'm': color=6;break;
			case 'y': color=7;
		}
	}
	jx=0;jy=0;
	for(j=0;j<(mx>my?mx:my);j++) {
		for(i=0;i<n;i++) {
			xf[i]=xx(i+jx);
			yf[i]=yy(i+jy);
		}
		if(color==-1) {
			cpgsci(c+1);
			c++;
			c%=7;
		} else cpgsci(color);
		if(marker) cpgpt(n,xf,yf,marker);
		if(line_style) {
			cpgsls(line_style);
			cpgline(n,xf,yf);
		}
		if(mx>1) jx+=n;
		if(my>1) jy+=n;
	}
	cpgsls(1);
	cpgsci(1);
	xlog=0;ylog=0;

}

void figure::plot(const matrix &y,const char *line) {

	matrix x;
	
	if(y.nrows()==1&y.ncols()>1) 
		x=vector(1,y.ncols(),y.ncols());
	else
		x=vector_t(1,y.nrows(),y.nrows());
	plot(x,y,line);

}

void figure::axis(double sx0,double sx1,double sy0,double sy1,int sjust) {

	x0=sx0;
	x1=sx1;
	y0=sy0;
	y1=sy1;
	just=sjust;
	axis_set=1;

}

void figure::caxis(double sz0,double sz1) {
	
	z0=sz0;
	z1=sz1;
	caxis_set=1;

}

void figure::hold(int state) {

	hold_state=state;
	
}

void figure::draw(int state) {

	cpgslct(id);
	if(state&&!draw_state) cpgebuf();
	if(!state&&draw_state) cpgbbuf();
	draw_state=state;
	
}

void figure::subplot(int nr,int nc) {

	cpgslct(id);
	cpgsubp(nc,nr);
	
}

void figure::next() {

	cpgslct(id);
	cpgpage();
	
}

void figure::clear() {

	cpgslct(id);	
	cpgeras();
	
}

void figure::semilogx(const matrix &x,const matrix &y,const char *line) {

	xlog=1;
	plot(log10(x),y,line);
	
}

void figure::semilogx(const matrix &y,const char *line) {

	matrix x;
	
	if(y.nrows()==1&y.ncols()>1) 
		x=vector(1,y.ncols(),y.ncols());
	else
		x=vector_t(1,y.nrows(),y.nrows());
	semilogx(x,y,line);

}

void figure::semilogy(const matrix &x,const matrix &y,const char *line) {

	ylog=1;
	plot(x,log10(y),line);
	
}

void figure::semilogy(const matrix &y,const char *line) {

	matrix x;
	
	if(y.nrows()==1&y.ncols()>1) 
		x=vector(1,y.ncols(),y.ncols());
	else
		x=vector_t(1,y.nrows(),y.nrows());
	semilogy(x,y,line);

}

void figure::loglog(const matrix &x,const matrix &y,const char *line) {

	xlog=1;
	ylog=1;
	plot(log10(x),log10(y),line);
	
}

void figure::loglog(const matrix &y,const char *line) {

	matrix x;
	
	if(y.nrows()==1&y.ncols()>1) 
		x=vector(1,y.ncols(),y.ncols());
	else
		x=vector_t(1,y.nrows(),y.nrows());
	loglog(x,y,line);

}

void figure::label(const char *xlabel,const char *ylabel,const char *title) {

	cpgslct(id);
	cpglab(xlabel,ylabel,title);
	
}

void figure::pcolor(const matrix &x,const matrix &y,const matrix &z) {

	int i,j,k,c,cmax;
	double zz;
	float xf[4],yf[4];
	
	if(!hold_state) {
		if(!axis_set) {
			x0=min(x);x1=max(x);y0=min(y);y1=max(y);just=0;
			if(x0==x1) {x0=x0-1;x1=x1+1;}
			if(y0==y1) {y0=y0-1;y1=y1+1;}
		}
		if(!caxis_set) {
			z0=min(z);z1=max(z);
			if(z0==z1) {z0=z0-1;z1=z1+1;}
		}
	}
	axis_set=0;
	caxis_set=0;
	cpgslct(id);
	if(draw_state) cpgbbuf();
	cpgsch(1);
	cpgqcol(&i,&cmax);
	if(!hold_state) 
		cpgenv(x0,x1,y0,y1,just,10*xlog+20*ylog);
	for(j=0;j<z.ncols()-1;j++) {
		for(i=0;i<z.nrows()-1;i++) {
			if(x.nrows()==1) {
				xf[0]=x(j);xf[1]=x(j+1);xf[2]=x(j+1);xf[3]=x(j);
			} else {
				xf[0]=x(i,j);xf[1]=x(i,j+1);xf[2]=x(i+1,j+1);xf[3]=x(i+1,j);
			}
			if(y.ncols()==1) {
				yf[0]=y(i);yf[1]=y(i);yf[2]=y(i+1);yf[3]=y(i+1);
			} else {
				yf[0]=y(i,j);yf[1]=y(i,j+1);yf[2]=y(i+1,j+1);yf[3]=y(i+1,j);
			}
			zz=(z(i,j)-z0)/(z1-z0);
			zz=zz>1?1:zz;
			zz=zz<0?0:zz;
			c=round(zz*(cmax-16)+16);
			cpgsci(c);
			cpgpoly(4,xf,yf);
		}
	}
	if(draw_colorbar) {
		cpgsch(1);
		cpgsci(1);
		if(z0>1e36) z0=1e36;
		if(z0<-1e36) z0=-1e36;
		if(z1>1e36) z1=1e36;
		if(z1<-1e36) z1=-1e36;
		cpgwedg("RI",1,3,z0,z1,"");
		draw_colorbar=0;
	}
	if(draw_state) cpgebuf();
	cpgsci(1);
	xlog=0;ylog=0;
}

void figure::pcolor(const matrix &z) {

	matrix x,y;
	
	x=vector(1,z.ncols(),z.ncols());
	y=vector_t(1,z.nrows(),z.nrows());
	pcolor(x,y,z);

}

void figure::colorbar(int set) {
	draw_colorbar=set;
}

void figure::contour(const matrix &x,const matrix &y,const matrix &z,const matrix &cc,const char *line) {

	double xf[4],yf[4];
	int ii[4],jj[4];
	
	int ncontours=cc.nrows()*cc.ncols();
	if(!hold_state) {
		if(!axis_set) {
			x0=min(x);x1=max(x);y0=min(y);y1=max(y);just=0;
			if(x0==x1) {x0=x0-1;x1=x1+1;}
			if(y0==y1) {y0=y0-1;y1=y1+1;}
		}
		if(!caxis_set) {
			z0=min(z);z1=max(z);
			if(z0==z1) {z0=z0-1;z1=z1+1;}
		}
	}
	axis_set=0;
	caxis_set=0;
	cpgslct(id);
	cpgsch(1);
	if(draw_state) cpgbbuf();
	if(!hold_state) 
		cpgenv(x0,x1,y0,y1,just,10*xlog+20*ylog);
		
	int line_style=1,color=-1;
	for(int i=0;i<strlen(line);i++) {
		switch(line[i]) {
			case '-': line_style=1;break;
			case '=': line_style=2;break;
			case ';': line_style=3;break;
			case ':': line_style=4;
		}
	}
	for(int i=0;i<strlen(line);i++) {
		switch(line[i]) {
			case 'w': color=0;break;
			case 'k': color=1;break;
			case 'r': color=2;break;
			case 'g': color=3;break;
			case 'b': color=4;break;
			case 'c': color=5;break;
			case 'm': color=6;break;
			case 'y': color=7;
		}
	}
	cpgsls(line_style);
	int cmax;
	{int i;cpgqcol(&i,&cmax);}
	int c;
	for(int n=0;n<ncontours;n++) {
		if(color==-1) {
			double zz;
			zz=(cc(n)-z0)/(z1-z0);
			zz=zz>1?1:zz;
			zz=zz<0?0:zz;
			c=round(zz*(cmax-16)+16);
			cpgsci(c);
		} else cpgsci(color);
		for(int j=0;j<z.ncols()-1;j++) {
			for(int i=0;i<z.nrows()-1;i++) {
				if(x.nrows()==1) {
					xf[0]=x(j);xf[1]=x(j+1);xf[2]=x(j+1);xf[3]=x(j);
				} else {
					xf[0]=x(i,j);xf[1]=x(i,j+1);xf[2]=x(i+1,j+1);xf[3]=x(i+1,j);
				}
				if(y.ncols()==1) {
					yf[0]=y(i);yf[1]=y(i);yf[2]=y(i+1);yf[3]=y(i+1);
				} else {
					yf[0]=y(i,j);yf[1]=y(i,j+1);yf[2]=y(i+1,j+1);yf[3]=y(i+1,j);
				}
				ii[0]=i;ii[1]=i;ii[2]=i+1;ii[3]=i+1;
				jj[0]=j;jj[1]=j+1;jj[2]=j+1;jj[3]=j;
				int ncruces=0;
				for(int k=0;k<4;k++) {
					double za,zb;
					za=z(ii[k],jj[k])-cc(n);
					zb=z(ii[(k+1)%4],jj[(k+1)%4])-cc(n);
					if((za<=0&&zb>0)||(za>=0&&zb<0)) {
						ncruces++;
						if(ncruces<=2) {
							float xx,yy;
							xx=xf[k]-(xf[(k+1)%4]-xf[k])*za/(zb-za);
							yy=yf[k]-(yf[(k+1)%4]-yf[k])*za/(zb-za);
							if(ncruces==1) 
								cpgmove(xx,yy);
							else
								cpgdraw(xx,yy);
						}
					}
				}
			}
		}
	}
	if(draw_state) cpgebuf();
	xlog=0;ylog=0;
	cpgsls(1);
	cpgsci(1);
}

void figure::contour(const matrix &x,const matrix &y,const matrix &z,int ncontours,const char *line) {

	matrix cc(ncontours);
	
	cc=vector(min(z),max(z),ncontours);
	contour(x,y,z,cc,line);
}

void figure::contour(const matrix &z,int ncontours,const char *line) {
	matrix x,y;
	
	x=vector(1,z.ncols(),z.ncols());
	y=vector_t(1,z.nrows(),z.nrows());
	contour(x,y,z,ncontours,line);
}
void figure::contour(const matrix &x,const matrix &y,const matrix &z,const char *line) {
	contour(x,y,z,10,line);
}
void figure::contour(const matrix &z,const char *line) {
	contour(z,10,line);
}

#else

static int msg=0;

figure::figure(const char *device) {
	if(!msg) {
		fprintf(stderr,"Graphics functionality disabled at compile time. Recompile with USE_PGPLOT=1\n");
		msg=1;
	}
}
figure::~figure() {}
void figure::plot(const matrix &x,const matrix &y,const char *line) {}
void figure::plot(const matrix &y,const char *line) {}
void figure::axis(double sx0,double sx1,double sy0,double sy1,int sjust) {}
void figure::caxis(double sz0,double sz1) {}
void figure::hold(int state) {}
void figure::draw(int state) {}
void figure::subplot(int nr,int nc) {}
void figure::next() {}
void figure::clear() {}
void figure::semilogx(const matrix &x,const matrix &y,const char *line) {}
void figure::semilogx(const matrix &y,const char *line) {}
void figure::semilogy(const matrix &x,const matrix &y,const char *line) {}
void figure::semilogy(const matrix &y,const char *line) {}
void figure::loglog(const matrix &x,const matrix &y,const char *line) {}
void figure::loglog(const matrix &y,const char *line) {}
void figure::label(const char *xlabel,const char *ylabel,const char *title) {}
void figure::pcolor(const matrix &z) {}
void figure::pcolor(const matrix &x,const matrix &y,const matrix &z) {}
void figure::colorbar(int set) {}
void figure::contour(const matrix &x,const matrix &y,const matrix &z,int ncontours,const char *line) {}
void figure::contour(const matrix &z,int ncontours,const char *line) {}
void figure::contour(const matrix &x,const matrix &y,const matrix &z,const matrix &contours,const char *line) {}
void figure::contour(const matrix &z,const matrix &contours,const char *line) {}
void figure::contour(const matrix &x,const matrix &y,const matrix &z,const char *line) {}
void figure::contour(const matrix &z,const char *line) {}

#endif


