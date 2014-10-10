#ifndef _GRAPHICS_H
#define _GRAPHICS_H

#include "config.h"
#include "matrix.h"

class figure {
	int id;
	char dev[32];
	int hold_state;
	int axis_set,caxis_set;
	double x0,x1,y0,y1,z0,z1;
	int xlog,ylog;
	int just;
	int draw_state;
	int draw_colorbar;
public:
	explicit figure(const char *device="/NULL");
	~figure();
	void plot(const matrix &x,const matrix &y,const char *line="");
	void plot(const matrix &y,const char *line="");
	void axis(double sx0,double sx1,double sy0,double sy1,int sjust=0);
	void caxis(double sz0,double sz1);
	void hold(int state);
	void draw(int state);
	void subplot(int nr,int nc);
	void next();
	void clear();
	void semilogx(const matrix &x,const matrix &y,const char *line="");
	void semilogx(const matrix &y,const char *line="");
	void semilogy(const matrix &x,const matrix &y,const char *line="");
	void semilogy(const matrix &y,const char *line="");
	void loglog(const matrix &x,const matrix &y,const char *line="");
	void loglog(const matrix &y,const char *line="");
	void label(const char *xlabel,const char *ylabel,const char *title);
	void pcolor(const matrix &z);
	void pcolor(const matrix &x,const matrix &y,const matrix &z);
	void colorbar(int set=1);
	void contour(const matrix &x,const matrix &y,const matrix &z,int ncontours,const char *line="");
	void contour(const matrix &z,int ncontours,const char *line="");
	void contour(const matrix &x,const matrix &y,const matrix &z,const matrix &contours,const char *line="");
	void contour(const matrix &z,const matrix &contours,const char *line="");
	void contour(const matrix &x,const matrix &y,const matrix &z,const char *line="");
	void contour(const matrix &z,const char *line="");
};

#endif
