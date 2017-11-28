      subroutine masubi(drelpr,nt,xd,yd,zd,xi,yi,zi,pdx,pdy,
     +                  nout,iwk,wk,ier)
      implicit double precision (a-h,o-z)
c it carries out a smooth bivariate interpolation when the data
c points projection in the  x-y  plane is irregularly distributed
c in the plane: it estimates the interpolant function at prescribed
c rectangular grid points by means of the 9-parameters discretized
c version of the nielson triangular interpolant. the required derivative
c values are calculated by minimizing suitable tension's functional
c and the triangulation of the convex hull of the x-y data is realized
c by using lawson's local optimization procedure.
c if extrapolation is requested (iex parameter) it adds to the data
c points some additional data, evaluated by means of the shepard's
c method, so that the new convex hull of the data points contains
c the rectangular domain where the surface has to be reconstructed.
c
c     History:
c
c     10.6.1992: creation, derived from ACM,TOMS ALG.677 masub.f
c
c     25.8.1992: reduce of iwk-      & wk-arrays to
c                          iwk(3*nt) & wk(2*n) respectively
c                where nt = # of triangles (1568) and n = # od datapoints
c                (nxi*nyi = 17*50 = 850)
c
c      6.9.1992: use modified interpolation routine intrpi to get
c                partial derivative values pdx,pdy
c
c     2.10.1992: inserted line 'data itpv/0/', because some compiler
c                (e.g. apollo/domain_os) don't initiate the variable
c                itpv to 0 with the fortran-implicit rule, 
c                therefore itpv has any arbitrary value at the very begin !
c
c     12.3.1993: included error-argument ier
c                (and in s/r ordgri)
c                included nt (nr. of triangles) into
c                the argument-list, and removed the
c                common /jpoint/
c                common /ipoint/
c                statements
c
c      7.4.1993: included save-statement and changed
c                initial value for itpv to /1/
c
c     Last modification:
c      7.4.1993
c
c the input parameters are
c     drelpr = double relative precision,
c     nt  = nr. of triangels
c     xd  = array of dimension n containing the  x  coordinates of
c           the nd data points,
c     yd  = array of dimension n containing the  y  coordinates of
c           the nd data points,
c     zd  = array of dimension n containing the  z  coordinates of
c           the nd data points,
c     xi  = the x coordinate of the rectangular grid point,  where the
c           surface has to be reconstructed,
c     yi  = the y oordinates of the rectangular grid point,  where the
c           surface has to be reconstructed,
c     nout= logical unit number for the standard output unit of the
c           system.
c
c the output parameters are
c     xd,yd,zd = arrays of dimension n, containing the x, y and z
c                coordinates of the nd data points and,if extrapolation
c                is requested (iex=1), the x, y and z coordinates of
c                n-nd additional points evaluated by means of shepard's
c                method. if extrapolation is not requested (iex=0) then
c                n=nd and xd, yd, zd are not changed,
c     zi  = interpolated  z  values at the rectangular grid point is
c           to be stored.
c           if extrapolation is not requested (iex=0), the z values at
c           the rectangular grid points outside the convex hull of the
c           data points, are set to 1.0e6.
c     pdx = interpolated partial derivative value in x-direction
c     pdy = interpolated partial derivative value in y-direction
c     ier = if ier.gt.0 then an error occured (no valid result)
c
c the other parameters are
c     iwk = integer array of dimension 3*nt
c           used internally as work area (indices of triangle-vertexes)
c     wk  = array of dimension 2*n used internally as work area.
c           (n=# of datapoints [850]) (partial derivative values)
c
c declaration statement.
      dimension xd(*),yd(*),zd(*),iwk(*),wk(*)
      integer iti
c
      data itpv /1/
c
      save
c
      toll=1.0d0-drelpr
c
c sort rectangular grid point according to its belonging to the
c triangle.
      ier=0
      call ordgri(xd,yd,nt,iwk(1),xi,yi,iti,itpv,ier)
c
c     zi=99.d0
      if(ier.eq.0)
     +  call intrpi(drelpr,toll,xd,yd,zd,iwk(1),wk(1),iti,itpv,
     +                          xi,yi,zi,pdx,pdy)
      return
c
      end
