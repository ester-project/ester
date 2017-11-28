      subroutine masub(drelpr,ic,iex,nd,xd,yd,zd,tp,nxi,nyi,xi,yi,zi,
     *                 ndim,nout,iwk,wk)
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
c
c the input parameters are
c     drelpr = double relative precision,
c     ic  = flag of computation (must be 1,2,3 or 4),
c         = 1 for the first call and for new iex,nd,xd-yd,
c         = 2 for new zd,tp,
c         = 3 for new tp,
c         = 4 for new nxi,nyi,xi-yi,
c     iex = flag of computation (must be 0 or 1),
c         = 1 (one) allows the extrapolation, 0 (zero) does not;
c           this parameter must be the same in the following call,
c     nd  = number of data points (must be 4 or greater);
c           given iex and nd we define n as
c                       n = nd + iex*(2*int(nd/25)+4),
c     xd  = array of dimension n containing the  x  coordinates of
c           the nd data points,
c     yd  = array of dimension n containing the  y  coordinates of
c           the nd data points,
c     zd  = array of dimension n containing the  z  coordinates of
c           the nd data points,
c     tp  = tension parameter (must be greater than or equal to zero).
c     nxi = number of rectangular grid points in the x coordinate (must
c           be 1 or greater),
c     nyi = number of rectangular grid points in the y coordinate (must
c           be 1 or greater),
c     xi  = array of dimension nxi containing,in ascending order, the x
c           coordinates of the rectangular grid points,  where the
c           surface has to be reconstructed,
c     yi  = array of dimension nyi containing,in ascending order, the y
c           coordinates of the rectangular grid points,  where the
c           surface has to be reconstructed,
c     ndim= declared row dimension of the array containing zi,
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
c     zi  = doubly-dimensioned array of dimension nxi*nyi, where the
c           interpolated  z  values at the rectangular grid points are
c           to be stored.
c           if extrapolation is not requested (iex=0), the z values at
c           the rectangular grid points outside the convex hull of the
c           data points, are set to 1.0e6.
c
c the other parameters are
c     iwk = integer array of dimension  23*n-31+max0(8*n+25,nxi*nyi)
c           used internally as work area,
c     wk  = array of dimension  19*n-20  used internally as work area.
c
c the first call to this subroutine and the call with a new nd value
c and/or new contents of the xd and yd arrays must be made with ic=1.
c the call with ic=2 must be preceded by another call with the same
c iex,nd,nxi and nyi values and with the same contents of the xd,yd,
c xi and yi arrays.
c the call with ic=3 must be preceded by another call with the same
c iex,nd,nxi and nyi values and with the same contents of the xd,yd,
c zd,xi and yi arrays.
c the call with ic=4 must be preceded by another call with the same
c iex and nd values and with the same contents of the xd,yd and zd
c arrays.
c iwk and wk arrays must not be disturbed between the call with ic not
c equal to  1  and the preceding one.
c subroutine masub calls the extrp,ctang,adjac,pdste,pdmin,ordgr and
c intrp subroutines.
c
c declaration statement.
      dimension xd(*),yd(*),zd(*),xi(*),yi(*),zi(ndim,*),iwk(*),wk(*)
      toll=1.0d0-drelpr
c  error check.
      if(ic.lt.1.or.ic.gt.4)go to 1000
      if(iex.lt.0.or.iex.gt.1)go to 1000
      if(nd.lt.4)go to 1000
      if(nxi.lt.1.or.nyi.lt.1)go to 1000
      if(ic.gt.1)go to 10
      iwk(1)=nd
      iwk(2)=iex
      iwk(3)=nxi
      iwk(4)=nyi
      iwk(6)=nd
        if(iex.eq.0)go to 50
        wk(1)=xi(1)
        wk(2)=xi(nxi)
        wk(3)=yi(1)
        wk(4)=yi(nyi)
        go to 40
10    if(iwk(1).ne.nd) go to 1000
      if(iwk(2).ne.iex)go to 1000
      if(ic.gt.3) go to 20
      if(iwk(3).ne.nxi)go to 1000
      if(iwk(4).ne.nyi)go to 1000
        if(ic.eq.2.and.iex.eq.1)go to 40
        go to 50
20    iwk(3)=nxi
      iwk(4)=nyi
        if(iex.eq.0)go to 50
        if(xi(1).lt.wk(1).or.xi(nxi).gt.wk(2))go to 1000
        if(yi(1).lt.wk(3).or.yi(nyi).gt.wk(4))go to 1000
        go to 50
c adds some external points for the extrapolation.
40    call extrp(nd,xd,yd,zd,ic,iwk,wk(1),wk(2),wk(3),wk(4))
50    ndp=iwk(6)
c allocation storage areas in the iwk and wk arrays.
      jipt=8
      jipl=6*ndp-7
      jiwp=12*ndp-7
      jind=13*ndp-7
      jngp=19*ndp-22
      jigp=23*ndp-32
      ipd=5
      ial=2*ndp+5
      ibe=5*ndp-1
      iga=8*ndp-7
      iei=11*ndp-13
      ials=14*ndp-19
      ibes=ials+ndp
      igas=ibes+ndp
      izx=igas+ndp
      izy=izx+ndp
      if(ic.gt.1)go to 60
c triangulates the x-y plane.
      call ctang(ndp,xd,yd,nout,nt,iwk(jipt),iwk(jipl),iwk(jind),
     *           iwk(jiwp),wk(ipd))
      iwk(5)=nt
c constructs the adjacencies monodimensional array.
      call adjac(nt,iwk(jipt),ndp,iwk(jipl),iwk(jiwp-1))
60    nt=iwk(5)
      if(ic.gt.3)go to 70
      if(ic.eq.3)go to 65
c estimates partial derivatives at all data points.
      call pdste(ndp,xd,yd,zd,nt,iwk(jipt),wk(ipd),iwk(jind))
65    call pdmin(drelpr,ndp,xd,yd,zd,nout,iwk(jipl),iwk(jiwp-1),
     *     iwk(jind),tp,wk(ipd),wk(ial),wk(ibe),wk(iga),wk(iei),
     *     wk(ials),wk(ibes),wk(igas),wk(izx),wk(izy))
      if(ic.gt.1)go to 80
c sorts rectangular grid points according to their belonging to the
c triangles.
70    call ordgr(xd,yd,nt,iwk(jipt),nxi,nyi,xi,yi,iwk(jngp),
     *           iwk(jigp))
80    do 85 i=1,nxi
      do 85 j=1,nyi
85    zi(i,j)=1.d6
      itpv=0
      imax=0
      i1min=nxi*nyi+1
      do 120 kngp=1,nt
        iti=kngp
        jwngp=jngp-1+kngp
        ngp0=iwk(jwngp)
        if(ngp0.eq.0)go to 100
        imin=imax+1
        imax=imax+ngp0
        do 90 kigp=imin,imax
          jwigp=jigp+kigp-1
          izi=iwk(jwigp)
          iyi=(izi-1)/nxi+1
          ixi=izi-nxi*(iyi-1)
          call intrp(drelpr,toll,xd,yd,zd,iwk(jipt),wk(ipd),iti,itpv,
     *               xi(ixi),yi(iyi),zi(ixi,iyi))
90      continue
100     jwngp=jngp+2*nt-kngp
        ngp1=iwk(jwngp)
        if(ngp1.eq.0)go to 120
        i1max=i1min-1
        i1min=i1min-ngp1
        do 110 kigp=i1min,i1max
          jwigp=jigp+kigp-1
          izi=iwk(jwigp)
          iyi=(izi-1)/nxi+1
          ixi=izi-nxi*(iyi-1)
          call intrp(drelpr,toll,xd,yd,zd,iwk(jipt),wk(ipd),iti,itpv,
     *               xi(ixi),yi(iyi),zi(ixi,iyi))
110     continue
120   continue
      return
c error exit.
1000  write(nout,130)iex,ic,nd
130   format(1x,34himproper input parameter value(s)./
     *       4hiex=,i4,4h ic=,i4,4h nd=,i4)
      return
      end
