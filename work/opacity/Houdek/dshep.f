      double precision function shep(ncp,x,y,z,x0,y0)
      implicit double precision (a-h,o-z)
c
c extrapolate the z0-value at coordinate (x0,y0) following the shepard method.
c
c the input parameters are
c     ncp  = number of data points,
c     x,y,z = arrays of dimension  ncp  containing the  x,y  and  z
c             coordinates of the data points
c
c the output parameter is
c     shep =  extrapolated value at the coordinate x0,y0
c
c declaration statement.
      double precision x0,y0
      dimension x(ncp),y(ncp),z(ncp)
c
      intrinsic dmax1,dabs
c
c statement functions.
c
c     dinf(u1,v1,u2,v2)=dmax1(dabs(u1-u2),dabs(v1-v2))
      dsqf(u1,v1,u2,v2)=((u2-u1)**2+(v2-v1)**2)**2
c
c  estimates the function value z0 at the new external point (x0,y0).
c
        anum=0.d0
        aden=0.d0
          do 260 i=1,ncp
            r4=dsqf(x0,y0,x(i),y(i))
              if(r4.eq.0.d0)go to 260
                anum=anum+z(i)/r4
                aden=aden+1/r4
260       continue
        shep=anum/aden
      return
      end
