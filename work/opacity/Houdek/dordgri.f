      subroutine ordgri(xd,yd,nt,ipt,xi,yi,it0,ntp,ier)
      implicit double precision (a-h,o-z)
c
c it organizes grid points for surface reconstruction by
c sorting them according to their belonging to the triangles.
c
c modified version for only 1 given coordinate pair xi,yi
c
c     History:
c     10.6.1992: creation,derived from ACM,TOMS Alg. 677, ordgr.f
c                search triangle# with the algorithm 
c		 ntp,ntp+1,ntp-1,ntp+2,ntp-2,....
c     25.8.1992: reduce of iwk-      & wk-arrays to
c                          iwk(3*nt) & wk(2*n) respectively
c                where nt = # of triangles (1568) and n = # of datapoints
c                (17*50 = 850)
c      8.9.1992: define new variable ntpr instead of ntp for the case
c                when the searched triangle# it0 is 1, then the ntp-value
c                is also 1 instead of 0 !
c
c     12.3.1993: included ier argument for recognising
c                if desired interpolation point (xi,yi) lies
c                outside triangulation-domain
c
c      7.4.1993: use data statement in s/r masubi for initial value 
c                of ntp (ntp=1 at the very first call -> dismiss 
c                change of 8.9.1992)
c                included if statement for ntp to check if ntp <= nt !
c                ntp > nt happens, when the previous table had more
c                table-points (e.g. transition from Kurucz -> Livermore tab.)
c
c     Last modification:
c      7.4.1993
c
c the input parameters are
c     xd,yd = array of dimension n containing the x and y coordinates
c             of the data points, where n is the number of the data
c             points,
c     nt  = number of triangles,
c     ipt = integer array of dimension 3*nt containing the indices of
c           the vertexes of the triangles,
c     xi,yi = containing the x and y
c             coordinate of the grid point,respectively.
c     ntp = previously used triangle number
c
c the output parameters are
c     it0 = searched triangle number in which the point xi,yi lies
c     ier = if ier.gt.0 then the desired points lies outside the
c           triangulation domain
c
c declaration statement.
c
      dimension xd(*),yd(*),ipt(*)
c
c statement function.
      side(u1,v1,u2,v2,u3,v3) = (u1-u3)*(v2-v3)-(v1-v3)*(u2-u3)
c
      if(ntp.gt.nt)ntp=nt
      idelr=nt-ntp
      idell=ntp-1
      idelm=max0(idelr,idell)
      i=0
c
c determines grid points inside the data area.
c search-order: ntp,ntp+1,ntp-1,ntp+2,ntp-2,....
c
c     ntp
c
        it0=ntp
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        ip2 = ipt(it0t3-1)
        ip3 = ipt(it0t3)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = dmin1(x1,x2,x3)
        xmx = dmax1(x1,x2,x3)
        ymn = dmin1(y1,y2,y3)
        ymx = dmax1(y1,y2,y3)
        if (xi.ge.xmn .and. xi.le.xmx) then
          if (yi.lt.ymn .or. yi.gt.ymx) go to 90
          if (side(x1,y1,x2,y2,xi,yi)) 90, 50, 50
   50     if (side(x2,y2,x3,y3,xi,yi)) 90, 70, 70
   70     if (side(x3,y3,x1,y1,xi,yi)) 90, 80, 80
   80     return
        endif
c
c     from ntp -> nt (ntp+i)
c
   90 i=i+1
      if(i.le.idelr)then
        it0=ntp+i
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        ip2 = ipt(it0t3-1)
        ip3 = ipt(it0t3)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = dmin1(x1,x2,x3)
        xmx = dmax1(x1,x2,x3)
        ymn = dmin1(y1,y2,y3)
        ymx = dmax1(y1,y2,y3)
        if (xi.ge.xmn .and. xi.le.xmx) then
          if (yi.lt.ymn .or. yi.gt.ymx) go to 190
          if (side(x1,y1,x2,y2,xi,yi)) 190, 150, 150
  150     if (side(x2,y2,x3,y3,xi,yi)) 190, 170, 170
  170     if (side(x3,y3,x1,y1,xi,yi)) 190, 180, 180
  180     return
        endif
      endif
c
c     from ntp -> 1 (ntp-i)
c
  190 if(i.le.idell)then
        it0=ntp-i
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        ip2 = ipt(it0t3-1)
        ip3 = ipt(it0t3)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = dmin1(x1,x2,x3)
        xmx = dmax1(x1,x2,x3)
        ymn = dmin1(y1,y2,y3)
        ymx = dmax1(y1,y2,y3)
        if (xi.ge.xmn .and. xi.le.xmx) then
          if (yi.lt.ymn .or. yi.gt.ymx) go to 290
          if (side(x1,y1,x2,y2,xi,yi)) 290, 250, 250
  250     if (side(x2,y2,x3,y3,xi,yi)) 290, 270, 270
  270     if (side(x3,y3,x1,y1,xi,yi)) 290, 280, 280
  280     return
        endif
      endif
  290 if(i.lt.idelm)goto 90
c
c     desired interpolation point (xi,yi) lies outside the
c     triangulation-area
c
      ier=1
      return
c
      end
