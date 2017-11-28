      subroutine ordgr(xd,yd,nt,ipt,nxi,nyi,xi,yi,ngp,igp)
      implicit double precision (a-h,o-z)
c it organizes grid points for surface reconstruction by
c sorting them according to their belonging to the triangles.
c
c the input parameters are
c     xd,yd = array of dimension n containing the x and y coordinates
c             of the data points, where n is the number of the data
c             points,
c     nt  = number of triangles,
c     ipt = integer array of dimension 3*nt containing the indices of
c           the vertexes of the triangles,
c     nxi = number of grid points in the x coordinates,
c     nyi = number of grid points in the y coordinates,
c     xi,yi = array of dimension nxi and nyi containing the x and y
c             coordinates of the grid points,respectively.
c
c the output parameters are
c     ngp = integer array of dimension 2*nt where the number of grid
c           points belonging to each triangle is to be stored,
c     igp = integer array of dimension nxi*nyi where the indices of the
c           grid points are to be stored according to their belonging
c           to the triangles considered in ascending order numbers.
c
c declaration statement.
c
      dimension xd(*),yd(*),ipt(*),xi(*),yi(*),ngp(*),igp(*)
c
c statement function.
      side(u1,v1,u2,v2,u3,v3) = (u1-u3)*(v2-v3)-(v1-v3)*(u2-u3)
      nt0 = nt
c
c preliminary processing.
      nxi0 = nxi
      nyi0 = nyi
      nxinyi = nxi0*nyi0
c
c determines grid points inside the data area.
      jngp0 = 0
      jngp1 = 2*nt0 + 1
      jigp0 = 0
      jigp1 = nxinyi + 1
      do 160 it0=1,nt0
        ngp0 = 0
        ngp1 = 0
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
        insd = 0
        do 20 ixi=1,nxi0
          if (xi(ixi).ge.xmn .and. xi(ixi).le.xmx) go to 10
          if (insd.eq.0) go to 20
          iximx = ixi - 1
          go to 30
   10     if (insd.eq.1) go to 20
          insd = 1
          iximn = ixi
   20   continue
        if (insd.eq.0) go to 150
        iximx = nxi0
   30   do 140 iyi=1,nyi0
          yii = yi(iyi)
          if (yii.lt.ymn .or. yii.gt.ymx) go to 140
          do 130 ixi=iximn,iximx
            xii = xi(ixi)
            l = 0
            if (side(x1,y1,x2,y2,xii,yii)) 130, 40, 50
   40       l = 1
   50       if (side(x2,y2,x3,y3,xii,yii)) 130, 60, 70
   60       l = 1
   70       if (side(x3,y3,x1,y1,xii,yii)) 130, 80, 90
   80       l = 1
   90       izi = nxi0*(iyi-1) + ixi
            if (l.eq.1) go to 100
            ngp0 = ngp0 + 1
            jigp0 = jigp0 + 1
            igp(jigp0) = izi
            go to 130
  100       if (jigp1.gt.nxinyi) go to 120
            do 110 jigp1i=jigp1,nxinyi
              if (izi.eq.igp(jigp1i)) go to 130
  110       continue
  120       ngp1 = ngp1 + 1
            jigp1 = jigp1 - 1
            igp(jigp1) = izi
  130     continue
  140   continue
  150   jngp0 = jngp0 + 1
        ngp(jngp0) = ngp0
        jngp1 = jngp1 - 1
        ngp(jngp1) = ngp1
  160 continue
      return
      end
