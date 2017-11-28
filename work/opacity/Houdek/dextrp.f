      subroutine extrp(nd,x,y,z,kc,iwk,a,b,c,d)
      implicit double precision (a-h,o-z)
c it adds some external points to the data point set and estimates
c the z coordinate at these points following the shepard method.
c
c the input parameters are
c     nd  = number of data points,
c     x,y,z = arrays of dimension  nd  containing the  x,y  and  z
c             coordinates of the data points,
c     kc  = flag of computation,
c     a,b,c,d = extreme of the rectangular grid points.
c
c the output parameters are
c     x,y,z = arrays of dimension  ndp  containing the  x,y  and  z
c             coordinates of the new set of data points where
c                 ndp = nd+2*int(nd/25)+4,
c     a,b,c,d = extreme of the rectangular region containing the data
c               point set and the rectangular grid points.
c
c the other parameter is
c     iwk = integer array of dimension 4*nd+7 used as work area.
c
c declaration statement.
      dimension x(*),y(*),z(*),iwk(*),dist(5),ipc0(5)
c statement functions.
      dinf(u1,v1,u2,v2)=dmax1(dabs(u1-u2),dabs(v1-v2))
      dsqf(u1,v1,u2,v2)=((u2-u1)**2+(v2-v1)**2)**2
      ndp=iwk(6)
      ncp=iwk(7)
      if(kc.eq.2)go to 200
c  estimates the smallest rectangle containing the data point set
c  and the rectangular grid points.
      ja=0
      ia=7
      jb=0
      ib=ia+nd
      jc=0
      ic=ib+nd
      jd=0
      id=ic+nd
      do 50 i=1,nd
        if(x(i)-a)16,18,20
16        ja=0
18        ja=ja+1
          iwk(ia+ja)=i
          a=x(i)
          go to 30
20      if(b-x(i))26,28,30
26        jb=0
28        jb=jb+1
          iwk(ib+jb)=i
          b=x(i)
30      if(y(i)-c)36,38,40
36        jc=0
38        jc=jc+1
          iwk(ic+jc)=i
          c=y(i)
          go to 50
40      if(d-y(i))46,48,50
46        jd=0
48        jd=jd+1
          iwk(id+jd)=i
          d=y(i)
50    continue
c  estimates the number of points and where they have to be adjointed.
      n=nd/25
c
c     just for test (setting a=0.)
c
c      a=0.d0
c
      bma=b-a
      dmc=d-c
      hx=bma
      hy=dmc
      nx=1
      ny=1
      if(n.eq.0)go to 75
        do 70 i=1,n
          if(hx.gt.hy)go to 60
            ny=ny+1
            hy=dmc/dfloat(ny)
          go to 70
60          nx=nx+1
            hx=bma/dfloat(nx)
70      continue
      hx=bma/dfloat(nx)
      hy=dmc/dfloat(ny)
c  adds the new external points and checks that they are not
c  coincident with the old ones.
75    ndp=nd+1
      yp=c
      do 100 i=1,ny
        if(ja.eq.0)go to 80
          do 77 j=1,ja
            if( abs(yp-y(iwk(ia+j)) ).le.1.d-10)go to 90
77        continue
80      x(ndp)=a
        y(ndp)=yp
        ndp=ndp+1
90      yp=yp+hy
100   continue
      xp=a
      do 130 i=1,nx
        if(jd.eq.0)go to 110
          do 105 j=1,jd
            if( abs(xp-x(iwk(id+j)) ).le.1.d-10)go to 120
105       continue
110     x(ndp)=xp
        y(ndp)=d
        ndp=ndp+1
120     xp=xp+hx
130   continue
      yp=d
      do 160 i=1,ny
        if(jb.eq.0)go to 140
          do 135 j=1,jb
            if( abs(yp-y(iwk(ib+j)) ).le.1.d-10)go to 150
135       continue
140     x(ndp)=b
        y(ndp)=yp
        ndp=ndp+1
150     yp=yp-hy
160   continue
      xp=b
      do 190 i=1,nx
        if(jc.eq.0)go to 170
          do 165 j=1,jc
            if( abs(xp-x(iwk(ic+j)) ).le.1.d-10)go to 180
165       continue
170     x(ndp)=xp
        y(ndp)=c
        ndp=ndp+1
180     xp=xp-hx
190   continue
      ndp=ndp-1
      iwk(6)=ndp
      ncp=5
      if(nd.le.5)ncp=3
      iwk(7)=ncp
c  estimates the function value at the new external points.
c
200   n1=nd+1
        do 270 ip1=n1,ndp
          x0=x(ip1)
          y0=y(ip1)
          dmx=0.d0
            do 220 ip2=1,nd
              dm=dinf(x0,y0,x(ip2),y(ip2))
              dist(ip2)=dm
              ipc0(ip2)=ip2
                if(dm.le.dmx)go to 210
                  dmx=dm
                  jmx=ip2
210             if(ip2.ge.ncp)go to 230
220         continue
230       ip2p1=ip2+1
            do 250 ip2=ip2p1,nd
              dm=dinf(x0,y0,x(ip2),y(ip2))
                if(dm.ge.dmx)go to 250
                  dist(jmx)=dm
                  ipc0(jmx)=ip2
                  dmx=0.d0
                    do 240 j1=1,ncp
                      if(dist(j1).le.dmx)go to 240
                        dmx=dist(j1)
                        jmx=j1
240                 continue
250         continue
        anum=0.d0
        aden=0.d0
          do 260 j2=1,ncp
            ip2=ipc0(j2)
            r4=dsqf(x0,y0,x(ip2),y(ip2))
              if(r4.eq.0)go to 260
                anum=anum+z(ip2)/r4
                aden=aden+1/r4
260       continue
        z(ip1)=anum/aden
270   continue
      return
      end
