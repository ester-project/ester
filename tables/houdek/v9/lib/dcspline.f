      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 45H6:~K.
c------------------------------------------------------------
      SUBROUTINE splint(xa,ya,y2a,n,x,y,ydx)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      DOUBLE PRECISION ydx
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     +  ((a**3-a)*y2a(klo) +
     +   (b**3-b)*y2a(khi)
     +  )*(h**2)/6.d0
c
c     included derivative ydx=(dy)/(dx) hg: 28.12.1993
c
      ax= -1.d0/h
      bx= -ax
      ydx= ax*ya(klo)+bx*ya(khi)+
     +     ((3.d0*a*a-1.d0)*ax*y2a(klo)+
     +      (3.d0*b*b-1.d0)*bx*y2a(khi)
     +     )*(h**2)/6.d0
c
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 45H6:~K.
c------------------------------------------------------------
      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2ax,y2ay)
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      INTEGER m,n,NN
      DOUBLE PRECISION x1a(m),x2a(n),y2ax(m,n),y2ay(m,n),ya(m,n)
      PARAMETER (NN=100)
CU    USES spline
      INTEGER j,k
      DOUBLE PRECISION y2tmp(NN),ytmp(NN)
c
c     calculate 2nd derivatives of the rows -> y2ax
c
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do 12 k=1,n
          y2ax(j,k)=y2tmp(k)
12      continue
13    continue
c
c     calculate 2nd derivatives of the columns -> y2ay
c
      do 23 k=1,n
        do 21 j=1,m
          ytmp(j)=ya(j,k)
21      continue
        call spline(x1a,ytmp,m,1.d30,1.d30,y2tmp)
        do 22 j=1,m
          y2ay(j,k)=y2tmp(j)
22      continue
23    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 45H6:~K.
