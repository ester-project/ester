      subroutine rat2d(n,m,ndim,mdim,x,y,u,ir,pp,eps,p,q,r,a,iflag,
     +                 dx,ax,bx,cx,rx,dy,ay,by,cy,ry)
      implicit double precision(a-h,o-z)
      real a
      dimension x(n),y(m),u(ndim,m),p(ndim,m),q(ndim,m),
     +          r(ndim,m),a(ndim,mdim,4,4),
     +          dx(n),ax(n),bx(n),cx(n),rx(n),
     +          dy(m),ay(m),by(m),cy(m),ry(m),
     +          b(4,4),c(4,4),d(4,4),e(4,4)
      iflag=0
      if(n.lt.2.or.m.lt.2)then
        iflag=1
        return
      end if
      if(pp.lt.(eps-1.d0))then
        iflag=3
        return
      endif
      n1=n-1
      n2=n-2
      m1=m-1
      m2=m-2
      zero=0.d0
      do 10 i=1,n1
         dx(i)=1.d0/(x(i+1)-x(i))
  10  continue
      do 20 j=1,m1
         dy(j)=1.d0/(y(j+1)-y(j))
  20  continue
      p1=pp+1.d0
      p2=pp+2.d0
      p3=pp+3.d0
      ga=p3*p1
      dk=1.d0/p1
      b(1,1)=p2*dk
      b(2,1)=-dk
      b(4,3)=-dk
      b(4,1)= dk
      b(2,3)= b(1,1)
      b(1,3)=-dk
      b(3,1)=-dk
      b(3,3)= dk
      do 40 j=1,3,2
         do 30 i=1,4
            e(i,j)=b(i,j)
  30     continue
  40  continue
      if(ir.eq.2) goto 70
      do 50 j=1,m
         p(1,j)=(u(2,j)-u(1,j))*dx(1)
         p(n,j)=(u(n,j)-u(n1,j))*dx(n1)
  50  continue
      do 60 i=1,n
         q(i,1)=(u(i,2)-u(i,1))*dy(1)
         q(i,m)=(u(i,m)-u(i,m1))*dy(m1)
  60  continue
      r(1,1)=((p(1,2)-p(1,1))*dy(1)+(q(2,1)-q(1,1))*dx(1))/2.d0
      r(1,m)=((p(1,m)-p(1,m1))*dy(m1)+(q(2,m)-q(1,m))*dx(1))/2.d0
      r(n,1)=((p(n,2)-p(n,1))*dy(1)+(q(n,1)-q(n1,1))*dx(n1))/2.d0
      r(n,m)=((p(n,m)-p(n,m1))*dy(m1)+(q(n,m)-q(n1,m))*dx(n1))/2.d0
  70  if(n.eq.2.and.m.eq.2) goto 130
      if(n.eq.2) goto 90
      do 80 i=1,n2
         if(i.lt.n2) ax(i)=dx(i+1)
         bx(i)=p2*(dx(i+1)+dx(i))
  80  continue
      call trdisa(n2,ax,bx,cx,eps,iflag)
      if(iflag.ne.0) return
      ivj=0
      call rtperm(ivj,m,1,n,ndim,n1,n2,p,u,p3,ax,bx,cx,dx,rx)
  90  if(m.eq.2) goto 110
      do 100 j=1,m2
         if(j.lt.m2) ay(j)=dy(j+1)
         by(j)=p2*(dy(j+1)+dy(j))
 100  continue
      call trdisa(m2,ay,by,cy,eps,iflag)
      if(iflag.ne.0) return
      ivj=1
      call rtperm(ivj,n,1,m,ndim,m1,m2,q,u,p3,ay,by,cy,dy,ry)
 110  if(n.eq.2) goto 120
      ivj=0
      call rtperm(ivj,m,m1,n,ndim,n1,n2,r,q,p3,ax,bx,cx,dx,rx)
 120  if(m.eq.2) goto 130
      ivj=1
      call rtperm(ivj,n,1,m,ndim,m1,m2,r,p,p3,ay,by,cy,dy,ry)
 130  do 210 i=1,n1
         i1=i+1
         call ratmat(1.d0/dx(i),pp,b)
         do 200 j=1,m1
            j1=j+1
            c(1,1)=u(i,j)
            c(1,2)=q(i,j)
            c(2,1)=p(i,j)
            c(2,2)=r(i,j)
            c(1,3)=u(i,j1)
            c(1,4)=q(i,j1)
            c(2,3)=p(i,j1)
            c(2,4)=r(i,j1)
            c(3,1)=u(i1,j)
            c(3,2)=q(i1,j)
            c(4,1)=p(i1,j)
            c(4,2)=r(i1,j)
            c(3,3)=u(i1,j1)
            c(3,4)=q(i1,j1)
            c(4,3)=p(i1,j1)
            c(4,4)=r(i1,j1)
            do 160 k1=1,4
               do 150 k2=1,4
                  sum=zero
                  do 140 k=1,4
                     sum=sum+b(k1,k)*c(k,k2)
 140              continue
                  d(k1,k2)=sum
 150           continue
 160        continue
            call ratmat(1.d0/dy(j),pp,e)
            do 190 k1=1,4
                 do 180 k2=1,4
                    sum=zero
                    do 170 k=1,4
                       sum=sum+d(k1,k)*e(k2,k)
 170                continue
                    a(i,j,k1,k2)=sngl(sum)
 180           continue
 190        continue
 200     continue
 210  continue
      return
      end
