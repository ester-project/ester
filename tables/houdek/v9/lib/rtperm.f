      subroutine rtperm(ivj,m,m1,n,ndim,n1,n2,p,u,p3,ax,bx,cx,dx,rx)
      implicit double precision (a-h,o-z)
      dimension p(ndim,*),u(ndim,*),ax(n),bx(n),cx(n),dx(n),rx(n)
      do 40 j=1,m,m1
         do 20 i=1,n1
            i1=i-1
            if(ivj.ne.0)then
              r2=p3*dx(i)*dx(i)*(u(j,i+1)-u(j,i))
              if(i.eq.1)goto 10
              rx(i1)=r1+r2
              if(i.eq.2) rx(i1)=rx(i1)-dx(1) *p(j,1)
              if(i.eq.n1)rx(i1)=rx(i1)-dx(n1)*p(j,n)
            else
              r2=p3*dx(i)*dx(i)*(u(i+1,j)-u(i,j))
              if(i.eq.1)goto 10 
              rx(i1)=r1+r2 
              if(i.eq.2) rx(i1)=rx(i1)-dx(1) *p(1,j)
              if(i.eq.n1)rx(i1)=rx(i1)-dx(n1)*p(n,j)
            endif
  10        r1=r2
  20     continue
         call trdisb(n2,ax,bx,cx,rx)
         do 30 i=2,n1
            if(ivj.eq.0)p(i,j)=rx(i-1)
            if(ivj.ne.0)p(j,i)=rx(i-1)
  30     continue
  40  continue
      return
      end
