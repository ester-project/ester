      function rbival(u,v,n,ndim,m,mdim,x,y,a,pp,iflag)
      implicit double precision (a-h,o-z)
      dimension x(n),y(m),a(ndim,mdim,4,4)
      data i,j /2*1/
      if(n.lt.2.or.m.lt.2)then
        iflag=1
        rbival=0.0d0
        return
      endif
      call inttwo(x,n,y,m,u,v,i,j,iflag)
      if(iflag.ne.0)then
         rbival=0.0d0
         return
      endif
      ux=(u-x(i))/(x(i+1)-x(i))
      vy=(v-y(j))/(y(j+1)-y(j))
      rbival=0.d0
      do 20 k=1,4
         h=gi(k,ux,pp)
         do 10 l=1,4
            rbival=rbival+a(i,j,k,l)*h*gi(l,vy,pp)
  10     continue
  20  continue
      return
      end
