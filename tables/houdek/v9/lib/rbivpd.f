      subroutine rbivpd(u,v,n,ndim,m,mdim,x,y,a,pp,z,zx,zy,iflag)
      implicit double precision (a-h,o-z)
      real a
      dimension x(n),y(m),a(ndim,mdim,4,4)
      data i,j /2*1/
      if(n.lt.2.or.m.lt.2)then
        iflag=1
        return
      endif
      call inttwo(x,n,y,m,u,v,i,j,iflag)
      if(iflag.ne.0)return
      ux=(u-x(i))/(x(i+1)-x(i))
      vy=(v-y(j))/(y(j+1)-y(j))
      z =0.d0
      zx=0.d0
      zy=0.d0
      do 20 k=1,4
         h=gi(k,ux,pp)
         hx=gid(k,ux,pp)
         do 10 l=1,4
            z =z +dble(a(i,j,k,l))*h*gi(l,vy,pp)
            zx=zx+dble(a(i,j,k,l))*hx*gi(l,vy,pp)/(x(i+1)-x(i))
            zy=zy+dble(a(i,j,k,l))*h*gid(l,vy,pp)/(y(j+1)-y(j))
  10     continue
  20  continue
      return
      end
