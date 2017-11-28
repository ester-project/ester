      function sbilin(x,n,y,m,u,ndim,v,w,iflag)
      implicit double precision(a-h,o-z)
      dimension x(n),y(m),u(ndim,m)
      data i,j /1,1/
      iflag=0
      if(n.lt.2.or.m.lt.2) then
        iflag=1
        sbilin=0.0d0
        return
      endif
      call inttwo(x,n,y,m,v,w,i,j,iflag)
      if(iflag.ne.0)then
        sbilin=0.0d0
        return
      endif
      dv=v-x(i)
      dw=w-y(j)
      dx=x(i+1)-x(i)
      dy=y(j+1)-y(j)
      h1=u(i,j)
      h2=u(i+1,j)
      h3=u(i,j+1)
      h4=u(i+1,j+1)
      sbilin=h1+(h2-h1)/dx*dv+(h3-h1)/dy*dw+(h4-h3-h2+h1)/(dx*dy)*dv*dw
      return
      end
