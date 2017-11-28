      function gid(i,x,pp)
      implicit double precision(a-h,o-z)
      if(i.eq.1)gid=-1.d0
      if(i.eq.2)gid= 1.d0
      if(i.eq.3)gid=-((x-1.d0)**2)*(2.d0*pp*x+3.d0+pp)/
     +              ((pp*x+1.d0)**2)
      if(i.eq.4)gid=-x*x*(-3.d0*pp+2.d0*pp*x-3.d0)/
     +              ((-pp+pp*x-1.d0)**2)
      return
      end
