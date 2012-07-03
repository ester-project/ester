      function gi(i,x,pp)
      implicit double precision(a-h,o-z)
      h=1.d0-x
      if(i.eq.1)gi=h
      if(i.eq.2)gi=x
      if(i.eq.3)gi=h*h*h/(pp*x+1.d0) 
      if(i.eq.4)gi=x*x*x/(pp*h+1.d0) 
      return
      end
