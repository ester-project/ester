      function pderix(x,y,x1,x2,x3,y1,y2,y3,aij,aik,si,sik,sij)
      implicit double precision (a-h,o-z)
c
c     code was generated with MapleVR2
c
c     this function calculates the analytical derived formula
c     for the partial derivative at the interpolated point (x,y)
c     with respect to x.
c
c     17.8.1993: hg, creation
c
c     input parameters:
c
c     x,y ......... the coordinates of the point 
c     x1,x2,x3 .... the x-coordinate of the triangle in which 
c                   the point (x,y) lies
c     y1,y2,y3 .... the y-coordinate of the triangle in which 
c                   the point (x,y) lies
c     aij,aik ..... coefficients alpha_ij and alpha_ik
c                   (see G.M. Nielson, Math. of Comp. Vol.40,Nr.161
c                        p.253-271, January 1983)
c     si .......... z-value at triangle-vertex Vi 
c                   (see G.M. Nielsen p.258)
c     sij,sik ..... derivatives along triangle sides 
c                   (see G.M. Nielson p.258) 
c
c     output:
c
c     pderix ...... partial derivative with respect to x
c
      t1 = x3*y
      t2 = x*y3
      t3 = x2*y
      t4 = y2*x
      t5 = x2*y3
      t6 = y2*x3
      t7 = t1-t2-t3+t4+t5-t6
      t8 = x2*y1
      t9 = x1*y3
      t10 = y1*x3
      t11 = y2*x1
      t12 = -t8+t5-t9+t10-t6+t11
      t13 = t12**2
      t14 = 1/t13
      t15 = t7*t14
      t16 = 1/t12
      t19 = -y3+y2
      t22 = t7**2
      t24 = 1/t13/t12
      t25 = t22*t24
      t27 = t13**2
      t28 = 1/t27
      t30 = x1*y
      t31 = y1*x
      t32 = -t30+t9-t10+t1-t2+t31
      t33 = t32*t28*t7
      t34 = -t8+t11+t3-t30+t31-t4
      t37 = t32*t14
      t39 = -t15*t32+t15*t34-t37*t34
      t40 = 1/t39
      t41 = t34*t40
      t42 = t34*t16
      t43 = aij*t42
      t44 = t32*t16
      t45 = aik*t44
      t46 = t43-t45
      t50 = t22*t28
      t51 = -y3+y1
      t52 = t50*t51
      t55 = t50*t32
      t56 = y1-y2
      t57 = t56*t40
      t60 = t39**2
      t62 = 1/t60*t34
      t63 = t19*t14
      t71 = -t63*t32-t15*t51+t63*t34+t15*t56-t34*t14*t51-t37*t56
      t75 = t56*t16
      t76 = aij*t75
      t77 = t51*t16
      t78 = aik*t77
      t84 = t7*t24
      t88 = 3.D0*t43-t44-t42
      t107 = -3.D0*t45+t42+t44
      t123 = (2*t19*(3.D0-2.D0*t7*t16)*t15-2.D0*t25*t19-12.D0*t19*t46*t4
     #1*t33-6.D0*t46*t41*t52-6.D0*t46*t57*t55+6.D0*t71*t46*t62*t55-6.D0*
     #(t76-t78)*t41*t55)*si+(2*t19*t34*t84+t25*t56-2*t19*t88*t41*t33-t88
     #*t41*t52-t88*t57*t55+t71*t88*t62*t55-(3.D0*t76-t77-t75)*t41*t55)*s
     #ik+(-2*t19*t32*t84-t25*t51-2*t19*t107*t41*t33-t107*t41*t52-t107*t5
     #7*t55+t71*t107*t62*t55-(-3.D0*t78+t75+t77)*t41*t55)*sij
      pderix=t123
c
      return
      end
