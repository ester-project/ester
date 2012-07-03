      subroutine condux( ro6, t8, zz, iel, cond )
      implicit none

C==============================================
C define parameters                            
C==============================================

C==============================================
C define common blocks
C==============================================
      common /output/ gam,rs,s1,s2,s,fact
C==============================================
C define arguments
C==============================================
      integer iel
      double precision cond
      double precision ro6
      double precision t8
      double precision zz

C==============================================
C define local variables
C==============================================
      double precision a
      double precision fact
      double precision gam
      double precision rs
      double precision s
      double precision s1
      double precision s2
      double precision xx

C**********************************************
C executable statements of routine
C**********************************************
      cond = -1.d0
      a = zz*2.d0
      gam = 0.2275d0*zz**2/t8*(ro6/a)**(1.d0/3.d0)
      xx = 0.45641d0*dlog(gam)-1.31636d0
      rs = 0.01388d0*(a/zz)**(1.d0/3.d0)/ro6**(1.d0/3.d0)
      if (rs .gt. 1.0d1) then
        rs = 1.0d1
      endif
      call sfa( rs,xx,iel,s1,s2 )
      fact = 1.018d0*(zz/a)**(2.d0/3.d0)*ro6**(2.d0/3.d0)
      s = s1-fact/(1.d0+fact)*s2
      if (s .lt. 0.1d0) then
        s = 0.1d0
      endif
      cond = 2.363d+17*(ro6*t8/a)/((1.d0+fact)*s)
c
      end


