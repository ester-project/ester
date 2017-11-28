      subroutine itohin( xa, ya, za )
      implicit none

C==============================================
C define parameters                            
C==============================================
      integer nel
      parameter ( nel = 10 )

C==============================================
C define common blocks
C==============================================
      common /itohcz/ a, atw, y, x, z
      double precision a(nel)
      double precision atw(nel)
      double precision x(nel)
      double precision y(nel)
      double precision z(nel)

C==============================================
C define arguments
C==============================================
      double precision xa
      double precision ya
      double precision za

C==============================================
C define local variables
C==============================================
      double precision a1
      double precision aa(nel)
      double precision at
      double precision atw_m
      integer i
      double precision x_r
      double precision xt
      double precision y_r
      double precision yt
      double precision z_r

C==============================================
C define data
C==============================================
      data z/1.,2.,6.,8.,10.,12.,14.,16.,20.,26./
      data atw/1.008,4.003,12.011,15.9994,20.179,24.305,28.0855,32.06,
     $40.,55.847/
      data a/12.00,11.00,8.55,8.87,8.07,7.58,7.55,7.21,6.36,7.50/

C**********************************************
C executable statements of routine
C**********************************************
      a1 = a(3)
      do i = 3, nel
        aa(i) = 10.d0**(a(i)-a1)
      end do
      at = 0.d0
      do i = 3, nel
        at = at+aa(i)
      end do
      yt = 0.d0
      do i = 3, nel
        y(i) = aa(i)/at
        yt = yt+y(i)
      end do
      atw_m = 0.d0
      do i = 3, nel
        atw_m = atw_m+atw(i)*y(i)
      end do
      if (za .gt. 0.d0) then
        z_r = 1.d0/(1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/(za*atw(2)))
        x_r = z_r*xa*atw_m/(za*atw(1))
        y_r = z_r*ya*atw_m/(za*atw(2))
      else if (xa .gt. 0.d0) then
        z_r = 0.d0
        x_r = 1.d0/(1.d0+ya*atw(1)/(xa*atw(2)))
        y_r = x_r*ya*atw(1)/(xa*atw(2))
      else
        z_r = 0.d0
        x_r = 0.d0
        y_r = 1.d0
      endif
      yt = z_r+y_r+x_r
c     write(*,1003) x_r,y_r,z_r
c1003 format(' x,y,z rel. number abundances > ',3f10.5)
      y(1) = x_r
      y(2) = y_r
      do i = 3, nel
        y(i) = y(i)*z_r
      end do
      atw_m = 0.d0
      do i = 1, nel
        atw_m = atw_m+atw(i)*y(i)
      end do
      xt = 0.d0
      do i = 1, nel
        x(i) = y(i)*atw(i)*yt/atw_m
        xt = xt+x(i)
      end do
c     write(*,1001) atw_m,yt,xt
c1001 format(' mean atomic weight ',f8.4,'; sum of number/mass abnds. ',
c    $2f8.3)
c     do i = 1, nel
c       write(*,1002) z(i),atw(i),y(i),x(i)
c     end do
c1002 format(f6.2,f7.2,1p,2e14.5)
c
      end


