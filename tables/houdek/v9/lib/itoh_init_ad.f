      subroutine aditohin( xa, ya, za, adxa, adya, adza )
      implicit none

C==============================================
C define parameters                            
C==============================================
      integer nel
      parameter ( nel = 10 )

C==============================================
C define common blocks
C==============================================
      common /aditohcz/ ady
      double precision ady(nel)

      common /itohcz/ a, atw, y, x, z
      double precision a(nel)
      double precision atw(nel)
      double precision x(nel)
      double precision y(nel)
      double precision z(nel)

C==============================================
C define arguments
C==============================================
      double precision adxa
      double precision adya
      double precision adza
      double precision xa
      double precision ya
      double precision za

C==============================================
C define local variables
C==============================================
      double precision adatw_m
      double precision adx_r
      double precision ady_r
      double precision adz_r
      double precision atw_m
      integer i
      integer i1
      double precision x_r
      double precision y_r
      double precision z_r

C==============================================
C define data
C==============================================
c     data z/1.,2.,6.,8.,10.,12.,14.,16.,20.,26./
c     data atw/1.008,4.003,12.011,15.9994,20.179,24.305,28.0855,32.06,
c    $40.,55.847/
c     data a/12.00,11.00,8.55,8.87,8.07,7.58,7.55,7.21,6.36,7.50/

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adatw_m = 0
      adx_r = 0
      ady_r = 0
      adz_r = 0

C----------------------------------------------
C ADJOINT STATEMENTS
C----------------------------------------------
      atw_m = 0.d0
      do i = 3, nel
        atw_m = atw_m+atw(i)*y(i)
      end do
      if (za .gt. 0.d0) then
        z_r = 1.d0/(1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/(za*atw(2)))
        y_r = z_r*ya*atw_m/(za*atw(2))
      else if (xa .gt. 0.d0) then
        z_r = 0.d0
        x_r = 1.d0/(1.d0+ya*atw(1)/(xa*atw(2)))
        y_r = x_r*ya*atw(1)/(xa*atw(2))
      else
        z_r = 0.d0
        y_r = 1.d0
      endif
      do i = nel, 3, -1
        y(2) = y_r
        do i1 = 3, i-1
          y(i1) = y(i1)*z_r
        end do
        adz_r = adz_r+ady(i)*y(i)
        ady(i) = ady(i)*z_r
      end do
      ady_r = ady_r+ady(2)
      ady(2) = 0
      adx_r = adx_r+ady(1)
      ady(1) = 0
      if (za .gt. 0.d0) then
        z_r = 1.d0/(1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/(za*atw(2)))
        adatw_m = adatw_m+ady_r*(z_r*ya/(za*atw(2)))
        adya = adya+ady_r*(z_r*atw_m/(za*atw(2)))
        adz_r = adz_r+ady_r*(ya*atw_m/(za*atw(2)))
        adza = adza-ady_r*(z_r*ya*atw_m*atw(2)/(za*atw(2)*za*atw(2)))
        ady_r = 0
        adatw_m = adatw_m+adx_r*(z_r*xa/(za*atw(1)))
        adxa = adxa+adx_r*(z_r*atw_m/(za*atw(1)))
        adz_r = adz_r+adx_r*(xa*atw_m/(za*atw(1)))
        adza = adza-adx_r*(z_r*xa*atw_m*atw(1)/(za*atw(1)*za*atw(1)))
        adx_r = 0
        adatw_m = adatw_m-adz_r*(1.d0*(xa/(za*atw(1))+ya/(za*atw(2)))/
     $((1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/(za*atw(2)))*(1.d0+xa*atw_m/
     $(za*atw(1))+ya*atw_m/(za*atw(2)))))
        adxa = adxa-adz_r*(1.d0*(atw_m/(za*atw(1)))/((1.d0+xa*atw_m/(za*
     $atw(1))+ya*atw_m/(za*atw(2)))*(1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/
     $(za*atw(2)))))
        adya = adya-adz_r*(1.d0*(atw_m/(za*atw(2)))/((1.d0+xa*atw_m/(za*
     $atw(1))+ya*atw_m/(za*atw(2)))*(1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/
     $(za*atw(2)))))
        adza = adza-adz_r*(1.d0*((-(xa*atw_m*atw(1)/(za*atw(1)*za*atw(1)
     $)))-ya*atw_m*atw(2)/(za*atw(2)*za*atw(2)))/((1.d0+xa*atw_m/(za*
     $atw(1))+ya*atw_m/(za*atw(2)))*(1.d0+xa*atw_m/(za*atw(1))+ya*atw_m/
     $(za*atw(2)))))
        adz_r = 0
      else if (xa .gt. 0.d0) then
        x_r = 1.d0/(1.d0+ya*atw(1)/(xa*atw(2)))
        adx_r = adx_r+ady_r*(ya*atw(1)/(xa*atw(2)))
        adxa = adxa-ady_r*(x_r*ya*atw(1)*atw(2)/(xa*atw(2)*xa*atw(2)))
        adya = adya+ady_r*(x_r*atw(1)/(xa*atw(2)))
        ady_r = 0
        adxa = adxa+adx_r*(1.d0*(ya*atw(1)*atw(2)/(xa*atw(2)*xa*atw(2)))
     $/((1.d0+ya*atw(1)/(xa*atw(2)))*(1.d0+ya*atw(1)/(xa*atw(2)))))
        adya = adya-adx_r*(1.d0*(atw(1)/(xa*atw(2)))/((1.d0+ya*atw(1)/
     $(xa*atw(2)))*(1.d0+ya*atw(1)/(xa*atw(2)))))
        adx_r = 0
      endif

 1003 format(' x,y,z rel. number abundances > ',3f10.5)
 1001 format(' mean atomic weight ',f8.4,'; sum of number/mass abnds. ',
     $2f8.3)
 1002 format(f6.2,f7.2,1p,2e14.5)
      end


