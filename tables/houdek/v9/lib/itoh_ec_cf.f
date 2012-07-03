      subroutine itohec( xch, zch, t, rho, ec )
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
      double precision ec
      double precision rho
      double precision t
      double precision xch
      double precision zch

C==============================================
C define local variables
C==============================================
      double precision cond
      double precision conmed
      integer i
      integer iel
      double precision ro6
      double precision t8
      double precision xa
      double precision ya
      double precision za
      double precision zz

C**********************************************
C executable statements of routine
C**********************************************
      xa = xch
      za = zch
      ya = 1.d0-xch-zch
      call itohin( xa,ya,za )
      t8 = 10.d0**(t-8.d0)
      ro6 = 10.d0**(rho-6.d0)
      conmed = 0.d0
      do i = 1, nel
        zz = z(i)
        iel = i
        call condux( ro6,t8,zz,iel,cond )
        conmed = conmed+cond*y(i)
      end do
      ec = dlog10(3.024272d0)-4.d0+3.d0*t-rho-dlog10(conmed)
c
      end


