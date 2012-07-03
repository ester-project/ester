      subroutine aditohec( xch, zch, t, rho, adxch, adzch, adt, adrho, 
     $adec )
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
      double precision adec
      double precision adrho
      double precision adt
      double precision adxch
      double precision adzch
      double precision rho
      double precision t
      double precision xch
      double precision zch

C==============================================
C define local variables
C==============================================
      double precision adcond
      double precision adconmed
      double precision adro6
      double precision adt8
      double precision adxa
      double precision adya
      double precision adza
      double precision cond
      double precision conmed
      integer i
      integer iel
      integer ip1
      double precision ro6
      double precision t8
      double precision xa
      double precision y_h(nel)
      double precision ya
      double precision za
      double precision zz

C----------------------------------------------
C SAVE ARGUMENTS
C----------------------------------------------
      do ip1 = 1, nel
        y_h(ip1) = y(ip1)
      end do

C----------------------------------------------
C RESET GLOBAL ADJOINT VARIABLES
C----------------------------------------------
      call adzero

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adcond = 0
      adconmed = 0
      adro6 = 0
      adt8 = 0
      adxa = 0
      adya = 0
      adza = 0

C----------------------------------------------
C ADJOINT STATEMENTS
C----------------------------------------------
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
      adconmed = adconmed-adec*(1.d0/(conmed*dlog(10.d0)))
      adrho = adrho-adec
      adt = adt+3*adec
      adec = 0
      do i = nel, 1, -1
        zz = z(i)
        iel = i
        call condux( ro6,t8,zz,iel,cond )
        adcond = adcond+adconmed*y(i)
        ady(i) = ady(i)+adconmed*cond
        call adcondux( ro6,t8,zz,iel,adro6,adt8,adcond )
      end do
      adrho = adrho+adro6*10.d0**(rho-6.d0)*dlog(10.d0)
      adro6 = 0
      adt = adt+adt8*10.d0**(t-8.d0)*dlog(10.d0)
      adt8 = 0
      do ip1 = 1, nel
        y(ip1) = y_h(ip1)
      end do
      call aditohin( xa,ya,za,adxa,adya,adza )
      adxch = adxch-adya
      adzch = adzch-adya
      adya = 0
      adzch = adzch+adza
      adza = 0
      adxch = adxch+adxa
      adxa = 0

      end
c
      subroutine adzero
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

C==============================================
C define local variables
C==============================================
      integer ip1

      do ip1 = 1, nel
        ady(ip1) = 0
      end do
c
      end
