      subroutine adsfa( rs, xx, iel, adrs, adxx, ads1, ads2 )
      implicit none

C==============================================
C define parameters                            
C==============================================

C==============================================
C define common blocks
C==============================================
      common /itohcc/ a, b, c, d, e, f
      double precision a(4,10)
      double precision b(3,10)
      double precision c(3,10)
      double precision d(4,10)
      double precision e(3,10)
      double precision f(3,10)

C==============================================
C define arguments
C==============================================
      integer iel
      double precision adrs
      double precision ads1
      double precision ads2
      double precision adxx
      double precision rs
      double precision xx

C==============================================
C define local variables
C==============================================
      double precision adv1
      double precision adv2
      double precision adv3
      double precision adv4
      double precision adv5
      double precision adv6
      integer k
      double precision v1
      double precision v2
      double precision v3
      double precision v4
      double precision v5
      double precision v6

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adv1 = 0
      adv2 = 0
      adv3 = 0
      adv4 = 0
      adv5 = 0
      adv6 = 0

C----------------------------------------------
C ADJOINT STATEMENTS
C----------------------------------------------
      v1 = 0.d0
      v2 = 0.d0
      v3 = 0.d0
      v4 = 0.d0
      v5 = 0.d0
      v6 = 0.d0
      do k = 1, 4
        v1 = v1+a(k,iel)*xx**(k-1)
        v4 = v4+d(k,iel)*xx**(k-1)
        if (k .lt. 4) then
          v2 = v2+b(k,iel)*xx**(k-1)
          v5 = v5+e(k,iel)*xx**(k-1)
          v3 = v3+c(k,iel)*xx**(k-1)
          v6 = v6+f(k,iel)*xx**(k-1)
        endif
      end do
      adrs = adrs+ads2*v4*(v5+2*v6*rs)
      adv4 = adv4+ads2*(1+v5*rs+v6*rs*rs)
      adv5 = adv5+ads2*v4*rs
      adv6 = adv6+ads2*v4*rs*rs
      ads2 = 0
      adrs = adrs+ads1*v1*(v2+2*v3*rs)
      adv1 = adv1+ads1*(1+v2*rs+v3*rs*rs)
      adv2 = adv2+ads1*v1*rs
      adv3 = adv3+ads1*v1*rs*rs
      ads1 = 0
      do k = 4, 1, -1
        if (k .lt. 4) then
          adxx = adxx+adv6*f(k,iel)*(k-1)*xx**(k-1-1)
          adxx = adxx+adv3*c(k,iel)*(k-1)*xx**(k-1-1)
          adxx = adxx+adv5*e(k,iel)*(k-1)*xx**(k-1-1)
          adxx = adxx+adv2*b(k,iel)*(k-1)*xx**(k-1-1)
        endif
        adxx = adxx+adv4*d(k,iel)*(k-1)*xx**(k-1-1)
        adxx = adxx+adv1*a(k,iel)*(k-1)*xx**(k-1-1)
      end do
c
      end
