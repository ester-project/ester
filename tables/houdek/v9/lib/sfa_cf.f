      subroutine sfa( rs, xx, iel, s1, s2 )
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
      double precision rs
      double precision s1
      double precision s2
      double precision xx

C==============================================
C define local variables
C==============================================
      integer k
      double precision v1
      double precision v2
      double precision v3
      double precision v4
      double precision v5
      double precision v6

C**********************************************
C executable statements of routine
C**********************************************
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
      s1 = v1*(1.d0+v2*rs+v3*rs*rs)
      s2 = v4*(1.d0+v5*rs+v6*rs*rs)
c
      end


