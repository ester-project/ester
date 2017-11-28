      subroutine adcondux( ro6, t8, zz, iel, adro6, adt8, adcond )
      implicit none

C==============================================
C define parameters                            
C==============================================

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      integer iel
      double precision adcond
      double precision adro6
      double precision adt8
      double precision ro6
      double precision t8
      double precision zz

C==============================================
C define local variables
C==============================================
      double precision a
      double precision adfact
      double precision adgam
      double precision adrs
      double precision ads
      double precision ads1
      double precision ads2
      double precision adxx
      double precision fact
      double precision gam
      double precision rs
      double precision s
      double precision s1
      double precision s2
      double precision xx

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adfact = 0
      adgam = 0
      adrs = 0
      ads = 0
      ads1 = 0
      ads2 = 0
      adxx = 0

C----------------------------------------------
C ADJOINT STATEMENTS
C----------------------------------------------
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
      adfact = adfact-adcond*(2.363d+17*ro6*t8/a*s/((1.d0+fact)*s*(1.d0+
     $fact)*s))
      adro6 = adro6+adcond*(2.363d+17*(t8/a)/((1.d0+fact)*s))
      ads = ads-adcond*(2.363d+17*ro6*t8/a*(1.d0+fact)/((1.d0+fact)*s*
     $(1.d0+fact)*s))
      adt8 = adt8+adcond*(2.363d+17*(ro6/a)/((1.d0+fact)*s))
      adcond = 0
      s = s1-fact/(1.d0+fact)*s2
      if (s .lt. 0.1d0) then
        ads = 0
      endif
      adfact = adfact-ads*(1/(1.d0+fact)-fact/((1.d0+fact)*(1.d0+fact)))
     $*s2
      ads1 = ads1+ads
      ads2 = ads2-ads*(fact/(1.d0+fact))
      ads = 0
      adro6 = adro6+0.66666666666667d0*1.018d0*adfact*(zz/a)**
     $0.66666666666667d0*ro6**(-0.33333333333333d0)
      adfact = 0
      call adsfa( rs,xx,iel,adrs,adxx,ads1,ads2 )
      rs = 0.01388d0*(a/zz)**(1.d0/3.d0)/ro6**(1.d0/3.d0)
      if (rs .gt. 1.0d1) then
        adrs = 0
      endif
      adro6 = adro6-adrs*(0.33333333333333d0*0.01388d0*(a/zz)**
     $0.33333333333333d0*ro6**(-0.66666666666667d0)/(ro6**
     $0.33333333333333d0*ro6**0.33333333333333d0))
      adrs = 0
      adgam = adgam+0.45641d0*adxx*(1.d0/gam)
      adxx = 0
      adro6 = adro6+0.33333333333333d0*adgam*0.2275d0*zz**2/t8/a*(ro6/a)
     $**(-0.66666666666667d0)
      adt8 = adt8-adgam*0.2275d0*zz**2/(t8*t8)*(ro6/a)**
     $0.33333333333333d0
      adgam = 0

      end


