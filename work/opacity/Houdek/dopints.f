c     Birational spline interplation on a regular grid
c     s/r rat2d,rbival
c     from
c     Helmuth Spaeth
c     Zweidimensionale Spline-Interpolations-Algorithmen
c     Oldenbourg Verlag Muenschen Wien 1991
c
      subroutine opints(xval,zval,tlgv,rlgv,opalg,opx,opy,
     +                  oph,opz,iexp,ier)
      implicit double precision (a-h,o-z)
c
c     Version WITH extrapolation-domain checking
c
c     History:
c     7.11.1993: creation, derived from s/r opintc
c
c     13/11/95: modified for OPAL95 tables
c
c     21/05/97: use log10(Z) instead of Z
c     23/05/97: use different algorithm for selecting tables
c               in Z and X (select min(|Z_tabs-Z_val|) + iorder tables)
c
c     17/08/97: include electron conduction (new variable lec [logical])
c
c     Last modification:
c     21/08/97
c
c  input parameters
c
c  xval ............. Hydrogen mass fraction     [  0.0 <= x  <= 1.0  ]
c  zval ............. mass fraction of heavy elements [0.0001 <=z <= 0.04 ]
c  tlgv ............. Temperature log10(t) [log10(6000)<= tlgv<= log10(100.d6)]
c  rlgv ............. log10(rho/(t6**3))         [ -7.0 <= rlg<= 1.0  ]
c  
c  output parameter
c
c  opalg ............ log10(opacity_value)
c  opx   ............ partial derivative value d(log(kappa))/d(log(rho))
c  opy   ............ partial derivative value d(log(kappa))/d(log(t))
c  oph   ............ partial derivative value d(log(kappa))/d(log(X)) 
c  opz   ............ partial derivative value d(log(kappa))/d(log(Z))
c  iexp  ............ if point lies in extrapolation area iexp > 0
c
c     setup array-dimensions
c
      parameter(ntab=7,nzva=13)
      parameter(nxir=19,nyir=70,nyis=5)
      parameter(mt=30)
c     parameter(ndat=1786)
c     parameter(ndat=(nyir-nyis+1+mt)*nxir)
c     nyii=nyir-nyis+1+mt
c     nt=2*(nxir+nyii-2)+2*(nxir*nyii-2*(nxir+nyii-2)-1)
c     parameter(nt=3348,ial=2*ndat+5)
c     parameter(nt=3348)             ! nt=3024 for OPAL95+ALEX94
c
c---  dimension for s/2 rat2d
c
      parameter(nyif=nyir+mt-nyis+1)
c
      logical ttop,tbot,rlef,rrig,lfi
      double precision xval,zval,tlgv,rlgv,opalg,opx,opy,oph,opz
      real ra
      logical*1 lec
      integer iexp
c     dimension zz(iorder+1),px(iorder+1),py(iorder+1)
c     dimension op(iorder+1),ox(iorder+1),oy(iorder+1),oh(iorder+1)
      dimension zz(7),px(7),py(7)
      dimension op(7),ox(7),oy(7),oh(7)
      dimension xshep(7),yshep(7),zzshep(7),pxshep(7),pyshep(7)
c
c---- values from s/r opinit used as argument for s/r  masubi
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyifi,mdi,
     +                nti,iali
c
c---- birational splines for s/r rat2d and rbivpd
c
      common /birasp/ pp,tls(nyif),
     +                ra(nxir,nyif,4,4,ntab,nzva)
c
c---- common /opadat/ and
c---- common /valdat/ are used for 
c---- extrapolation domain checking
c
      common /opadat/ opa(nxir,nyir,ntab,nzva),
     +                rlg(nxir,ntab,nzva),
     +                tlg(nyir,ntab,nzva)
      common /valdat/ ivalo(nyir,ntab,nzva)
c
c---- common /tablex/ & common /tablez/ are
c---- used for univariate interpolation
c---- in X and Z respectively
c
      common /tablex/ xtab(ntab),iorder,lec
      common /tablez/ ztab(nzva),ztabl(nzva)
c
c---- common for array-consistency-check of common-arrays
c---- between s/r opinit and s/r opintc
d     common /arrdim/ nxiv,nyiv,ndatv,ntabv,nzvav,niwkv,nwkv
c
      data lfi /.true./
c
      if(zval.lt.1.0d-7)then
         print '(/,x,a)','***WARNING*** opints: setting Z=1.d-7'
         zval=1.0d-7
         zvall=-7.0d0
      else
         zvall=dlog10(zval)
      endif
c      print*, zval, xval
c---- common array-dim consistency-check
d     if(nxir.ne.nxiv)goto 9710
d     if(nyir.ne.nyiv)goto 9720
d     if(ndat.ne.ndatv)goto 9730
d     if(ntab.ne.ntabv)goto 9740
d     if(nzva.ne.nzvav)goto 9750
d     if(niwkv.ne.niwkl)goto 9760
d     if(nwkv.ne.nwkl)goto 9770
c---- reset iexp
      iexp=0
c---- define ln -> log10
      um = log(10.d0)
c
c---- iorder = order of the interpolation-polynom in x and z
c---- (hydrogen mass fraction and mass fraction of heavy element)
c---- -> iorder+1 opacity-tables for interpolation in x and z
c
c---- evaluate Z-domain (mass fraction of heavy element):
c---- search in z-table ztab (bisection-method)
      jlo=1
c      print*,'iorder+1', iorder+1
      call hunt(ztab,nzvai,zval,jlo)
      if(jlo.lt.iorder)then
          nzvab=1
          nzvae=nzvab+iorder
      elseif(jlo.gt.(nzvai-(iorder-1)))then
          nzvab=nzvai-iorder
          nzvae=nzvab+iorder
      else
          if(dabs(zval-ztab(jlo)).gt.dabs(zval-ztab(jlo+1)))jlo=jlo+1
          nzvab=jlo-int(float(iorder)/2.d0+0.5d0)
          nzvae=nzvab+iorder
      endif
c
c---- evaluate X-domain (hydrogen mass fraction):
c---- search in x-table xtab (bisection-method)
      jlo=1
      call hunt(xtab,ntabi,xval,jlo)
      if(jlo.lt.iorder)then
          ntabb=1
          ntabe=ntabb+iorder
      elseif(jlo.gt.(ntabi-(iorder-1)))then
          ntabb=ntabi-iorder
          ntabe=ntabb+iorder
      else
          if(dabs(xval-xtab(jlo)).gt.dabs(xval-xtab(jlo+1)))jlo=jlo+1
          ntabb=jlo-int(float(iorder)/2.d0+0.5d0)
          ntabe=ntabb+iorder
      endif
c
      if(lfi)then
         lfi=.false.
c         print '(/,x,a,i1,a,f8.4,a,f8.4)',
c     .'opints: using ',iorder+1,' ztabs ',ztab(nzvab),' ->',ztab(nzvae)
c         print '(  x,a,i1,a,f8.4,a,f8.4,/)',
c     .'opints: using ',iorder+1,' xtabs ',xtab(ntabb),' ->',xtab(ntabe)
      endif
c
c---- interpolate at xi=rlgv ,yi=tlgv to get z=opalg
      xi=rlgv
      yi=tlgv
      jlo = 1
      ier=0
c
      do 4010 l=nzvab,nzvae
        do 4001 itab=ntabb,ntabe
c          check if point lies in the extrapolation area
           call hunt(tlg(1,itab,l),nyiri,yi,jlo)
           if(jlo.ge.nyiri)jlo=nyiri-1
           if(xi.gt.rlg(ivalo(jlo+1,itab,l),itab,l)) iexp=iexp+1
           ier=0
           call rbivpd(xi,yi,nxiri,nxir,nyifi,nyif,
     +                 rlg(1,itab,l),tls,
     +                 ra(1,1,1,1,itab,l),pp,
     +                 zi,pdx,pdy,ier)
           if(ier.gt.0)goto 9810
           zz(itab-ntabb+1)=zi
           px(itab-ntabb+1)=pdx
           py(itab-ntabb+1)=pdy
 4001   continue
c
c----   now interpolate between the iorder+1 opacity-tables
c----   for the given X-value xval using Akima-Interpolation
        call uvip3d(3,iorder+1,xtab(ntabb),zz,1,xval,op(l-nzvab+1),
     +                                               oh(l-nzvab+1))
        call uvip3p(3,iorder+1,xtab(ntabb),px,1,xval,ox(l-nzvab+1))
        call uvip3p(3,iorder+1,xtab(ntabb),py,1,xval,oy(l-nzvab+1))
 4010 continue
c---- now interpolate between the iorder+1 opacity-tables
c---- for the given Z-value zval using Akima-Interpolation
      call uvip3d(3,iorder+1,ztabl(nzvab),op,1,zvall,opalg,opz)
      call uvip3p(3,iorder+1,ztabl(nzvab),ox,1,zvall,opx)
      call uvip3p(3,iorder+1,ztabl(nzvab),oy,1,zvall,opy)
      call uvip3p(3,iorder+1,ztabl(nzvab),oh,1,zvall,oph)
c
c
c---- transform oph -> d(ln(cap))/d(ln(h))
c               opz -> d(ln(cap))/d(ln(z))
c
c     opz = opz * zval * um    ! we interpolate in log10(Z) already
      oph = oph * xval * um
c
c---- transform opx -> d(log(cap))/d(log(rho)) = capr(rho,t)
c               opy -> d(log(cap))/d(log(t))   = capt(rho,t)
c 
c     capr(rho,t)=opx
c     capt(rho,t)=opy - 3 *opx
      opy        =opy-3.d0*opx
c
c     the transformation from capx(rho,t) -> capx(p,t)
c     requires 
c               alpha ......d(ln(rho))/d(ln(p))
c             - delta ......d(ln(rho))/d(ln(t))
c               
c     rho being the density,
c     hence it has to be performed in the calling program
c
c     capp(p,t)=capr(rho,t)*  alpha
c     capt(p,t)=capr(rho,t)*(-delta) + capt(rho,t)
c
c     electron conduction ?
c     ---------------------
      if(lec)then             
         rholgv = rlgv + 3.0d0*tlgv - 1.8d01      ! log(R) -> log(rho)
         call itohec(xval,zval,tlgv,rholgv,ecolg)
         if((ecolg-opalg).lt.3.0d0)then
            opa10 = 1.0d1**opalg
            eco10 = 1.0d1**ecolg
            ope10 = opa10 + eco10
            opalg = opalg + ecolg - dlog10(ope10) ! opa = opa//eco
            ecox = 0.0d0
            ecoz = 0.0d0
            ecot = 0.0d0
            ecor = 0.0d0
            ecoe = 1.0d0
            call aditohec(xval,zval,tlgv,rholgv,
     .                    ecox,ecoz,ecot,ecor,ecoe)
            ecox = ecox * xval * um   !dlog(eco)/dX->dln(eco)/dln(X)
            ecoz = ecoz * zval * um   !dlog(eco)/dZ->dln(eco)/dln(Z)
            oph  = oph*eco10/ope10 + ecox*opa10/ope10
            opz  = opz*eco10/ope10 + ecoz*opa10/ope10
            opy  = opy*eco10/ope10 + ecot*opa10/ope10
            opx  = opx*eco10/ope10 + ecor*opa10/ope10
         endif
      endif
c
      return
c
c
d9710 print *,'opints: error in common-array nxir'
d     ier=2
d     return
d9720 print *,'opints: error in common-array nyir'
d     ier=2
d     return
d9730 print *,'opints: error in common-array ndat'
d     ier=2
d     return
d9740 print *,'opints: error in common-array ntab'
d     ier=2
d     return
d9750 print *,'opints: error in common-array nzva'
d     ier=2
d     return
d9760 print *,'opints: error in common-array niwkl'
d     ier=2
d     return
d9770 print *,'opints: error in common-array nwkl'
d     ier=2
d     return
c
c     extrapolate desired point
c     estimate opacity by extrapolation using Shepard's method
c
 9810 print *,'opints: **Warning** interpolation point outside',
     +'table domain, apply extrapolation (Shepard''s method)'
c
      ier=1
c
c     first check at which side of the table lies the desired point
c
      nexp=3
      rlef=.false.
      rrig=.false.
      ttop=.false.
      tbot=.false.
c
      if(xi.lt.rlg(1,6,10))     rlef =.true.
      if(xi.gt.rlg(nxiri,6,10)) rrig =.true.
      if(yi.lt.tlg(1,6,10))     ttop =.true.
      if(yi.gt.tlg(nyiri,6,10)) tbot =.true.
c
c     point lies outside the temperature range
c
      if(tbot.or.ttop)then
c----   search in rlg (bisection-method)
        jlo=1
        call hunt(rlg(1,6,10),nxiri,xi,jlo)
        if(jlo.lt.nexp)then
          nrlgb=1
c         nrlge=nrlgb+nexp
        elseif(jlo.gt.(nxiri-(nexp-1)))then
          nrlgb=nxiri-nexp
c         nrlge=nrlgb+nexp
        else
          nrlgb=jlo-1
c         nrlge=nrlgb+nexp
        endif
      endif
c             
c     point lies outside the density (R) range
c
      if(rlef.or.rrig)then
c----   search in tlg (bisection-method)
        jlo=1
        call hunt(tls,nyifi,yi,jlo)
        if(jlo.lt.nexp)then
          ntlgb=1
c         ntlge=ntlgb+nexp
        elseif(jlo.gt.(nyifi-(nexp-1)))then
          ntlgb=nyifi-nexp
c         ntlge=ntlgb+nexp
        else
          ntlgb=jlo-1
c         ntlge=ntlgb+nexp
        endif
       endif              
c
c     now execute the extrapolation
c
      do 5010 l=nzvab,nzvae
        do 5001 itab=ntabb,ntabe
          do 5020 ii=1,nexp+1
             if(ttop.or.tbot)then
               xii=rlg(nrlgb+ii-1,itab,l)
               if(ttop)yii=tls(1)
               if(tbot)yii=tls(nyifi)
             endif
             if(rlef.or.rrig)then
               yii=tls(ntlgb+ii-1)
               if(rlef)xii=rlg(1,itab,l)
               if(rrig)xii=rlg(nxiri,itab,l)
             endif
             ier=0
             call rbivpd(xii,yii,nxiri,nxir,nyifi,nyif,
     +                 rlg(1,itab,l),tls,
     +                 ra(1,1,1,1,itab,l),pp,
     +                 zii,pdx,pdy,ier)
             if(ier.gt.0)
     +          print *,'opints: ***Error*** in extrapolation'
             zzshep(ii)=zii
             pxshep(ii)=pdx
             pyshep(ii)=pdy
             xshep(ii) =xii
             yshep(ii) =yii
 5020     continue
          zz(itab-ntabb+1)=shep(nexp+1,xshep,yshep,zzshep,xi,yi)
          px(itab-ntabb+1)=shep(nexp+1,xshep,yshep,pxshep,xi,yi)
          py(itab-ntabb+1)=shep(nexp+1,xshep,yshep,pyshep,xi,yi)
 5001   continue
c
c----   now interpolate between the iorder+1 opacity-tables
c----   for the given X-value xval using Akima-Interpolation
        call uvip3d(3,iorder+1,xtab(ntabb),zz,1,xval,op(l-nzvab+1),
     +                                               oh(l-nzvab+1))
        call uvip3p(3,iorder+1,xtab(ntabb),px,1,xval,ox(l-nzvab+1))
        call uvip3p(3,iorder+1,xtab(ntabb),py,1,xval,oy(l-nzvab+1))
 5010 continue
c---- now interpolate between the iorder+1 opacity-tables
c---- for the given Z-value zval using Akima-Interpolation
      call uvip3d(3,iorder+1,ztabl(nzvab),op,1,zvall,opalg,opz)
      call uvip3p(3,iorder+1,ztabl(nzvab),ox,1,zvall,opx)
      call uvip3p(3,iorder+1,ztabl(nzvab),oy,1,zvall,opy)
      call uvip3p(3,iorder+1,ztabl(nzvab),oh,1,zvall,oph)
c
c     transform the results
c
c     opz = opz * zval * um  ! we interpolate in log10(Z) already 
      oph = oph * xval * um
      opy = opy-3.d0*opx
c
c     electron conduction ?
c     ---------------------
      if(lec)then             
         rholgv = rlgv + 3.0d0*tlgv - 1.8d01      ! log(R) -> log(rho)
         call itohec(xval,zval,tlgv,rholgv,ecolg)
         if((ecolg-opalg).lt.3.0d0)then
            opa10 = 1.0d1**opalg
            eco10 = 1.0d1**ecolg
            ope10 = opa10 + eco10
            opalg = opalg + ecolg - dlog10(ope10) ! opa = opa//eco
            ecox = 0.0d0
            ecoz = 0.0d0
            ecot = 0.0d0
            ecor = 0.0d0
            ecoe = 1.0d0
            call aditohec(xval,zval,tlgv,rholgv,
     .                    ecox,ecoz,ecot,ecor,ecoe)
            ecox = ecox * xval * um   !dlog(eco)/dX->dln(eco)/dln(X)
            ecoz = ecoz * zval * um   !dlog(eco)/dZ->dln(eco)/dln(Z)
            oph  = oph*eco10/ope10 + ecox*opa10/ope10
            opz  = opz*eco10/ope10 + ecoz*opa10/ope10
            opy  = opy*eco10/ope10 + ecot*opa10/ope10
            opx  = opx*eco10/ope10 + ecor*opa10/ope10
         endif
      endif
c
      return
c
      end
