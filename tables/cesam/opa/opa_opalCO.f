
c***********************************************************************

	SUBROUTINE opa_opalCO(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)

c routine private du module mod_opa

c Routine de calcul de l'opacité basée sur les tables d'opacité
c de type 1 de Livermore.
c il n'y a pas interpolation en Z, mais correction pour C et O par
c rapport au Z de la table OPAL utilisée

c Auteurs : initialisation en F77 de L.Piau Univ. Bruxelles
c P.Morel, Département J.D. Cassini, O.C.A., CESAM2k	

c entrées :
c	xchim: comp. chim. en H + C & O
c	t6 : température millions de K
c	r  : densité cgs/(t6**3)

c sorties :
c	kapa : opacité gr / cm2
c	dkapdt : kappa / d T
c	dkapdr : kappa / d densite              
c       dkapdx : kappa / d X

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY: ab_min, nchim, nom_elem
	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xh
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out) :: dkapdt, dkapdr, dkapdx, kappa

	INTEGER, PARAMETER :: ipr=20, mx=5, mc=8, mo=8, nrm=19, nrb=1,
	1 nre=19, nr=nre+1-nrb, ntabs=60, ntm=70, ntb=1, nt=ntm+1-ntb

	REAL (kind=dp), SAVE, DIMENSION(mx,mc,mo,nt,nr) :: co=1.d35
	REAL (kind=dp), SAVE, DIMENSION(mx,mc,nt,nr) :: diag=1.d35
	REAL (kind=dp), SAVE, DIMENSION(mx,mo,nt,nr) :: diago=1.d35	
	REAL (kind=dp), SAVE, DIMENSION(mx,nt,nr) :: opl	
	REAL (kind=dp), SAVE, DIMENSION(mx,mc) :: cxdf, cxf, oxdf, oxf,
	1 xcdf, xodf
	REAL (kind=dp), SAVE, DIMENSION(mx,ntabs) :: x, xca, xoa, y, zz
	REAL (kind=dp), SAVE, DIMENSION(nt,nr) :: opk, opk2
	REAL (kind=dp), SAVE, DIMENSION(ntm,nrm) :: cof
	REAL (kind=dp), SAVE, DIMENSION(100,nr) :: coff
	REAL (kind=dp), SAVE, DIMENSION(85,ipr) :: f, fx, fy, fxy
	REAL (kind=dp), SAVE, DIMENSION(3,mx) :: a		
	REAL (kind=dp), SAVE, DIMENSION(mc) :: cx, cxd, xc, xcd,
	1 xcs=(/ 0.d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.d0 /), xod
	REAL (kind=dp), SAVE, DIMENSION(mo) :: ox, oxd, xo,
	1 xos=(/ 0.d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,0.6d0,1.d0 /)
	REAL (kind=dp), SAVE, DIMENSION(mx) :: dfsx, xx	
	REAL (kind=dp), SAVE, DIMENSION(nr) :: alr, dfsr
	REAL (kind=dp), SAVE, DIMENSION(nrm) :: alrf 		
	REAL (kind=dp), SAVE, DIMENSION(nt) :: alt, dfs, t6list  
	REAL (kind=dp), SAVE, DIMENSION(5) :: xa=(/ 0.d0,0.03d0,0.1d0,
	1 0.35d0,0.7d0 /)	
	REAL (kind=dp), SAVE, DIMENSION(100) :: t6arr
	REAL (kind=dp), SAVE, DIMENSION(4) :: q, h
	REAL (kind=dp), SAVE, DIMENSION(3) :: b
	REAL (kind=dp), PARAMETER :: dx=1.d-05, unpdx=1.d+00+dx
	REAL (kind=dp), SAVE :: cxx, dkap, oxx, rle, rls, tmax, xcdp,
	1 xodp, xxc0, xxco, xxo0, z_table, zzz
	REAL (kind=dp) :: dopacr, dopact, opact, r, t6, xhh, xxc, xxo

	INTEGER, SAVE, DIMENSION(mx,mc) :: n=0
	INTEGER, SAVE, DIMENSION(mx,ntabs) :: itab
	INTEGER :: i	
	INTEGER, SAVE, DIMENSION(nrm) :: nta=(/ (70,i=1,14),69,64,60,58,
	1 57 /)
	INTEGER, SAVE, DIMENSION(101) :: index = (/ 1,2,2,3,3,3,3,3,3,3,4,
	1 4,4,4,4, 4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,
	2 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
	3 7, 7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7 /)
	
	INTEGER, SAVE, DIMENSION(5) :: isoh, isoc, isoo
	INTEGER, SAVE :: ip, iq, k1, k2, k3, k4, l1, l2, l3, l4,
	1 m, mf, mxzero, nc, no, nrl, nset, ntemp, nsm, nrlow, nrhigh,
	2 nisoh, nisoc, nisoo 

        LOGICAL, PARAMETER :: lisse=.TRUE. !On smooth .TRUE./.FALSE.
        LOGICAL, SAVE, DIMENSION(mx) :: itime=.TRUE. 	
        LOGICAL, SAVE :: init=.TRUE., sortie=.FALSE.
	
	CHARACTER (len=9), PARAMETER, DIMENSION(5) :: fich_opalCO =
	1 (/ 'COX00.tab','COX03.tab','COX10.tab','COX35.tab','COX70.tab'/)
	     
c-----------------------------------------------------------------

2000	FORMAT(8es10.3)

c appel à opa_yveline pour la zone externe et l'atmosphère T < 1eV	 
	IF(t < 1.5d4)THEN
	 CALL opa_yveline(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx) ; RETURN	 	 
	ENDIF	 

c--------------initialisations---------------------------------
	IF(init)THEN
	 init=.FALSE.
	 
c identification des indices de H, He, C, O
	 DO i=1,nchim
	  IF(nom_elem(i)(1:2) == ' H')THEN
	   nisoh=nisoh+1 ; isoh(nisoh)=i
	  ELSEIF(nom_elem(i)(1:2) == 'C1')THEN
	   nisoc=nisoc+1 ; isoc(nisoc)=i
	  ELSEIF(nom_elem(i)(1:2) == 'O1')THEN
	   nisoo=nisoo+1 ; isoo(nisoo)=i
	  ENDIF	  	 	  	
	 ENDDO
c	 PRINT*,nom_elem	 
c	 PRINT*,nisoh, nisoc, nisoo
c	 PRINT*,isoh ; PRINT*,isoc ; PRINT*,isoo
c	 PAUSE'iso'	 

c définitions de z_table, C/Z, O/Z suivant les groupes de fichiers
c d'opacité utilisés (donnés dans l'entête des tables)
c	 avec COX.. Z=0.02	 
	 z_table=0.02d0 ; xxc0=z_table*0.173285d0 ; xxo0=z_table*0.482272d0
	 
	ENDIF	
c--------------initialisations (fin)---------------------------------	

        t6=t*1.d-06 ; r=ro/t6**3

c abondances de H, He, C et O
        xhh=SUM(xh(isoh(1):isoh(nisoh))) 	!fraction totale de H
        xxc=SUM(xh(isoc(1):isoc(nisoc)))	!fraction totale de C 
        xxo=SUM(xh(isoo(1):isoo(nisoo))) 	!fraction totale de O

c xxc, xxO : excès de C et O dans Z par rapport aux C et O des tables	
	xxc=xxc-xxc0 ; xxo=xxo-xxo0
	CALL opac(z_table,xhh,xxc,xxo,t6,r) ; IF(sortie)STOP
	kappa=10.d0**opact ; dkapdt=kappa*dopact/t
	dkapdr=kappa*dopacr/ro ; dkapdx=0.d0	
	IF(xhh > ab_min(1))THEN
	 xhh=xhh*unpdx
	 CALL opac(z_table,xhh,xxc,xxo,t6,r)
	 dkapdx=(10.d0**opact-kappa)/xhh/dx
	ENDIF	 
	
	RETURN
	
	CONTAINS
	
c**********************************************************************

	SUBROUTINE opac(z_i,xh_i,xxc_i,xxo_i,t6,r)	
	
	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: z_i, xh_i, xxc_i, xxo_i, t6, r
	
	REAL (kind=dp) :: cmod, dixr, slr, slt, somme, xh, xhe, xxc, xxci,
	1 xxo, xxoi, xxx, z
	
	INTEGER :: i, iadvance, ihi, ilo, imd, ir, is, istep1, it, iw, kmin,
	1 k1in, k3s, l3s, mfin, mf2, mg, mh, mi, ntd, ntlimit
	
c----------------------------------------------------------------

	xh=xh_i ; xxc=xxc_i ; xxo=xxo_i ; z=z_i
	
c normalisation
	somme=z+xh+xxc+xxo
	IF(somme > 1.d0)THEN
	 xxc=xxc/somme ; xxo=xxo/somme ; z=z/somme
	ENDIF

c refers to goto to compute opacity composition derivative
	IF(nr < 6)THEN
	 WRITE(*,65)nre,nrb ; WRITE(2,65)nre,nrb ; STOP
65	 FORMAT('opa_opalCO: Too few R values; nre+1-nrb < 6, nre=',i3,
	1 ', nrb=',i3)
	ENDIF	
	
	IF(xh > 1.d-6 .AND. mx < 4)THEN
	 WRITE(*,64) ; WRITE(2,64) ; STOP
64	 FORMAT('opa_opalCO: xh not equal to zero: To run this case it
	1 is necessary',/,'to recompile with PARAMETER (mx=5)')
	ENDIF
	
c..... set-up C/O axis points
	xxco=xxc+xxo
	IF(z+xh+xxco-1.d-6 > 1.d0 )THEN
	 WRITE(*,61)z,xh,xxc,xxo,z+xh+xxco ; WRITE(2,61)z,xh,xxc,xxo,z+xh+xxco
	 STOP	
61	 FORMAT('opa_opalCO: Mass fractions exceed unity',/,'z=',es10.3,
	1 ', xh=',es10.3,', xxc=',es10.3,', xxo=',es10.3,', z+xh+xxco=',es10.3)
	ENDIF	
	
	zzz=z+0.001d0 ; xxci=xxc ; xxoi=xxo

c..... convert xxc and xxo to logarithmic shifted by Z
	cxx=LOG10(zzz+xxc) ; oxx=LOG10(zzz+xxo) ; xxx=LOG10(0.005d0+xh)
	slt=LOG10(t6) ; slr=LOG10(r)
 
c..... set x indices
	ilo=2 ; ihi=mx	
	B3: DO
	 IF(ihi-ilo <= 1)THEN
	  EXIT B3
	 ELSE
	  imd=(ihi+ilo)/2
          IF(xh <= xa(imd)+1.d-7)THEN
           ihi=imd
          ELSE
           ilo=imd
          ENDIF
	 ENDIF
	ENDDO B3
	i=ihi ; mf=i-2 ; mg=i-1 ; mh=i ; mi=i+1 ; mf2=mi ; istep1=1
	IF (mx > 1)THEN
	 istep1=mx-1
	 IF((xh <= xa(2)+1.d-7) .OR. (xh >= xa(istep1)-1.d-7))mf2=mh
	ENDIF
	IF (mx == 1 .OR. xh < 1.d-6)THEN
	 mf=1 ; mg=1 ; mh=1 ; mi=2 ; mf2=1
	ENDIF
	IF(itime(1))THEN 
	 alr(1)=-8.d0+(nrb-1)*0.5d0 
	 DO i=2,nr  
	  alr(i)=alr(i-1)+0.5d0  
	 ENDDO  
	 alt(1)=-2.25d0+(ntb-1)*0.05d0  
	 DO i=ntb+1,46
	  alt(i)=alt(i-1)+.05d0
	 ENDDO
	 ntd=47
	 IF(ntb+1 > 47)ntd=ntb+1
	 DO i=ntd,68
	  alt(i)=alt(i-1)+.1d0
	 ENDDO
	 DO i=68,70
	  alt(i)=alt(i-1)+.2d0
	 ENDDO
	 DO i=1,nt
	  t6list(i)=10.d0**alt(i)
	 ENDDO
	ENDIF   
	ilo=2 ; ihi=nr
	B4: DO	
	 IF(ihi-ilo <= 1)THEN
	  EXIT B4
	 ELSE
	  imd=(ihi+ilo)/2
	  IF(slr <= alr(imd)+1.d-7)THEN
	   ihi=imd
	  ELSE
	   ilo=imd
	  ENDIF
	 ENDIF
	ENDDO B4
	i=ihi ; l1=i-2 ; l2=i-1 ; l3=i ; l4=l3+1 ; ilo=2 ; ihi=nt
	
	B5: DO
	 IF(ihi-ilo <= 1)THEN	
	  EXIT B5
	 ELSE
	  imd=(ihi+ilo)/2
	  IF(t6 <= t6list(imd)+1.d-7)THEN
	   ihi=imd
	  ELSE
	   ilo=imd
	  ENDIF
	 ENDIF
	ENDDO B5
	i=ihi ; k1=i-2 ; k2=i-1 ; k3=i ; k4=k3+1 ; k3s=k3+ntb-1 
	l3s=l3+nrb-1

c-----set-up indices to skip when low T&R data missing for x=0.
	kmin=0 ; k1in=k1 ; iadvance=0 ; mfin=mf
	IF(mfin == 1 .AND. co(1,1,1,k1,l1) > 9.d0)THEN! data missing
	 DO i=1,6
	  IF(co(1,1,1,i,l1) > 9.d0)THEN
	   IF(xh < 0.1d0)THEN
	    kmin=i+1
	   ELSEIF(iadvance == 0)THEN  ! sfift x index to avoid x=0.
	    iadvance=iadvance+1 ; mf=mf+1 ; mg=mg+1 ; mh=mh+1
	    mi=mi+1 ; mf2=mf2+1
	   ENDIF
	  ENDIF
	 ENDDO
	 IF(iadvance == 0 .AND. k1 <= kmin .AND. slt <= alt(kmin))THEN
	  k1=kmin
	  IF(co(1,1,1,kmin,l1+1) < 9.d0 .AND. slr+.01d0 > alr(l1+1))THEN
	   l1=l1+1 ; kmin=0 ; k1=k1in
	   DO i=1,6
	    IF(co(1,1,1,i,l1) > 9.d0)kmin=i+1
	   ENDDO
	   IF(kmin /= 0 .AND. k1in < kmin)k1=kmin
	  ENDIF
	 ENDIF
	 IF(slt+.001d0 < alt(k1))THEN
	  WRITE(*,70)xh,slt,slr ; WRITE(2,70)xh,slt,slr
70	  FORMAT('op_opalCO: data not available for x=',es10.3,/,
	1 'logT6=',es10.3,', logR=',es10.3)
	  sortie=.TRUE. ; RETURN	
	 ENDIF
	 l2=l1+1 ; l3=l2+1 ; l4=l3+1 ; l3s=l3+nrb-1 ; k2=k1+1 ; k3=k2+1
	 k4=k3+1 ; k3s=k3+ntb-1
	ENDIF	
c-----end of check for missing data

	DO m=mf,mf2
	 IF(mx >= 4)THEN
c..... C and O  fractions determined by the ray through the origin that
c        also passes through the point (xc,xo). Specific interpolation 
c        values determined by tabulated x values;i.e. xa(m).  inter-
c        polation along the ray gives log (kappa(xc,xo)).  (Advantage
c        of method: keeps indices within table boundaries)
c        Subtract Z to prevent out-of-range C+O values for small x
	  IF(1.d0-xh-z > 1.d-6)THEN
	   cmod=(1.d0-xa(m)-z)/(1.d0-xh-z)
	  ELSE
	   cmod=0.d0
	  ENDIF
	  xxc=cmod*xxci ; xxo=cmod*xxoi ; cxx=LOG10(zzz+xxc)
	  oxx=LOG10(zzz+xxo)
	 ENDIF
c        ythism=z+xa(m)+xxc+xxo

	 B1: DO i=1,mc
          xhe=1.d0-xa(m)-z ! xhe est une variable intermediaire
	  nc=i ; no=i ; xc(i)=xcs(i) ; xo(i)=xos(i)
	  IF(xcs(i) > xhe)THEN
	   xc(i)=xhe ; xo(i)=xhe ; EXIT B1
	  ENDIF
	 ENDDO B1

	 IF(itime(m))THEN
	  itime(m)=.FALSE. ; mxzero=0
	  DO i=1,mx
	   xx(i)=LOG10(0.005d0+xa(i))
	   IF(xa(i) == 0.d0)mxzero=i
	  ENDDO
	  
c.....    this is the first time throught this m. Calculate the decadic
c         log of the perimeter points shifted by Z+0.001(to avoid diver. 
c         at origin); m refers to xa(m); the hydrogen table value.
c         note that the nc-th elements are sometimes needed!
	  DO i=1,nc
           oxf(m,i)=LOG10(zzz+xo(i)) ; cxf(m,i)=LOG10(zzz+xc(i))
	   xcdf(m,i)=-xo(i)+xo(no) ; xodf(m,i)=-xc(i)+xc(nc)
	   cxdf(m,i)=LOG10(zzz+xcdf(m,i)) ; oxdf(m,i)=LOG10(zzz+xodf(m,i))
	  ENDDO

c         note that the nc-th elements are sometimes needed!
	  DO i=1,nc
	   ox(i)=oxf(m,i) ; cx(i)=cxf(m,i) ; xcd(i)=xcdf(m,i)
	   xod(i)=xodf(m,i) ; cxd(i)=cxdf(m,i) ; oxd(i)=oxdf(m,i)
	  ENDDO

c.....    READ the data files
	  CALL readco
	 ENDIF
	 DO i=1,nc
	  ox(i)=oxf(m,i) ; cx(i)=cxf(m,i) ; xcd(i)=xcdf(m,i)
	  xod(i)=xodf(m,i) ; cxd(i)=cxdf(m,i) ; oxd(i)=oxdf(m,i)
	 ENDDO

c.....   Determine log R and log T6 grid points to use in the
c        interpolation.
	 IF(slt < alt(1) .OR. slt > alt(nt))THEN
	  WRITE(*,62)slt,alt(1),alt(nt) ; WRITE(2,62)slt,alt(1),alt(nt)
	  sortie=.TRUE. ; RETURN
62	  FORMAT('opa_opalCO: slt outside of table range',/,'slt=',es10.3,
	1  ' en dehors de: [',es10.3,',',es10.3,']')	
	 ENDIF
	 	 
	 IF(slr < alr(1) .OR. slr > alr(nr))THEN
	  WRITE(*,60)slr,alr(1),alr(nr) ; WRITE(2,60)slr,alr(1),alr(nr)
60	  FORMAT('opa_opalCO: slr outside of table range',/,'slr=',es10.3,
	1 ' en dehors de: [',es10.3,',',es10.3,']')
	  sortie=.TRUE. ; RETURN	
	 ENDIF	 
	 
	 IF (m == mf)THEN  !  calculate table indices
	  IF(mf2 /= mxzero .AND. k3s > ntm)THEN
	   WRITE(*,67)mf2,mxzero,k3s,ntm ; WRITE(2,67)mf2,mxzero,k3s,ntm
67	   FORMAT('opa_opalCO: mf2 /= mxzero avec k3s > ntm',/,
	1  'mf2=',i3,', mxzer=',i3,', k3s=',i3,', ntm=',i3)   	  
	   sortie=.TRUE. ; RETURN
	  ENDIF 
	 ENDIF	 
	   
	 DO i=14,18
          IF(l3s > i .AND. k3s > nta(i+1))THEN
	   WRITE(*,68)l3s,i,k3s,nta(i+1) ; WRITE(2,68)l3s,i,k3s,nta(i+1)
68	   FORMAT('opa_opalCO: l3s > i .AND. k3s > nta(i+1)',/,
	1  'l3s=',i4,', i=',i3,', k3s=',i3,', nta(i+1)=',i3)	    	   
	   sortie=.TRUE. ; RETURN	
	  ENDIF
	 ENDDO	 	 
	  
	 ip=3 ; iq=3 ; ntlimit=nta(l3s)
	 IF(k3s == ntlimit .OR. .NOT.lisse)THEN
          ip=2 ; iq=2
	 ENDIF
	 IF(t6 <= t6list(2)+1.d-7 .OR. .NOT.lisse)ip=2
	 IF(l3 == nr .OR. .NOT.lisse)THEN ! right edge of full table
	  iq=2 ; ip=2
	 ENDIF
	 IF(slr <= alr(2)+1.d-7 .OR. .NOT.lisse)iq=2
	 xodp=MAX(-xxc+xc(nc),0.d0) ; xcdp=MAX(-xxo+xo(no),0.d0) ; is=0
	 CALL cointerp(xxc,xxo)
	ENDDO

	IF(zz(mg,1)/= zz(mf,1) .OR. zz(mh,1) /= zz(mf,1)
	1 .OR. z /= zz(mf,1))THEN
	 WRITE(*,66)z ; WRITE(2,66)z ; STOP 
66	 FORMAT('opa_opalCO: Z does not match Z in files used ,z=',es10.3)
	ENDIF
	xxc=xxci   ! restores input value; necessary IF STOP replaced 
c                  with RETURN
	xxo=xxoi   ! restores input value
	is=0 ; iw=1
	DO ir=l1,l1+iq
	 B2: DO it=k1,k1+ip
	  IF(mx == 1 .OR. mf2 == 1)THEN
	   opk(it,ir)=opl(mf,it,ir) ; CYCLE B2
          ENDIF
	  opk(it,ir)=quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir),
	1 opl(mh,it,ir),xx(mf),xx(mg),xx(mh))
          is=1
	 ENDDO B2
	ENDDO

	IF(mi == mf2)THEN  ! interpolate between quadratics
	 is=0 ; iw=1 ; dixr=(xx(mh)-xxx)*dfsx(mh)
	 DO ir=l1,l1+iq
	  DO it=k1,k1+ip
	   opk2(it,ir)=quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir),
	1  opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
	   opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.-dixr) ; is=1
          ENDDO
	 ENDDO
c        interpolate x between two overlapping quadratics
	ENDIF
	is=0

c.....  completed H,C,O interpolation. Now interpolate T6 and log R on a
c       4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure
c       mixes overlapping quadratics to obtain smoothed derivatives.

	CALL t6rinterp(slr,slt)
	
	RETURN

	END SUBROUTINE opac

c************************************************************************

	SUBROUTINE cointerp(xxc,xxo)
      
c The purpose of this SUBROUTINE is to interpolate in C and O abundances.

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: xxc, xxo
	REAL (kind=dp) :: cxdp, fac, oxdp

	INTEGER :: i, ie, iei, iej, ir, is, it, ix, iw, i1, i2, i3,
	1 jx, j1, j2, j3, m1, m2, w
       
c--------------------------------------------------------------

	is=0       
	IF(xxco < 1.d-6)THEN
	 DO ir=l1,l1+iq
	  DO it=k1,k1+ip
	   opl(m,it,ir)=co(m,1,1,it,ir)
	  ENDDO
	 ENDDO
	 RETURN
	ENDIF
	
c	include boundaries that could later cause division by 0!
	IF(xxc > xcd(3)-1.d-6)THEN
	 oxdp=LOG10(zzz+xodp)
c	 handle possibility that xodp=0
	 fac=MAX(MIN((oxx-ox(1))/MAX(oxdp-ox(1),1.d-6),1.d0),0.d0)
	 DO ir=l1,l1+iq
	  B1: DO it=k1,k1+ip
c	   interpolation in region c1
c	   include boundaries that could later cause division by 0!
	   IF(xxc > xcd(2)-1.d-6)THEN
	    iw=1
	    a(1,m)=quad(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir),
	1   diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
	    iw=iw+1
	    a(2,m)=quad(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir),
	1   diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
	    DO w=1,2
	     b(w)=a(w,m)
	    ENDDO
c	    handle possibility that xodp=0
	    opl(m,it,ir)=b(1)+(b(2)-b(1))*fac ; is=1 ; CYCLE B1
	   ENDIF
c          interpolation in region c2
	   iw=1
	   a(1,m)=quad(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir),
	1  diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
	   iw=iw+1
	   a(2,m)=quad(is,iw,cxx,co(m,n(m,2)-2,2,it,ir),co(m,n(m,2)-1,2,
	1  it,ir),diag(m,2,it,ir),cx(n(m,2)-2),cx(n(m,2)-1),cxd(2))
	   iw=iw+1
	   a(3,m)=quad(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir),
	1  diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
	   DO w=1,3
	    b(w)=a(w,m)
	   ENDDO
	   iw=iw+1
	   opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(1),ox(2),oxdp)
	   is=1
	  ENDDO B1
	 ENDDO
	 IF(is == 1)RETURN
	ENDIF

c interpolation in region c3 to c6
	is=0
	IF(nc >= 5)THEN
	 DO i=4,nc-1
c         DO not go beyond middle (where c3-c6 overlaps o3-o6), and
	  IF(xxc > xcd(i)-1.d-6 .AND. xxo > xo(i-1)-1.d-6 .AND.
	1  xcd(i-1) > xc(i-1))THEN
	   DO ir=l1,l1+iq
	    DO it=k1,k1+ip
	     oxdp=LOG10(zzz+xodp) ; iw=1 ; m1=i-1 ; m2=i-2
	     a(1,m)=quad(is,iw,cxx,co(m,n(m,m2)-2,m2,it,ir),
	1    co(m,n(m,m2)-1,m2,it,ir),diag(m,m2,it,ir),cx(n(m,m2)-2),
	2    cx(n(m,m2)-1),cxd(m2))
	     iw=iw+1
	     a(2,m)=quad(is,iw,cxx,co(m,n(m,m1)-2,m1,it,ir),
	1    co(m,n(m,m1)-1,m1,it,ir),diag(m,m1,it,ir),cx(n(m,m1)-2),
	2    cx(n(m,m1)-1),cxd(m1))
	     iw=iw+1
	     a(3,m)=quad(is,iw,cxx,diag(m,m2,it,ir),diag(m,m1,it,ir),
	1    diag(m,i,it,ir),cxd(m2),cxd(m1),cxd(i))
	     DO w=1,3
	      b(w)=a(w,m)
	     ENDDO
	     iw=iw+1
	     opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(i-2),
	1    ox(i-1),oxdp)
	     is=1
	    ENDDO
	   ENDDO
	   IF(is == 1)RETURN
	  ENDIF
	 ENDDO
	ENDIF
	IF(is == 1)RETURN

c	include boundaries that could later cause division by 0!
	IF(xxo > xod(3)-1.d-6)THEN
	 cxdp=LOG10(zzz+xcdp)
c	 handle possibility that xcdp=0
	 fac=MAX(MIN((cxx-cx(1))/MAX(cxdp-cx(1),1.d-6),1.),0.)
	 DO ir=l1,l1+iq
	  B2: DO it=k1,k1+ip
c	   interpolation in region  o1
c	   include boundaries that could later cause division by 0!
	   IF(xxo > xod(2)-1.d-6)THEN
	    iw=1
	    a(1,m)=quad(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir),
	1   diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
	    iw=iw+1
	    a(2,m)=quad(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir),
	1   diago(m,no-3,it,ir),oxd(1),oxd(2),oxd(3))
	    DO w=1,2
	     b(w)=a(w,m)
	    ENDDO
c	    handle possibility that xcdp=0
	    opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
	    is=1 ; CYCLE B2
	   ENDIF
c	   interpolation in region  o2

	   iw=1
	   a(1,m)=quad(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir),
	1  diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
	   iw=iw+1
	   a(2,m)=quad(is,iw,oxx,co(m,2,n(m,2)-2,it,ir),
	1  co(m,2,n(m,2)-1,it,ir),diago(m,no-2,it,ir),ox(n(m,2)-2),
	2  ox(n(m,2)-1),oxd(2))
	   iw=iw+1
	   a(3,m)=quad(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir),
	1  diago(m,nc-3,it,ir),oxd(1),oxd(2),oxd(3))
	   DO w=1,3
	    b(w)=a(w,m)
	   ENDDO
	   iw=iw+1
	   opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(1),cx(2),cxdp)
	   is=1
	  ENDDO B2
	 ENDDO
	 IF(is == 1)RETURN
	ENDIF
c	interpolation in region  o3 to o6
	is=0
	IF(no >= 5)THEN
	 DO i=4,no-1
c         DO not go beyond middle (where o3-o6 overlaps c3-c6), and
	  IF(xxo > xod(i)-1.d-6 .AND. xxc > xc(i-1)-1.d-6 .AND.
	1  xod(i-1) > xo(i-1)-1.d-6)THEN
	   DO ir=l1,l1+iq
	    DO it=k1,k1+ip
	     cxdp=LOG10(zzz+xcdp) ; iw=1 ; m2=i-2 ; m1=i-1
	     a(1,m)=quad(is,iw,oxx,co(m,m2,n(m,m2)-2,it,ir),
	1    co(m,m2,n(m,m2)-1,it,ir),diago(m,no-m2,it,ir),ox(n(m,m2)-2),
	2    ox(n(m,m2)-1),oxd(m2))
	     iw=iw+1
	     a(2,m)=quad(is,iw,oxx,co(m,m1,n(m,m1)-2,it,ir),
	1    co(m,m1,n(m,m1)-1,it,ir),diago(m,no-m1,it,ir),ox(n(m,m1)-2),
	2    ox(n(m,m1)-1),oxd(m1))
	     iw=iw+1
	     a(3,m)=quad(is,iw,oxx,diago(m,no-m2,it,ir),
	1    diago(m,no-m1,it,ir),diago(m,no-i,it,ir),oxd(m2),oxd(m1),
	2    oxd(i))
	     DO w=1,3
	      b(w)=a(w,m)
	     ENDDO
	     iw=iw+1
	     opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(m2),cx(m1),
	1    cxdp)
	     is=1
	    ENDDO
	   ENDDO
	   IF(is == 1)RETURN
	  ENDIF
	 ENDDO
	ENDIF
	IF(is == 1)RETURN

c.....find index of C grid.

	ie=100*NINT(xxc)+1 ; iei=index(ie)+1	
	
c	must also allow index = nc, to avoid extrapolation
	IF(iei > nc)iei=nc
	IF(iei > 3)THEN
	 i1=iei-2 ; i2=iei-1 ; i3=iei
	ELSE
	 i1=1 ; i2=2 ; i3=3
	ENDIF

c.....find index of O grid

	ie=100*NINT(xxo)+1 ; iej=index(ie)+1
c       must also allow index = no, to avoid extrapolation
	IF(iej > no)iej=no
	IF(iej > 3)THEN
	 j1=iej-2 ; j2=iej-1 ; j3=iej
	ELSE
	 j1=1 ; j2=2 ; j3=3
	ENDIF
	
c	lower-O part of grid: interpolate C before O
	IF(j3 < no .AND. i3 <= n(m,j3) .AND.
	1 (xxc < xcd(j3)+1.d-6 .OR. xxc >= xxo))THEN
	 DO ir=l1,l1+iq
	  DO it=k1,k1+ip
	   iw=0
	   DO jx=j1,j1+2
	    iw=iw+1
	    a(iw,m)=quad(is,iw,cxx,co(m,i1,jx,it,ir),co(m,i2,jx,it,ir),
	1   co(m,i3,jx,it,ir),cx(i1),cx(i2),MIN(cx(i3),cxd(jx)))
	   ENDDO
	   DO w=1,3
	    b(w)=a(w,m)
	   ENDDO
	   iw=iw+1
	   opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(j1),ox(j2),ox(j3))
	   is=1
	  ENDDO
	 ENDDO
c        ELSE: high-O part of grid: must interpolate O before C
	ELSE
	 DO ir=l1,l1+iq
	  DO it=k1,k1+ip
	   iw=0
	   DO ix=i1,i1+2
	    iw=iw+1
	    IF(j3 < n(m,ix))THEN
	     a(iw,m)=quad(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir),
	1    co(m,ix,j3,it,ir),ox(j1),ox(j2),ox(j3))
	    ELSE
	     a(iw,m)=quad(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir),
	1    diago(m,no-ix,it,ir),ox(j1),ox(j2),oxd(ix))
	    ENDIF
	   ENDDO
	   DO w=1,3
	    b(w)=a(w,m)
	   ENDDO
	   iw=iw+1
	   opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(i1),cx(i2),cx(i3))
	   is=1
	  ENDDO
	 ENDDO
	ENDIF
  
	RETURN
      
	END SUBROUTINE cointerp

c***********************************************************************

	SUBROUTINE t6rinterp(slr,slt)

c	The purpose of this SUBROUTINE is to interpolate in logT6 and logR

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL(kind=dp), INTENT(in) :: slr, slt

	REAL(kind=dp) :: dix, dix2, dkapq1, dkapq2, dkap1, dkap2, dopacrq,
	1 dopactq, opactq, opactq2, opact2
	
	INTEGER :: is, iu, iw, kx, lx
	
c-----------------------------------------------------------	

	is=0 ; iu=0
	DO kx=k1,k1+ip
	 iw=1 ; iu=iu+1
	 h(iu)=quad(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3),
	1 alr(l1),alr(l2),alr(l3))
	 IF(iq == 3)THEN
	  iw=2
	  q(iu)=quad(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4),
	1 alr(l2),alr(l3),alr(l4))
	 ENDIF
	 is=1
	ENDDO
	is=0 ; iw=1
c.....  k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
	opact=quad(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
	dopact=dkap ; dkap1=dkap
	IF(iq == 3)THEN
c.....   k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
	 opactq=quad(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3))
	 dkapq1=dkap
	ENDIF
	IF(ip == 3)THEN
c.....   k and Dlog(k)/Dlog(T6) in lower-left 3x3.
	 opact2=quad(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4))
	 dkap2=dkap
c.....   k and Dlog(k)/Dlog(T6) smoothed in left 3x4
	 dix=(alt(k3)-slt)*dfs(k3)
	 dopact=dkap1*dix+dkap2*(1.d0-dix)
	 opact=opact*dix+opact2*(1.d0-dix)
	 IF(iq == 3)THEN
c.....    k and Dlog(k)/Dlog(T6) in upper-right 3x3.
	  opactq2=quad(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
	  dkapq2=dkap ; dopactq=dkapq1*dix+dkapq2*(1.d0-dix)
	  opactq=opactq*dix+opactq2*(1.d0-dix)
	 ENDIF
	ENDIF

	iu=0
	DO lx=l1,l1+iq
	 iw=1; iu=iu+1
	 h(iu)=quad(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx),
	1 alt(k1),alt(k2),alt(k3))
	 IF(ip == 3)THEN
	  iw=2
	  q(iu)=quad(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx),
	1 alt(k2),alt(k3),alt(k4))
	 ENDIF
	 is=1
	ENDDO

	is=0 ; iw=1
c.....  k and Dlog(k)/Dlog(R) in lower-left 3x3
	dopacr=dkap
	IF(ip == 3)THEN
c.....   k and Dlog(k)/Dlog(R) in upper-left 3x3.
	 dopacrq=dkap
	ENDIF
	IF(iq == 3)THEN
c.....   k and Dlog(k)/Dlog(R) in lower-right 3x3.
	 opact2=quad(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4))
	 dix2=(alr(l3)-slr)*dfsr(l3)
	 dopacr=dopacr*dix2+dkap*(1.d0-dix2)
c.....   k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
	 dopact=dopact*dix2+dopactq*(1.d0-dix2)
	 opact=opact*dix2+opactq*(1.d0-dix2)
	ENDIF
	IF(ip == 3)THEN
	 IF(iq == 3)THEN
c.....	  k and Dlog(k)/Dlog(R) in upper-right 3x3.
c.....	  Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
	  dopacrq=dopacrq*dix2+dkap*(1.d0-dix2)
	 ENDIF
	 dopacr=dopacr*dix+dopacrq*(1.d0-dix)
	ENDIF
	IF(opact > 1.d+15)THEN
	 WRITE(*,1) ; WRITE(2,1) ; STOP
1	 FORMAT('opa_opalCO, t6rinterp interpolation indices out of range',
	1 'opact=',es10.3)	
	ENDIF
	IF(opact > 9.d0)THEN
	 opact=30.d0 ; dopact=99.d0 ; dopacr=99.d0
	ENDIF
	
	RETURN
      
	END SUBROUTINE t6rinterp

c************************************************************************
#if 0
	SUBROUTINE readco
	
c..... The purpose of this subroutine is to read the data tables

	USE mod_donnees, ONLY : nom_chemin	
	USE mod_kind
	
	IMPLICIT NONE
	
	REAL(kind=dp) :: altin, dum

	INTEGER :: i, int, istep, j, k, kk, l, ll
	
	CHARACTER (len=1), DIMENSION(250) :: dumarra
	CHARACTER (len=80) :: chain

	LOGICAL, SAVE :: isett6=.TRUE.	
	LOGICAL :: ok
  
c------------------------------------------------------------- 

2000	FORMAT(8es10.3)
	
	DO j=1,nc-1
	 DO i=1,nc
	  IF(xcd(j) >= xc(i))THEN
	   n(m,j)=i+1
	   IF(xcd(j) < xc(i)+1.d-6)n(m,j)=i
	  ENDIF
	 ENDDO
	ENDDO
c	PRINT*,nom_chemin
c	PRINT*,m
c	PRINT*,fich_opalCO(m)

	chain=TRIM(nom_chemin)//fich_opalCO(m)
	
c	PRINT*,chain
c	PAUSE'm'

	INQUIRE(file=TRIM(chain),exist=ok)
	IF(.NOT.ok)THEN
	 INQUIRE(file=TRIM(chain)//'.gz',exist=ok)
	 IF(ok)THEN
	  CALL SYSTEM('gunzip '//TRIM(chain)//'.gz')
	 ELSE
	  WRITE(*,5)TRIM(chain) ; WRITE(2,5)TRIM(chain) ; STOP
5	  FORMAT('Arret, Fichier inconnu compressé ou non: ',a)
	 ENDIF
	ENDIF	 
	OPEN(29,file=TRIM(chain),form='formatted',status='old')
	WRITE(*,6)TRIM(chain) ; WRITE(2,6)TRIM(chain)
6	FORMAT('Lecture du fichier d''opacités: ',a)                 
	READ(29,'(a)')dumarra(1:240)
	int=0
	DO j=1,no-1
	 DO i=1,n(m,j)
	  int=int+1 ; READ(29,1)dum
1	  FORMAT(f10.5)	  
	  READ(29,2)itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),
	1 xoa(m,int)
c	  PRINT*,m,j,int,itab(m,int),no,n(m,j) ; PAUSE'itab1'
2	  FORMAT(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)	
	  xca(m,int)=MIN(xca(m,int),1.d0-x(m,int)-zz(m,int)-xoa(m,int))
	  READ(29,1)dum,dum,dum ; READ(29,3)alrf(1:nrm)	; READ(29,1)dum
3	  FORMAT(4x,f6.1,18f7.1)

	  DO k=1,ntm
	   READ(29,4)altin,cof(k,1:nrm)	; coff(k,1:nrm)=cof(k,1:nrm)
4	   FORMAT(f4.2,19f7.3)
	  ENDDO
	  IF(isett6)THEN
	   t6arr(1:ntm)=t6list(1:ntm) ; isett6=.FALSE.
	  ENDIF
	  IF(lisse)THEN
	   IF(nrm /= nr .OR. ntm /= nt)THEN 
	    WRITE (*,'("Not set up to smooth data with reduced T-Rho grid")')
            STOP
	   ENDIF  
	   tmax=10.d0   ! modified
	   nset=67    ! 65 in earlier version
	   rls=-8.d0 ; nsm=1; rle=1.0d0 ; nrlow=1 ; nrhigh=2*(rle-rls)+1
	   CALL opaltab    !modified
	  ENDIF

	  ll=1
	  DO kk=nrb,nre
	   alr(ll)=alrf(kk)
	   DO k=1,nt
	    IF(lisse)THEN
	     IF(m == 1 .AND. k <= 9)THEN
	      co(m,i,j,k,ll)=cof(k+ntb-1,kk)
	     ELSE
	      co(m,i,j,k,ll)=coff(k+ntb-1,kk)
	     ENDIF
	    ELSE
	     co(m,i,j,k,ll)=coff(k+ntb-1,kk)
	    ENDIF
	   ENDDO
	   ll=ll+1
	  ENDDO
	 ENDDO
	ENDDO

	IF(x(m,1) /= xa(m))THEN
	 WRITE(*,'(" x in the codata? file does not match xa(m)")')
	 STOP
	ENDIF

	DO i=1,nc-1
	 DO k=1,nt
	  DO l=1,nr
	   diag(m,i,k,l)=co(m,n(m,i),i,k,l)
	  ENDDO
	 ENDDO
	ENDDO

	DO j=1,no-1
	 int=int+1 ; READ(29,1)dum
         READ(29,2)itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),
	1 xoa(m,int)
c	 PRINT*,m,int,itab(m,int),no ; PAUSE'itab2'	
	 READ(29,1)dum,dum,dum ; READ(29,3)alrf(1:nrm) ; READ(29,1) dum
         DO k=1,ntm
          READ(29,4)dum,cof(k,1:nrm)
c	  set up to smooth final "diago" opacity tables
	  coff(k,1:nrm)=cof(k,1:nrm)
	 ENDDO
	 
c	 smooth final "diago" opacity tables too!
	 IF(lisse)THEN
          tmax=10.d0   ! modified
          nset=67 !65 in earlier version
          rls=-8.d0 ; nsm=1 ; rle=1.0d0 ; nrlow=1 ; nrhigh=2*(rle-rls)+1
	  CALL opaltab    !modified
	  DO k=3,ntemp-2
	   DO ll=nrlow,nrhigh
c	    Following skip required because, due to missing data,
c	    the x=0  low T data cannot be smoothed
	    IF(m == 1 .AND. k <= 9)THEN
	     cof(k,ll)=cof(k,ll)
	    ELSE
	     cof(k,ll)=coff(k,ll)
	    ENDIF
	   ENDDO  
	  ENDDO
	    
	  ll=1
	  DO kk=nrb,nre
	   DO k=1,nt
	    diago(m,j,k,ll)=cof(k+ntb-1,kk)
	   ENDDO
	   ll=ll+1
	  ENDDO
	 ENDIF
	ENDDO
	WRITE(*,7)TRIM(chain) ; WRITE(2,7)TRIM(chain)
7	FORMAT('Fermeture et compression du fichier d''opacités: ',a)                 
    
	CLOSE(unit=29)
	CALL SYSTEM('gzip '//TRIM(chain))	
	 
	DO i=2,nt
	 dfs(i)=1.d0/(alt(i)-alt(i-1))
	ENDDO
	DO i=2,nr
	 dfsr(i)=1.d0/(alr(i)-alr(i-1))
	ENDDO
	istep=-1
	IF(mx > 1 )THEN
	 istep=1
	 DO i=2,mx,istep
	  dfsx(i)=1.d0/(xx(i)-xx(i-1))
	 ENDDO
	ENDIF

	RETURN
      
	END SUBROUTINE readco
#endif

c************************************************************************

	REAL(kind=dp) FUNCTION quad(ic,i,x,y1,y2,y3,x1,x2,x3)
      
c..... this function performs a quadratic interpolation.

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: x, y1, y2, y3, x1, x2, x3
	INTEGER, INTENT(in) :: i, ic		
	
	REAL (kind=dp), DIMENSION(30) :: xx12, xx13, xx23, xx1sq, xx1pxx2  	
	REAL (kind=dp), DIMENSION(3) :: xx, yy
	REAL (kind=dp) :: c1, c2, c3	 

c------------------------------------------------------------------

	xx(1)=x1 ; xx(2)=x2 ; xx(3)=x3 ; yy(1)=y1 ; yy(2)=y2 ; yy(3)=y3
	IF(ic == 0)THEN
	 xx12(i)=1.d0/(xx(1)-xx(2)) ; xx13(i)=1.d0/(xx(1)-xx(3))
	 xx23(i)=1.d0/(xx(2)-xx(3)) ; xx1sq(i)=xx(1)*xx(1)
	 xx1pxx2(i)=xx(1)+xx(2)
	ENDIF
	c3=(yy(1)-yy(2))*xx12(i) ; c3=c3-(yy(2)-yy(3))*xx23(i)
	c3=c3*xx13(i) ; c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3
	c1=yy(1)-xx(1)*c2-xx1sq(i)*c3 ; dkap=c2+(x+x)*c3
	quad=c1+x*(c2+x*c3)
      
	RETURN
      
	END FUNCTION quad

c*********************************************************************

	SUBROUTINE spline(x,y,n,y2)

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: x, y
	INTEGER, INTENT(in) :: n	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: y2	

	REAL (kind=dp), DIMENSION(n) :: u
	REAL (kind=dp) :: p, qn, un, sig, ypn, yp1
	
	INTEGER :: i, k

c------------------------------------------------------

c first derivatives at end points using cubic fit

	yp1=((y(3)-y(1))*(x(2)-x(1))**2
	1 -(y(2)-y(1))*(x(3)-x(1))**2)/
	2 ((x(3)-x(1))*(x(2)-x(1))*(x(2)-x(3)))
        ypn=((y(n-2)-y(n))*(x(n-1)-x(n))**2
	1 -(y(n-1)-y(n))*(x(n-2)-x(n))**2)/
	2 ((x(n-2)-x(n))*(x(n-1)-x(n))*(x(n-1)-x(n-2)))

	y2(1)=-0.5d0
	u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	DO i=2,n-1
	 sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	 p=sig*y2(i-1)+2.d0
	 y2(i)=(sig-1.d0)/p
	 u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
	1 /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
	ENDDO
	qn=0.5d0
	un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
	DO k=n-1,1,-1
	 y2(k)=y2(k)*y2(k+1)+u(k)
	ENDDO

	RETURN
      
	END SUBROUTINE spline
      
c*********************************************************************

	SUBROUTINE splint(xa,ya,n,y2a,x,y,yp)
      
	USE mod_kind
      
	IMPLICIT NONE
      
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xa, ya, y2a	
	REAL (kind=dp), INTENT(in) :: x 
	INTEGER, INTENT(in) :: n
	
	REAL (kind=dp), INTENT(out) :: y, yp  
	REAL (kind=dp) :: a, b, h
	INTEGER :: k, khi, klo       
 
c------------------------------------------------------------------

2000	FORMAT(8es10.3) 
      
	klo=1 ; khi=n
1	IF(khi-klo > 1)THEN
	 k=(khi+klo)/2
	 IF(xa(k) > x)THEN
	  khi=k
	 ELSE
	  klo=k
	 ENDIF
	 GOTO 1
	ENDIF
	h=xa(khi)-xa(klo)
	IF(h == 0.d0)THEN
	 PRINT*,'STOP, dans splint Bad xa input.khi, klo=',khi,klo
	 WRITE(*,2000)h,xa(khi),xa(klo) ; STOP
	ENDIF	 
	a=(xa(khi)-x)/h ; b=(x-xa(klo))/h
	y=a*ya(klo)+b*ya(khi)+
	1 ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
	yp=0.05d0*((-ya(klo)+ya(khi))/h+(-(3.d0*a**2-1)*y2a(klo)
	1 +(3.d0*b**2-1)*y2a(khi))*h/6.d0)
     
	RETURN
      
	END SUBROUTINE splint
      
c*********************************************************************

	SUBROUTINE fity

c this routine makes spline fits for f and fx, and obtains fy and fxy

	USE mod_kind
      
	IMPLICIT NONE
	
	REAL (kind=dp), DIMENSION(ipr) :: a, ad, b, bd
	REAL (kind=dp) :: apn, ap1, bpn, bp1
	INTEGER :: i, j
	
c----------------------------------------------------------------------	

	DO i=1,nset   ! modified
	 a(1:nrl)=f(i,1:nrl) ; b(1:nrl)=fx(i,1:nrl)
	 CALL getd(a,nrl,ad,ap1,apn) ; CALL getd(b,nrl,bd,bp1,bpn)
	 fy(i,1)=ap1 ; fy(i,nrl)=apn ; fxy(i,1)=bp1 ; fxy(i,nrl)=bpn
	 DO j=2,nrl-1
	  fy(i,j)= -a(j)+a(j+1)-2.d0*ad(j)-ad(j+1)
	  fxy(i,j)=-b(j)+b(j+1)-2.d0*bd(j)-bd(j+1)
	 ENDDO
	ENDDO

	RETURN
      
	END SUBROUTINE fity 
      
c*********************************************************************

	SUBROUTINE fitx
      
c  this routine is used only after smoothing.
c  its function is to recompute fx using smoothed f.

	USE mod_kind
      
	IMPLICIT NONE

	REAL (kind=dp), DIMENSION(85) :: a, d
	REAL (kind=dp) :: apn, ap1	
	INTEGER :: i, j
	
c----------------------------------------------------------------------	

	DO j=1,nrl
	 a(1:nset)=f(1:nset,j)
	 CALL getd(a,nset,d,ap1,apn)  ! modified
	 fx(1,j)=ap1 ; fx(nset,j)=apn   ! modified
	 DO i=2,nset-1  ! modified
	  fx(i,j)=-a(i)+a(i+1)-2.d0*d(i)-d(i+1)
	 ENDDO
	ENDDO

	RETURN
      
	END SUBROUTINE fitx

c************************************************************************

	SUBROUTINE getd(f,n,d,fp1,fpn)

c  simplified code for spline coefficients, for case of intervals
c  of unity.

	USE mod_kind
      
	IMPLICIT NONE
      
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: f
	INTEGER, INTENT(in) :: n
	      
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: d     
	REAL (kind=dp), INTENT(out) :: fp1,fpn         
      
	REAL (kind=dp), DIMENSION(85) :: t
	INTEGER :: j
	
c----------------------------------------------------------------------	

	fp1=(-11.d0*f(1)+18.d0*f(2)-9.d0*f(3)+2.d0*f(4))/6.d0
	fpn=(11.d0*f(n)-18.d0*f(n-1)+9.d0*f(n-2)-2.d0*f(n-3))/6.d0
	d(1)=-0.5d0 ; t(1)=0.5d0*(-f(1)+f(2)-fp1)
	DO j=2,n-1
	 d(j)=-1.d0/(4.d0+d(j-1))
	 t(j)=-d(j)*(f(j-1)-2.d0*f(j)+f(j+1)-t(j-1))
	ENDDO
	d(n)=(fpn+f(n-1)-f(n)-t(n-1))/(2.d0+d(n-1))
	DO j=n-1,1,-1
	 d(j)=d(j)*d(j+1)+t(j)
	ENDDO

	RETURN
      
	END SUBROUTINE getd

c*********************************************************************

	SUBROUTINE interp(flt,flrho,g,dgdt,dgdrho,ierr)

c  given f,fx,fy and fxy on the grid points, this routine
c  does bi-cubic interpolations using methods described in
c  numerical recipes, pp. 118 to 120
c       input is flt=LOG10(t), flrho=LOG10(rho)
c       output is g=LOG10(ross)
c              dgdt=dg/d(LOG10(t))
c            dgdrho=dg/d(LOG10(rho))
c              ierr=.true. IF input flt, flrho are out-of-range,
c                          ELSE ierr=.FALSE.

c--------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: flt, flrho	
	REAL (kind=dp), INTENT(out) :: g, dgdt, dgdrho
	LOGICAL, INTENT(out) :: ierr
	
	REAL (kind=dp), DIMENSION(16) :: b	
	REAL (kind=dp) :: ff, ffx, ffy, flr, s, t, u, v, x, y
	INTEGER :: i, j

c--------------------------------------------------------------------

c  extreme limits allowed are:-
c     (3.800-0.0125) to (8.000+0.0125) for LOG10(t)
c     (rls-0.125) to (rle+0.1254) for LOG10(r)
c     (allowing for small extrapolations beyond tabular values)

c  function definitions for cubic expansion

	ff(s,t)=b( 1)+t*(b( 2)+t*(b( 3)+t*b( 4)))
	1 +s*(  b( 5)+t*(b( 6)+t*(b( 7)+t*b( 8)))
	2 +s*(  b( 9)+t*(b(10)+t*(b(11)+t*b(12)))
	3 +s*(  b(13)+t*(b(14)+t*(b(15)+t*b(16))) )))

	ffx(s,t)=   b( 5)+t*(b( 6)+t*(b( 7)+t*b( 8)))
	1 +s*(2.d0*(b( 9)+t*(b(10)+t*(b(11)+t*b(12))))
	2 +s*(3.d0*(b(13)+t*(b(14)+t*(b(15)+t*b(16)))) ))

	ffy(s,t)=   b( 2)+s*(b( 6)+s*(b(10)+s*b(14)))
	1 +t*(2.d0*(b( 3)+s*(b( 7)+s*(b(11)+s*b(15))))
	2 +t*(3.d0*(b( 4)+s*(b( 8)+s*(b(12)+s*b(16)))) ))

	ierr=.FALSE.

	x=20.d0*(flt-3.800d0)+1 ; flr=flrho+18.d0-3.d0*flt
	y=2.d0*( flr - rls )+1.d0

	IF(x < 2.d0)THEN
	 IF(x < 0.75d0)THEN
          ierr=.TRUE.
	 ELSE
	  i=1
	 ENDIF
	ELSEIF(x > 84.d0)THEN
	 IF(x > 85.25d0)THEN
	  ierr=.TRUE.
	 ELSE
	  i=84
	 ENDIF
	ELSE
	 i=x
	ENDIF
	u=x-i

	IF(y < 2.d0)THEN
	 IF(y < 0.75d0)THEN
	  ierr=.TRUE.
	 ELSE
	  j=1
	 ENDIF
	ELSEIF(y > nrl-1)THEN
	 IF(y > nrl+.25d0)THEN
	  ierr=.TRUE.
	 ELSE
	  j=nrl-1
	 ENDIF
	ELSE
	 j=y
	ENDIF
	v=y-j

	IF(ierr)THEN
	 g=9.999d0 ; dgdt=9.999d0 ; dgdrho=9.999d0 ; RETURN
	ENDIF

c  given functions and derivatives at grid points, compute coefficients.

	b(1)=f(i,j) ; b(2)=fy(i,j)
	b(3)=3.d0*(-f(i,j)+f(i,j+1))-2.d0*fy(i,j)-fy(i,j+1)
	b(4)=2.d0*(f(i,j)-f(i,j+1))+fy(i,j)+fy(i,j+1)

	b(5)=fx(i,j) ; b(6)=fxy(i,j)
	b(7)=3.d0*(-fx(i,j)+fx(i,j+1))-2.d0*fxy(i,j)-fxy(i,j+1)
	b(8)=2.d0*(fx(i,j)-fx(i,j+1))+fxy(i,j)+fxy(i,j+1)

	b(9)=3.d0*(-f(i,j)+f(i+1,j))-2.d0*fx(i,j)-fx(i+1,j)
	b(10)=3.d0*(-fy(i,j)+fy(i+1,j))-2.d0*fxy(i,j)-fxy(i+1,j)
	b(11)=9.d0*(f(i,j)-f(i+1,j)+f(i+1,j+1)-f(i,j+1))
	1 +6.d0*(fx(i,j)-fx(i,j+1)+fy(i,j)-fy(i+1,j))+4.d0*fxy(i,j)
	2 +3.d0*(fx(i+1,j)-fx(i+1,j+1)-fy(i+1,j+1)+fy(i,j+1))
	4 +2.d0*(fxy(i,j+1)+fxy(i+1,j))+fxy(i+1,j+1)
	b(12)=6.d0*(-f(i,j)+f(i+1,j)-f(i+1,j+1)+f(i,j+1))
	1 +4.d0*(-fx(i,j)+fx(i,j+1))
	2 +3.d0*(-fy(i,j)+fy(i+1,j)+fy(i+1,j+1)-fy(i,j+1))
	3 +2.d0*(-fx(i+1,j)+fx(i+1,j+1)-fxy(i,j)-fxy(i,j+1))
	4 -fxy(i+1,j)-fxy(i+1,j+1)

	b(13)=2.d0*(f(i,j)-f(i+1,j))+fx(i,j)+fx(i+1,j)
	b(14)=2.d0*(fy(i,j)-fy(i+1,j))+fxy(i,j)+fxy(i+1,j)
	b(15)=6.d0*(-f(i,j)+f(i+1,j)-f(i+1,j+1)+f(i,j+1))
	1 +4.d0*(-fy(i,j)+fy(i+1,j))
	2 +3.d0*(-fx(i,j)-fx(i+1,j)+fx(i+1,j+1)+fx(i,j+1))
	3 +2.d0*(fy(i+1,j+1)-fy(i,j+1)-fxy(i,j)-fxy(i+1,j))
	4 -fxy(i+1,j+1)-fxy(i,j+1)
	b(16)=4.d0*(f(i,j)-f(i+1,j)+f(i+1,j+1)-f(i,j+1))
	1 +2.d0*(fx(i,j)+fx(i+1,j)-fx(i+1,j+1)-fx(i,j+1)
	2 +fy(i,j)-fy(i+1,j)-fy(i+1,j+1)+fy(i,j+1))
	3 +fxy(i,j)+fxy(i+1,j)+fxy(i+1,j+1)+fxy(i,j+1)

c  get g=LOG10(ross), dgdt=d LOG10(ross)/d LOG10(t),
c      dgdrho=d LOG10(ross)/d LOG10(rho)

	g=ff(u,v) ; dgdt=20.d0*ffx(u,v)-6.d0*ffy(u,v)
	dgdrho=2.d0*ffy(u,v)

	RETURN
      
	END SUBROUTINE interp

c*********************************************************************

	SUBROUTINE smooth

c  this subroutine uses a 2-dimensional generalisation of the smoothing
c  techniques described on pp. 644 to 649 of numerical recipes.

c  consider the 25 points defined by
c       i+n, n=-2,-1,0,1,2 and j+m, m=-2,-1,0,1,2.
c  the function to be smoothed is fitted to a bi-cubic, involving
c  16 coefficients, using techniques of least-squares. the smoothed
c  function (temporarily stored in fxy) is given by the fitted value
c  at the point i and j.

c  the fitting is shifted for points close to boundaries.

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp),DIMENSION(6) :: gam=(/+0.0073469388d0,
	1 -0.0293877551d0,-0.0416326531d0,+0.1175510204d0,+0.1665306122d0,
	2 +0.2359183673d0 /)  
	REAL (kind=dp),DIMENSION(11) ::alp=(/
	1 -0.0844897959d0,-0.0048979592d0,+0.0073469388d0,+0.0012244898d0,
	2  0.3379591837d0,+0.0195918367d0,-0.0293877551d0,+0.4787755102d0,
	3  0.0277551020d0,-0.0416326531d0,-0.0069387755d0 /),bet=(/
	4 -0.0048979592d0,-0.0661224490d0,-0.0293877551d0,+0.0195918367d0,
	5  0.2644897959d0,+0.1175510204d0,-0.0783673469d0,+0.0277551020d0,
	6  0.3746938776d0,+0.1665306122d0,-0.1110204082d0 /)

	INTEGER :: i, j

c----------------------------------------------------------------

	DO i=3,nset-2 
	 j=1
	 fxy(i,j)=alp(1)*( f(i-2,j  )+f(i+2,j  ) )
	1 +alp(2)*( f(i-2,j+1)+f(i+2,j+1)+f(i-2,j+3)+f(i+2,j+3)
	2 +f(i-1,j+4)+f(i+1,j+4) )+alp(3)*( f(i-2,j+2)+f(i+2,j+2) )
	3 +alp(4)*( f(i-2,j+4)+f(i+2,j+4))+alp(5)*( f(i-1,j  )+f(i+1,j))
	4 +alp(6)*( f(i-1,j+1)+f(i+1,j+1)+f(i-1,j+3)+f(i+1,j+3) )
	5 +alp(7)*( f(i-1,j+2)+f(i+1,j+2))+alp(8)*  f(i,j  )
	6 +alp(9)*( f(i,j+1)+f(i ,j+3))+alp(10)*f(i,j+2)+alp(11)*f(i,j+4)
	 j=2
	 fxy(i,j)=bet(1)*( f(i-2,j-1)+f(i+2,j-1)+f(i-2,j+3)+f(i+2,j+3) )
	1 +bet(2)*( f(i-2,j  )+f(i+2,j  ) )+bet(3)*( f(i-2,j+1)+f(i+2,j+1))
	2 +bet(4)*( f(i-2,j+2)+f(i+2,j+2)+f(i-1,j-1)+f(i+1,j-1)
	3 +f(i-1,j+3)+f(i+1,j+3) )+bet(5)*( f(i-1,j  )+f(i+1,j  ) )
	4 +bet(6)*( f(i-1,j+1)+f(i+1,j+1) )+bet(7)*( f(i-1,j+2)+f(i+1,j+2))
	5 +bet(8)*( f(i  ,j-1)+f(i  ,j+3) )
	6 +bet(9)*f(i  ,j  ) +bet(10)*f(i  ,j+1) +bet(11)*f(i  ,j+2)
	 DO j=3,nrl-2
	  fxy(i,j)=gam(1)*( f(i-2,j-2)+f(i-2,j+2)+f(i+2,j-2)+f(i+2,j+2) )
	1 +gam(2)*( f(i-2,j+1)+f(i-2,j-1)+f(i-1,j-2)+f(i-1,j+2)
	2 +f(i+1,j-2)+f(i+1,j+2)+f(i+2,j-1)+f(i+2,j+1) )
	3 +gam(3)*( f(i-2,j  )+f(i+2,j  )+f(i  ,j-2)+f(i  ,j+2) )
	4 +gam(4)*( f(i-1,j-1)+f(i-1,j+1)+f(i+1,j-1)+f(i+1,j+1) )
	5 +gam(5)*( f(i-1,j  )+f(i  ,j-1)+f(i  ,j+1)+f(i+1,j  ) )
	6 +gam(6)*  f(i  ,j  )
	 ENDDO
	 j=nrl-1
	 fxy(i,j)= bet(1)*( f(i-2,j+1)+f(i+2,j+1)+f(i-2,j-3)+f(i+2,j-3) )
	1 +bet(2)*( f(i-2,j  )+f(i+2,j  ) )+bet(3)*( f(i-2,j-1)+f(i+2,j-1) )
	2 +bet(4)*( f(i-2,j-2)+f(i+2,j-2)+f(i-1,j+1)+f(i+1,j+1)
	3 +f(i-1,j-3)+f(i+1,j-3) )+bet(5)*( f(i-1,j  )+f(i+1,j  ) )
	4 +bet(6)*( f(i-1,j-1)+f(i+1,j-1) )+bet(7)*( f(i-1,j-2)+f(i+1,j-2) )
	5 +bet(8)*( f(i  ,j+1)+f(i  ,j-3) )
	6 +bet(9)*f(i  ,j  ) +bet(10)*f(i  ,j-1) +bet(11)*f(i  ,j-2)
	 j=nrl
	 fxy(i,j)=alp(1)*( f(i-2,j  )+f(i+2,j  ) )
	1 +alp(2)*( f(i-2,j-1)+f(i+2,j-1)+f(i-2,j-3)+f(i+2,j-3)
	2 +f(i-1,j-4)+f(i+1,j-4) )+alp(3)*( f(i-2,j-2)+f(i+2,j-2) )
	3 +alp(4)*( f(i-2,j-4)+f(i+2,j-4) )+alp(5)*( f(i-1,j  )+f(i+1,j  ) )
	4 +alp(6)*( f(i-1,j-1)+f(i+1,j-1)+f(i-1,j-3)+f(i+1,j-3) )
	5 +alp(7)*( f(i-1,j-2)+f(i+1,j-2) )+alp(8)* f(i,j  )
	6 +alp(9)*( f(i  ,j-1)+f(i  ,j-3) )+alp(10)*f(i,j-2)+alp(11)*f(i,j-4)
	ENDDO
	DO i=3,nset-2   ! modified
	 DO j=1,nrl
	  f(i,j)=fxy(i,j)
	 ENDDO
	ENDDO

	RETURN
      
	END SUBROUTINE smooth

c*********************************************************************

	SUBROUTINE opaltab

c  code for fitting and smoothing opal data. adapted from a code
c     written by mike seaton(obtained june 1993)
c
c     opal data.
c     assumes first t6=0.006, last t6=10.or 0.04). depending on position
c     in the table. 
c     uses rectangular array for variables t6 and LOG10(r)
c
c     (1) nsm=number of passes through smoothing filter.
c     use of nsm=1 or 2 is recommended.
c     no smoothing with nsm=0
c     (2) range for LOG10(r),
c     rls=first value, rle=last vale
c     (rls must be first valuye in table)

c  subroutine interp
c     after processing, data are in a form for use of
c               subroutine interp
c     which gives log(ross) and two first derivatives for any
c     values of log(t) and log(rho). see below for further
c     explanation.

c  output for the case of nsm>0.
c     interp is used to obtain smoothed data interpolated
c     back to the original opal mesh. two files are written.


c  the subroutines spline and splint are adapted from those give by
c  w.h. press, s.a. teulolsky, w.t. vettering and b.p. flannery,
c  "numerical recipes in fortran", 2nd edn., 1992, c.u.p.
c  other references are made to methods described in that book.

c--------------------------------------------------------------
 
	USE mod_kind

	IMPLICIT NONE
	
	INTEGER, PARAMETER :: ip=100
      
	REAL(kind=dp), DIMENSION(ip,ipr) :: rossl     
	REAL(kind=dp), DIMENSION(ip) :: u,v, v2
	REAL(kind=dp) :: dgdrho, dgdt, flrho, flr, flt, g, t6	

	INTEGER :: i, j, k, l, ns
	LOGICAL :: ierr
      
c-------------------------------------------------------------- 

	nrl=2*(rle-rls)+1
c	store LOG10(t) in u and LOG10(ross) in rossl
c	check first value of t6
	t6=t6arr(1)
	DO j=1,nrl
	 rossl(1,j)=coff(1,j)
	ENDDO

	IF(ABS(t6-.0056341325d0) < 1.d-8)THEN
	 u(1)=6.d0+LOG10(t6)
	ENDIF
	
C	set rossl up to t6=t6arr(nset)

	i=1
5	i=i+1 ; t6=t6arr(i)
	DO j=1,nrl
	 rossl(i,j)=coff(i,j)
	ENDDO
	u(i)=6.d0+LOG10(t6)
	IF(t6 < tmax)GOTO 5
	 ntemp=i
	IF(ntemp > ip)THEN
	 PRINT*,'REQUIRE PARAMETER ip OF AT LEAST ',ntemp ; STOP
	ENDIF

c     define variables
c         x=20.0*(LOG10(t)-3.80)+1
c         y=2.0*(LOG10(r)-rls)+1
c     use indices i=1 to nset and j=1 to nrl
c     x and y are such that, on mesh-point (i,j), x=i and y=j
c     obtain:-
c         f(i,j)=LOG10(ross)
c         fx(i,j)=df/dx
c         fy(i,j)=df/dy
c         fxy(i,j)=ddf/dxdy

c     first get f and fx, interpolating from opal t6 to
c     interval of 0.05 in LOG10(t).
	DO j=1,nrl
c       for each LOG10(r), store LOG10(ross) in v(i)
	 DO i=1,ntemp
	  v(i)=rossl(i,j)
	 ENDDO

c        get first derivatives at end points
c        get second derivatives for spline fit
	 CALL spline(u,v,ntemp,v2)

c        interpolate to LOG10(t)=flt, flt=3.8(0.05)8.0
	 DO i=1,nset ! modified
	  flt=3.75d0+0.05d0*i
	  CALL splint(u,v,ntemp,v2,flt,f(i,j),fx(i,j))
	 ENDDO
	ENDDO

C  option for smoothing
	IF(nsm > 0)THEN
	 DO ns=1,nsm
	  CALL smooth
	 ENDDO
	 CALL fitx
	ENDIF

C  get fy and fxy
	CALL fity

c  the arrays f, fx, fy and fxy are now stored
c  can now DO interpolations using
c       CALL interp(flt,flrho,g,dgdt,dgdrho,ierr)
c       input is flt=LOG10(t), flrho=LOG10(rho)
c       output is g=LOG10(ross)
c              dgdt=dg/d(LOG10(t))
c            dgdrho=dg/d(LOG10(rho))
c              ierr=.true. IF input flt, flrho are out-of-range,
c                          ELSE ierr=.FALSE.

c interpolate back to opal points
	IF(nsm > 0)THEN
	 DO l=1,nrl
	  coff(1,l)=rossl(1,l)
	 ENDDO
	 DO k=2,ntemp
	  flt=u(k)
	  DO l=nrlow,nrhigh
	   flr=rls+.5*(l-1)
	   flrho=flr-18.d0+3.d0*flt
	   CALL interp(flt,flrho,g,dgdt,dgdrho,ierr)
	   IF(ierr)THEN
	    STOP 'STOP -- opaltab: interp t/rho range error.'	   
	   ENDIF
	   v(l)=g
	  ENDDO
	  t6=t6arr(k)
	  DO l=nrlow,nrhigh
	   coff(k,l)=v(l)
	  ENDDO
	 ENDDO
	ENDIF
      
	RETURN

1000	FORMAT('SMOOTHED OPAL DATA')
1100	FORMAT('OPAL DATA, (SMOOTHED-ORIGINAL)')
2000	FORMAT(100a1)
2222	FORMAT(f8.3,20f7.3)
6000	FORMAT(/'FIRST T6=',es10.3,', SHOULD BE 0.006')
6003	FORMAT(/'!!! OUT-OF-RANGE !!!'/' FLT=',es10.3,', FLRHO=',es10.3,
	1 ', FLR=',es10.3)

	END SUBROUTINE opaltab

	END SUBROUTINE opa_opalCO
