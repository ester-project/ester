
c************************************************************************

	SUBROUTINE etat_opal5Z(pp,tt,xchim,deriv,
	1    ro,drop,drot,drox,u,dup,dut,dux,
	2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c	routine private du module mod_etat

c	interface de l'équation d'état OPAL 2005 avec CESAM2k
c	cas où on a effectué une tabulation pour Z fixé
c	autre cas: etat_opalX avec Z=0 et X dans [0,1]	
c	package OPAL_EOS
c	appel à EFF en cas de difficultés 

c	CESAM2k 

c       entree : 
c	p : pression 
c	t : temperature 
c	xchim : composition chimique 
c	deriv=.TRUE. : evite le calcul de certaines dérivées 

c       sortie : 
c	ro : densité et dérivées 
c	u : énergie interne et dérivées 


c------------------------------------------------------------------------

	USE mod_donnees, only : aradia, f_eos, nchim, nom_chemin, z0 
	USE mod_kind

	IMPLICIT NONE

 	REAL (kind=dp), PARAMETER :: aprop=83.14510d0, sigmacc=1.8914785d-3 

c      	iorder=9 gives all 1st and 2nd order data. See instructions in esac. 
c     	irad=0 does not add radiation; irad=1 adds radiation 

	INTEGER, PARAMETER :: iorder=9, irad=1, mx=5, mv=10, nr=169, 
	1    nt=197

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: pp, tt
	LOGICAL, INTENT(in) :: deriv

	REAL (kind=dp), INTENT(out) :: ro,drop,drot,drox,u,dup,dut,dux, 
	1    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx, 
	2    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1 

	REAL (kind=dp), SAVE, DIMENSION(mx,mv,nt,nr) ::  xz	 
	REAL (kind=dp), SAVE, DIMENSION(mx,nt,nr) ::  epl 
	REAL (kind=dp), SAVE, DIMENSION(mx,nr) :: rhogr 
	REAL (kind=dp), SAVE, DIMENSION(mx,6) :: frac		 
	REAL (kind=dp), SAVE, DIMENSION(nr,nt) :: amu_M,alogN_e, t6list 
	REAL (kind=dp), SAVE, DIMENSION(nt,nr) :: esk, esk2 
	REAL (kind=dp), SAVE, DIMENSION(mv) :: eos	 
	REAL (kind=dp), PARAMETER, DIMENSION(mx) :: xa=(/0.d0, 0.2d0, 
	1    0.4d0, 0.6d0, 0.8d0 /) 
	REAL (kind=dp), SAVE, DIMENSION(mx) :: dfsx, moles, tmass, xin, 
	1    xx, zz 
	REAL (kind=dp), SAVE, DIMENSION(nr) :: dfsr, rho 
	REAL (kind=dp), SAVE, DIMENSION(nt) :: dfs, t6a	 
	REAL (kind=dp), DIMENSION(nchim) :: Lxchim		 
	REAL (kind=dp), SAVE, DIMENSION(4) :: q, h	 
	REAL (kind=dp), PARAMETER :: dx=1.d-3, unpdx=1.d0+dx	 
	REAL (kind=dp), SAVE :: aradias3, esact 
	REAL (kind=dp) :: r, stor, stor0, dstor, x, ztab, t6, p12, p, t, 
	1    duro	 
	
	INTEGER, SAVE, DIMENSION(mx,nr) :: icycuse 
      	INTEGER :: i 
      	INTEGER, PARAMETER, DIMENSION(nt) :: nra=(/(169,i=1,26),168,167,
	1    167,166,165,164,163,163,162,161,160,159,159,158,130,130,129,128,
	2    128,126,126,125,125,124,122,121,120,120,119,119,(118,i=1,6),117,
	3    117,116,116,115,115,(114,i=1,4),(113,i=1,8),(111,i=1,22),
	4    (110,i=1,5),(109,i=1,6),108,108,(107,i=1,5),106,106,105,104,104,
	5    (103,i=1,8),(102,i=1,16),(100,i=1,21),(99,i=1,9),(98,i=1,6),
	6    (97,i=1,4),95,94,(87,i=1,6) /) 
	
      	INTEGER, PARAMETER, DIMENSION(10) :: index=(/ 1, 2, 3, 4, 5, 6, 7, 
	1    8, 9, 10 /) 
      	INTEGER, SAVE, DIMENSION(10) :: iri	 
	INTEGER, SAVE :: itime, l1,l2,l3,l4,k1,k2,k3,k4,ip,iq, m, mf  
	
	LOGICAL, SAVE :: init=.TRUE.	 
	LOGICAL :: ok 
	
c-----------------------------------------------------------------
	
	p=pp ; t=tt	 
      	IF(init)THEN 
	   init=.FALSE. ; aradias3=aradia/3.d0 
	   
	   WRITE(2,1) ; WRITE(*,1) 
 1	   FORMAT(/,'-------------Equation d''état-----------------',/, 
	1	'Equation d état OPAL5Z si T dans [ 1870 , 2.e8 ], EFF sinon',/, 
	2	'der. num.; (p,T)-->(ro,T) calculé par rhoofp',/)  
	ENDIF 
	
	t6=t*1.d-6 ; p12=p*1.d-12 ; x=xchim(1) ; ztab=z0 
	
      	ro=rhoofp(x,ztab,t6,p12,irad,ok) ! calculate density (ro) for fixed P. 
      	IF(ok)THEN      	 
	   CALL esac(x,ztab,t6,ro,iorder,irad,ok) ! calc EOS; use ro from rhoofp       	 
	   IF(ok)THEN 
	      drop=ro/p/eos(6)	!1 / dp/dro 
	      drot=-ro/t*eos(7)/eos(6) !- dp/dT / dp/dro	 
	      u=eos(2)*1.d12 ; duro=-p*(eos(7)-1.d0)/(ro**2) 
	      dup=duro*drop	!du/dT dro/dp 
	      dut=duro*drot+eos(5)*1.e6	!du/dro  dro/dT + du/dt	 
	      delta=eos(7)/eos(6) !Ki T/Ki ro 
	      gradad=1.d0/eos(9) ; cp=p/ro/t*eos(7)*(eos(9)/eos(8)+delta)	 
	      alfa=1.d0/eos(6)	!1/Ki ro 
	      beta=1.d0-aradias3*t**4/p ; gamma1=eos(8)	 
	      IF(deriv)THEN 
		 
c       dérivées / P 
		 
		 stor0=p ; stor=stor0*unpdx  
		 IF(stor < dx)stor=dx 
		 dstor=stor-stor0 ; p=stor ; p12=p*1.d-12 
		 r=rhoofp (x,ztab,t6,p12,irad,ok) ! calculate density (r) for fixed P. 
		 IF(ok)THEN 
		    CALL esac(x,ztab,t6,r,iorder,irad,ok) ! calc EOS; use r from rhoofp 
		    IF(ok)THEN 
		       drop=(r-ro)/dstor ; dup=(eos(2)*1.d12-u)/dstor 
		       deltap=(eos(7)/eos(6)-delta)/dstor 
		       dgradadp=(1.d0/eos(9)-gradad)/dstor 
		       dcpp=(p/r/t*eos(7)*(eos(9)/eos(8)+eos(7)/eos(6))-cp)/dstor 
		       p=stor0 ; p12=p*1.d-12 
		       
c       dérivées / T	 
		       
		       stor0=t ; stor=stor0*unpdx  
		       IF(stor < dx)stor=dx 
		       dstor=stor-stor0 ; t=stor ; t6=t*1.d-6	 
		       r=rhoofp(x,ztab,t6,p12,irad,ok) ! calculate density (r) for fixed P. 
		       IF(ok)THEN    	 
			  CALL esac(x,ztab,t6,r,iorder,irad,ok) ! calc EOS; use r from rhoofp 
			  IF(ok)THEN 
			     drot=(r-ro)/dstor ; dut=(eos(2)*1.d12-u)/dstor 
			     deltat=(eos(7)/eos(6)-delta)/dstor 
			     dgradadt=(1./eos(9)-gradad)/dstor 
			     dcpt=(p/r/t*eos(7)*(eos(9)/eos(8)+eos(7)/eos(6))-cp)/dstor 
			     t=stor0 ; t6=t*1.d-6         
			     
c       dérivées / X	 
			     
			     stor0=xchim(1) ; stor=stor0*unpdx 
			     IF(stor < dx)stor=dx 
			     dstor=stor-stor0 ; x=stor 
			     r=rhoofp (x,ztab,t6,p12,irad,ok) !calculate density (r) for fixed P. 
			     IF(ok)then 	 
				CALL esac(x,ztab,t6,r,iorder,irad,ok) !calc EOS use r from rhoofp 
				IF(ok)THEN 
				   drox=(r-ro)/dstor ; dux=(eos(2)*1.d12-u)/dstor 
				   deltax=(eos(7)/eos(6)-delta)/dstor 
				   dgradadx=(1./eos(9)-gradad)/dstor 
				   dcpx=(p/r/t*eos(7)*(eos(9)/eos(8)+eos(7)/eos(6))-cp)/dstor 
				   x=stor0 
				   RETURN !C bon avec deriv 
				ENDIF !dérivées / X 
			     ENDIF !dérivées / X 
			  ENDIF	!dérivées / T 
		       ENDIF	!dérivées / T 
		    ENDIF	!dérivées / P 
		 ENDIF		!dérivées / P 
	      ENDIF		!deriv  
	      RETURN		!C bon sans deriv   
	   ENDIF		!esac 
	ENDIF			!rhoofp 
	
c	en cas de Pb avec opal appel à EFF 
	
	Lxchim(1:nchim)=xchim(1:nchim) 
	CALL etat_eff(p,t,Lxchim, 
	1    ro,drop,drot,drox,u,dup,dut,dux, 
	2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx, 
	3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1) 
	
	RETURN 
	
	CONTAINS 
	
c*************************************************************************
	
	SUBROUTINE esac (xh,ztab,t6,r,iorder,irad,ok) 
	
c..... 	The purpose of this subroutine is to interpolate  
c      	the equation of state and its derivatives in X, T6, density 
c       izi=0 recalulate table indices to use; =1 keep previous 
	
c       xh=hydrogen mass fraction 
c       ztab is metal fraction of the EOSdata tables you are using. 
c       included only for purpose of preventing mismatch 
c       t6=T6=temperature in millions of degrees kelvin 
c       r=rho=Rho=density(g/cm**3) 
c..... 	to use esac in the calling routine. 
c      	This COMMON contains the interpolated EOS values for the EOS 
c       
c..... 	eos(i) are obtained from a quadradic interpolation at 
c      	fixed T6 at three values of Rho; followed by quadratic 
c      	interpolation along T6. Results smoothed by mixing 
c      	overlapping quadratics. 
c       definitions: 
c       
c       T6 is the temperature in units of 10**6 K 
c       
c       rho is the density in grams/cc 
c       R=Rho/T6**3 
	
c       eos(1) is the pressure in megabars (10**12dyne/cm**2)
c       eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
c       eos(3) is the entropy in units of energy/T6
c       eos(4) is dE/dRHO at constant T6
c       eos(5) is the specific heat, dE/dT6 at constant V.
c       eos(6) is dlogP/dlogRho at constant T6. 
c       Cox and Guil1 eq 9.82
c       eos(7) is dlogP/dlogT6 at conxtant Rho.
c       Cox and Guil1 eq 9.81
c       eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
c       eos(9) is gamma2/(gaamma2-1). Eqs. 9.88 Cox and Guili
	
c       iorder sets maximum index for eos(i);i.e., iorder=1 
c       gives just the pressure 
c       
c       irad  if =0 no radiation correction; if =1 adds radiation 
	
c       index(i),i=1,10  sets order in which the equation of state 
c       variables are stored in eos(i).  Above order corresponds 
c       to block data statement: 
c       data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/.  
c       If you, for example, only want to RETURN gamma1: set iorder=1  
c       and set: data (index(i),i=1,10)/8,2,3,4,5,6,7,1,9,10/ 
	
c--------------------------------------------------------------
	
	USE mod_kind 
	
	IMPLICIT NONE 
	
	LOGICAL, INTENT(out) :: ok 
	
	REAL (kind=dp), DIMENSION(7) :: frac	 
	REAL (kind=dp) :: slr,slt,moles,xh,ztab,t6,r,dixr,p0,z,sum1,sum2, 
	1    sum23,sum33,tmass,eground,fracz 
	
	INTEGER	:: iorder,irad,i,j,ilo,ihi,imd,ipu,iqu,mg,mh,mi,mf2,ms, 
	1    ir,it,is,iw,iv 
	
	LOGICAL :: init=.TRUE., logic1=.TRUE., logic2=.TRUE. 
	
	SAVE 
	
c--------------------------------------------------------------------
	
	IF(init)THEN 
	   init=.FALSE.	  
	   IF((irad /= 0) .AND. (irad /= 1))THEN 
	      WRITE (*,'(" Irad must be 0 or 1")') ; STOP 
	   ENDIF 
	ENDIF 
	
	slt=t6 ; slr=r 
	IF(itime /= 12345678)THEN 
	   itime=12345678 
	   DO i=1,10 
	      DO j=1,10 
		 IF(index(i) == j)iri(i)=j 
	      ENDDO 
	   ENDDO 
	   DO i=1,mx 
	      xx(i)=xa(i) 
	   ENDDO 
	   
c.....  read the data files 
	   
	   CALL r_opal_bin 
	   z=zz(1) 
	   IF(ABS(ztab-z) > 0.03d0)THEN 
	      WRITE(*,10)ztab,z ; WRITE(2,10)ztab,z 
 10	      FORMAT(/,'!!ATTENTION, Z=',es10.3, 
	1	   ' est très différent du Z de la table d''équation d''état, Z=', 
	2	   es10.3,/,'!!ATTENTION',/) 
	   ENDIF 
	   IF(z+xh-1.d-6 > 1.d0 )GOTO 61 
	ENDIF 
	
	ok=.TRUE. 
	
c..... 	Determine T6,rho grid points to use in the interpolation. 
	IF((slt > t6a(1)) .OR. (slt < t6a(nt)) .OR. 
	1    (slr < rho(1)) .OR. (slr > rho(nr)))THEN 
c       PAUSE'logic1' 
	   IF(logic1)THEN 
	      WRITE(*,21)slt,slr ; WRITE(2,21)slt,slr ;logic1=.FALSE. 
 21	      FORMAT(/,'!!!!!ATTENTION!!!!!!!!!!',/, 
	1	   'on a recontré au moins 1 fois le Pb. suivant dans la',/, 
	2	   'routine "esac" de livermore, on utilise MHD  ou EFF',/,  
	3	   'T6 ou LogR outside of table range, T6=',es10.3, 
	4	   ' LogR=',es10.3,/,'!!!!!ATTENTION!!!!!!!!!!',/)		 
	   ENDIF 
	   ok=.FALSE. ; RETURN 
	ENDIF 
	
c.......IF(izi == 0)THEN  ! freeze table indices if not 0 
	ilo=2 ; ihi=mx 
 8	IF(ihi-ilo > 1)THEN 
	   imd=(ihi+ilo)/2 
	   IF(xh <= xa(imd)+1.d-7)THEN 
	      ihi=imd 
	   ELSE 
	      ilo=imd 
	   ENDIF 
	   GOTO 8 
	ENDIF 
	i=ihi ; mf=i-2 ; mg=i-1 ; mh=i ; mi=i+1 ; mf2=mi 
	IF(xh < 1.d-6)THEN 
	   mf=1 ; mg=1 ; mh=1 ; mi=2 ; mf2=1 
	ENDIF 
	IF((xh <= xa(2)+1.d-7) .OR. (xh >= xa(mx-2)-1.d-7)) mf2=mh 
	ilo=2 ; ihi=nr 
 12	IF(ihi-ilo > 1)THEN 
	   imd=(ihi+ilo)/2 
	   IF(slr == rho(imd))THEN 
	      ihi=imd 
	      GOTO 13 
	   ENDIF 
	   IF(slr <= rho(imd))THEN 
	      ihi=imd 
	   ELSE 
	      ilo=imd 
	   ENDIF 
	   GOTO 12 
	ENDIF 
 13	i=ihi ; l1=i-2 ; l2=i-1 ; l3=i ; l4=l3+1 ; iqu=3 
	IF(l4 > nr) iqu=2 
	ilo=nt ; ihi=2 
 11	IF(ilo-ihi > 1)THEN 
	   imd=(ihi+ilo)/2 
	   IF(t6 == t6list(1,imd))THEN 
	      ilo=imd ; GOTO 14 
	   ENDIF 
	   IF(t6 <= t6list(1,imd))THEN 
	      ihi=imd 
	   ELSE 
	      ilo=imd 
	   ENDIF 
	   GOTO 11 
	ENDIF 
 14	i=ilo ; k1=i-2 ; k2=i-1 ; k3=i ; k4=k3+1 ; ipu=3 
	IF(k4 > nt) ipu=2 
	IF(k3 == 0)THEN 
	   WRITE (*,'(" ihi,ilo,imd",3i5)') 
	ENDIF 
	
c.......ENDIF 
	
c     	check to determine if interpolation indices fall within 
c     	table boundaries.  choose largest allowed size. 
	
	sum1=0.d0 ; sum2=0.d0 ; sum23=0.d0 ; sum33=0.d0 
	DO m=mf,mf2 
	   DO ir=l1,l1+1 
	      DO it=k1,k1+1 
		 sum1=sum1+xz(m,1,it,ir) 
	      ENDDO 
	   ENDDO 
	   DO ir=l1,l1+2 
	      DO it=k1,k1+2 
		 sum2=sum2+xz(m,1,it,ir) 
	      ENDDO 
	   ENDDO 
	   IF(ipu == 3)THEN 
	      DO ir=l1,l1+2 
		 DO it=k1,k1+ipu 
		    sum23=sum23+xz(m,1,it,ir) 
		 ENDDO 
	      ENDDO 
	   ELSE 
	      sum23=2.d+30 
	   ENDIF 
	   IF(iqu == 3)THEN 
	      DO ir=l1,l1+3 
		 DO it=k1,k1+ipu 
		    sum33=sum33+xz(m,1,it,ir) 
		 ENDDO 
	      ENDDO 
	   ELSE 
	      sum33=2.d+30 
	   ENDIF 
        ENDDO 
        iq=2 ; ip=2 
        IF(sum2 > 1.d+30)THEN 
	   IF(sum1 < 1.d+25 )THEN 
	      k1=k3-3 ; k2=k1+1; k3=k2+1 ; l1=l3-3 ; l2=l1+1 ; l3=l2+1 
	      GOTO 15 
	   ELSE 
c       PAUSE'logic2' 
	      IF(logic2)THEN 
		 WRITE(*,22)xh,t6*1.d6,r ; WRITE(2,22)xh,t6*1.d6,r	  
 22		 FORMAT(/,'!!!!!ATTENTION!!!!!!!!!!',/, 
	1	      'on a recontre au moins 1 fois le Pb. suivant dans la',/, 
	2	      'routine "esac" de livermore, on utilise EFF',/,  
	3	      'T6/log rho in empty region of table (65)',/, 
	4	      'X=',es10.3,' T=',es10.3,' ro=',es10.3,/, 
	5	      '!!!!!ATTENTION!!!!!!!!!!',/) 
		 logic2=.FALSE. 
	      ENDIF 
	      ok=.FALSE. ; RETURN 
	   ENDIF 
	ENDIF 
	IF(sum23 < 1.d+30) ip=3 ; IF(sum33 < 1.d+30)iq=3 
	IF(t6 >= t6list(1,2)+1.d-7)ip=2 
	IF(slr <= rho(2)+1.d-7)iq=2 
	IF((l3 == nr) .OR. (k3 == nt))THEN 
	   iq=2 ; ip=2 
	ENDIF 
	
 15	DO iv=1,iorder 
	   DO m=mf,mf2 
	      is=0 
	      DO ir=l1,l1+iq 
		 DO it=k1,k1+ip 
		    epl(m,it,ir)=xz(m,iv,it,ir) ; is=1 
		 ENDDO 
	      ENDDO 
	   ENDDO 
	   IF(zz(mg) /= zz(mf) .OR. zz(mh) /= zz(mf))THEN 
	      WRITE(*,'("Z does not match Z in EOSdata files you are using")') 
	      STOP 
	   ENDIF 
	   IF(z/= zz(mf))GOTO 66 
	   is=0 ; iw=1 
	   DO ir=l1,l1+iq 
	      B1: DO it=k1,k1+ip 
	      IF(mf2 == 1)THEN 
		 esk(it,ir)=epl(mf,it,ir) ; CYCLE B1 
	      ENDIF 
	      esk(it,ir)=quadeos(is,iw,xh,epl(mf,it,ir),epl(mg,it,ir), 
	1	   epl(mh,it,ir),xx(mf),xx(mg),xx(mh)) 
	      IF(esk(it,ir) > 1.d+20)THEN 
		 WRITE(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir,l3,k3,iq,ip 
		 WRITE(*,'(3es12.4)') (epl(ms,it,ir),ms=mf,mf+2) 
	      ENDIF 
	      is=1 
	   ENDDO B1 
	ENDDO 
	
	IF(mi == mf2)THEN	! interpolate between quadratics 
	   is=0 ; iw=1 
	   dixr=(xx(mh)-xh)*dfsx(mh) 
	   DO ir=l1,l1+iq 
	      DO it=k1,k1+ip 
		 esk2(it,ir)=quadeos(is,iw,xh,epl(mg,it,ir),epl(mh,it,ir), 
	1	      epl(mi,it,ir),xx(mg),xx(mh),xx(mi)) 
		 IF(esk(it,ir) > 1.d+20)THEN 
		    WRITE(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir,l3,k3,iq,ip 
		    WRITE(*,'(3es12.4)')  (epl(ms,it,ir),ms=mg,mg+2) 
		 ENDIF 
		 esk(it,ir)=esk(it,ir)*dixr+esk2(it,ir)*(1.d0-dixr) ; is=1 
	      ENDDO 
	   ENDDO 
	ENDIF 
	
	is=0 
	
c.....  completed X interpolation. Now interpolate T6 and rho on a 
c       4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure 
c       mixes overlapping quadratics to obtain smoothed derivatives. 
	
	CALL t6rinteos(slr,slt) ; eos(iv)=esact 
	ENDDO 
	
	p0=t6*r ; eos(iri(1))=eos(iri(1))*p0 ! interpolated in p/po 
	eos(iri(2))=eos(iri(2))*t6 ! interpolated in E/T6 
	tmass=gmass(xh,z,moles,eground,fracz,frac) 
	IF(irad == 1)THEN 
	   CALL radsub (t6,r,moles,tmass) 
	ELSE 
	   eos(iri(5))=eos(iri(5))*moles*aprop/tmass ! cv
	ENDIF 
	
	RETURN 
	
 61	WRITE(*,'(" Mass fractions exceed unity (61)")') ; STOP 
 62	WRITE(*,'(" T6/LogR outside of table range (62)")') ; STOP 
 65	WRITE(*,'("T6/log rho in empty region of table (65)")') 
	WRITE (*,'("xh,t6,r=", 3e12.4)') xh,t6,r ; STOP 
 66	WRITE(*,'("Z does not match Z in EOSdata* files you are using (66)")') 
	WRITE (*,'("mf,zz(mf)=",i5,e12.4)') mf,zz(mf) 
	WRITE (*,'("  iq,ip,k3,l3,xh,t6,r,z= ",4i5,4e12.4)')ip,iq,k3,l3, 
	1    xh,t6,r,z 
	STOP 
	
	END SUBROUTINE esac 
	
c*********************************************************************** 
 
	SUBROUTINE t6rinteos(slr,slt) 
	 
c     	The purpose of this subroutine is to interpolate in T6 and rho 
 
	USE mod_kind 
		 
	IMPLICIT NONE 
 
c	INTEGER, PARAMETER :: mx=5,mv=10,nr=169,nt=191 
	 
	INTEGER	:: iu,is,kx,iw 
	 
	REAL (kind=dp) :: slr,slt,dix,dix2,esact2,esactq,esactq2 
      
c---------------------------------------------------------------------	 
 
	iu=0 ; is=0 
 
	DO kx=k1,k1+ip 
	 iw=1 ; iu=iu+1 
	 h(iu)=quadeos(is,iw,slr,esk(kx,l1),esk(kx,l2),esk(kx,l3),rho(l1), 
	1 rho(l2),rho(l3)) 
	 IF(iq == 3)THEN 
	  iw=2 
	  q(iu)=quadeos(is,iw,slr,esk(kx,l2),esk(kx,l3),esk(kx,l4),rho(l2), 
	1 rho(l3),rho(l4)) 
	 ENDIF 
	 is=1 
	ENDDO 
	is=0 ; iw=1 
c..... 	eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2) 
	esact=quadeos(is,iw,slt,h(1),h(2),h(3),t6a(k1),t6a(k2),t6a(k3)) 
	IF(iq == 3)THEN 
c.....   eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2) 
	 esactq=quadeos(is,iw,slt,q(1),q(2),q(3),t6a(k1),t6a(k2),t6a(k3)) 
	ENDIF 
	IF(ip == 3)THEN 
c.....   eos(i) in lower-left 3x3. 
	 esact2=quadeos(is,iw,slt,h(2),h(3),h(4),t6a(k2),t6a(k3),t6a(k4)) 
c.....   eos(i) smoothed in left 3x4 
	 dix=(t6a(k3)-slt)*dfs(k3) 
	 esact=esact*dix+esact2*(1.-dix) 
c	 ENDIF   ! moved to loc a 
	 IF(iq == 3)THEN 
  
c.....    eos(i) in upper-right 3x3. 
	  esactq2=quadeos(is,iw,slt,q(2),q(3),q(4),t6a(k2),t6a(k3),t6a(k4)) 
	  esactq=esactq*dix+esactq2*(1.-dix) 
	 ENDIF 
	ENDIF  ! loc a 
c 
	IF(iq == 3)THEN 
	 dix2=(rho(l3)-slr)*dfsr(l3) 
	 IF(ip == 3)THEN 
c.....	  eos(i) smoothed in both log(T6) and log(R) 
	  esact=esact*dix2+esactq*(1-dix2) 
	 ENDIF 
	ENDIF 
	IF(esact > 1.d+15)THEN 
	 WRITE(*,'("Interpolation indices out of range; 
	1 please report conditions.")')  
	 STOP 
	ENDIF 
 
	RETURN 
	 
	END SUBROUTINE t6rinteos 
 
c********************************************************************** 
 
	SUBROUTINE r_opal_bin 
	 
c	adaptation à CESAM de la routine readco du package OPAL_EOS 
  
c..... The purpose of this subroutine is to read the data tables 
 
	USE mod_kind 
      	 
	IMPLICIT NONE 
		 
	INTEGER	:: i 
	CHARACTER (len=80) :: chain 
	LOGICAL :: ok 
 
c------------------------------------------------------------------- 
	 
2000	FORMAT(8es10.3) 
 
	IF(itime /= 12345678)THEN 
	 xz=1.d35 ; itime=12345678 
	ENDIF 
	chain=TRIM(nom_chemin)//TRIM(f_eos(1)) 
	INQUIRE(file=TRIM(chain),exist=ok) 
	IF(.NOT.ok)THEN 
	 WRITE(*,10)TRIM(chain) ; WRITE(2,10)TRIM(chain) ; STOP
10	 FORMAT('Fichier de données inconnu : ',a)	  
	ENDIF 
	chain=TRIM(nom_chemin)//TRIM(f_eos(1))	 
	WRITE(*,1)TRIM(chain) ; WRITE(2,1)TRIM(chain) 
1	FORMAT('données ASCII prises dans le fichier: ',a,/, 
	1 'lecture, et c''est long, de ce fichier')	 
	OPEN(unit=60,file=TRIM(chain),form='unformatted') 
	READ(60)xz,t6list,amu_M,alogN_e,rhogr,xin,zz,moles,tmass,frac,icycuse 
	CLOSE(unit=60)
 
	DO i=1,nt 
	 IF(t6list(1,i) /= 0.d0)t6a(i)=t6list(1,i) 
	ENDDO 
    
	DO i=2,nt 
	 dfs(i)=1.d0/(t6a(i)-t6a(i-1)) 
	ENDDO 
	    
	rho(1)=rhogr(1,1) 
	DO i=2,nr 
         rho(i)=rhogr(1,i) 
 	 dfsr(i)=1.d0/(rho(i)-rho(i-1)) 
	endDO   
	DO i=2,mx 
	 dfsx(i)=1.d0/(xx(i)-xx(i-1)) 
	ENDDO 
	 
	PRINT*,'fin de lecture des données EOS opal en binaire' 
       
	RETURN 
 
	END SUBROUTINE r_opal_bin 
 
c*********************************************************************** 
 
	REAL (kind=dp) FUNCTION quadeos(ic,i,x,y1,y2,y3,x1,x2,x3) 
       
c..... 	this function performs a quadratic interpolation. 
 
	USE mod_kind 
	 
	IMPLICIT NONE 
	 
	INTEGER	:: ic, i 
 
	REAL (kind=dp), DIMENSION(30) :: xx12, xx13, xx23, xx1sq, xx1pxx2 
	REAL (kind=dp), DIMENSION(3) ::	xx, yy      
	REAL (kind=dp) :: x, y1, y2, y3, x1, x2, x3, c1, c2, c3 
 
	xx(1)=x1 ; xx(2)=x2 ; xx(3)=x3 
	yy(1)=y1 ; yy(2)=y2 ; yy(3)=y3 
	IF(ic == 0)THEN 
	 xx12(i)=1.d0/(xx(1)-xx(2)) ; xx13(i)=1.d0/(xx(1)-xx(3)) 
	 xx23(i)=1.d0/(xx(2)-xx(3)) ; xx1sq(i)=xx(1)*xx(1) 
	 xx1pxx2(i)=xx(1)+xx(2) 
	ENDIF 
	c3=(yy(1)-yy(2))*xx12(i) ; c3=c3-(yy(2)-yy(3))*xx23(i) 
	c3=c3*xx13(i) ; c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3  
	c1=yy(1)-xx(1)*c2-xx1sq(i)*c3 
	quadeos=c1+x*(c2+x*c3) 
	 
	RETURN 
	 
	END FUNCTION quadeos 
	 
c******************************************************************8 
 
	REAL (kind=dp) FUNCTION gmass(x,z,amoles,eground,fracz,frac) 
	 
	USE mod_kind 
       
	IMPLICIT NONE 
	 
	REAL (kind=dp), PARAMETER, DIMENSION(7) :: 
	1 amas=(/ 0.00054858d0, 20.179d0,15.9994d0, 14.0067d0, 
	2 12.011d0, 4.0026d0, 1.0079d0 /), eion=(/ -3394.873554d0, 
	3 -1974.86545d0, -1433.92718d0, -993.326315d0, 
	4 -76.1959403d0, -15.29409d0, 0.d0 /)		 
	REAL (kind=dp), PARAMETER, DIMENSION(6) :: anum=(/ 10.d0, 8.d0, 
	1 7.d0, 6.d0, 2.d0, 1.d0/) 
	REAL (kind=dp), PARAMETER :: xc=0.247137766d0, xn=0.0620782d0, 
	1 xo=0.52837118d0, xne=0.1624188d0, 
	2 f=xc*amas(5)+xn*amas(4)+xo*amas(3)+xne*amas(2) 	 
	 
	REAL (kind=dp), DIMENSION(7) :: frac 
	 
	REAL (kind=dp) :: x,z,amoles,eground,fracz, 
	1 xc2,xn2,xo2,xh,xtot,xhe,xne2 
	 
	INTEGER	:: i 
		 
c------------------------------------------------------------------------ 
 
	fracz=z/f ; xc2=fracz*xc ; xn2=fracz*xn ; xo2=fracz*xo 
	xne2=fracz*xne ; xh=x/amas(7) ; xhe=(1.-x -z)/amas(6) 
	xtot=xh+xhe+xc2+xn2+xo2+xne2 ; frac(6)=xh/xtot ; frac(5)=xhe/xtot 
	frac(4)=xc2/xtot ; frac(3)=xn2/xtot ; frac(2)=xo2/xtot 
	frac(1)=xne2/xtot ; eground=0.d0 ; amoles=0.d0 
	DO i=1,6 
	 eground=eground+eion(i)*frac(i) 
	 amoles=amoles+(1.d0+anum(i))*frac(i) 
	ENDDO 
	gmass=0.d0 
	DO i=2,7 
	 gmass=gmass+amas(i)*frac(i-1) 
	ENDDO 
 
	RETURN 
       
	END FUNCTION gmass	 
			 
c*********************************************************************** 
 
	SUBROUTINE radsub(t6,density,moles,tmass) 
 
	USE mod_kind 
      	 
 	IMPLICIT NONE 
	 
 	REAL (kind=dp) :: t6,density,moles,tmass,molenak, 
	1 rat,pr,er,sr,st,chir,chitt,pt,et,cvtt,gam3pt,gam1t,gam2pt, 
	2 gam3pt_norad,gam1t_norad,gam2pt_norad	 
 
c---------------------------------------------------------------- 
 
	rat=sigmacc 
	molenak=moles*aprop  ! Mb-cc/unit T6 
	pr=4.d0/3.d0*rat*t6**4   ! Mb  
	er=3.d0*pr/density   ! Mb-cc/gm 
	sr=4.d0/3.d0*er/t6   ! Mb-cc/(gm-unit T6) 
 
c-----	Calculate EOS without radiation correction 
 
	pt=eos(iri(1)) ; et=eos(iri(2)) ; st=eos(iri(3)) 
	chir=eos(iri(6))*eos(iri(1))/pt  
	chitt=(eos(iri(1))*eos(iri(7)))/pt 
	cvtt=(eos(iri(5))*molenak/tmass) 
	gam3pt_norad=pt*chitt/(cvtt*density*t6) 
	gam1t_norad=chir+chitt*gam3pt_norad 
	gam2pt_norad=gam1t_norad/gam3pt_norad 
	 
c---- 	End  no radiation calculation 
 
c---- 	Calculate EOS with radiation calculation 
 
	pt=eos(iri(1))+pr ; et=eos(iri(2))+er ; st=eos(iri(3))+sr 
	chir=eos(iri(6))*eos(iri(1))/pt 
	chitt=(eos(iri(1))*eos(iri(7))+4.d0*pr)/pt 
c	gam1t(jcs,i)=(p(jcs,i)*gam1(jcs,i)+4.d0/3.d0*pr)/pt(jcs,i) 
c	gam2pt(jcs,i)=(gam2p(jcs,i)*p(jcs,i)+4.d0*pr)/pt(jcs,i) 
c	gam3pt(jcs,i)=gam1t(jcs,i)/gam2pt(jcs,i) 
	cvtt=(eos(iri(5))*molenak/tmass+4.*er/t6) 
	gam3pt=pt*chitt/(cvtt*density*t6) 
	gam1t=chir+chitt*gam3pt !eq 16.16 Landau_Lifshitz (Stat. Mech) 
	gam2pt=gam1t/gam3pt  
 
c-----	End Eos calculations with radiation 
 
c     	normalize cvt to 3/2 when gas is ideal,non-degenerate, 
c     	fully-ionized, and has no radiation correction 
c     	cvt=(eos(5)*molenak/tmass+4.*er/t6)/molenak 
	eos(iri(1))=pt ; eos(iri(2))=et ; eos(iri(3))=st 
	eos(iri(5))=cvtt ; eos(iri(6))=chir ; eos(iri(7))=chitt 
c-----	Add difference between EOS with and without radiation.  cvtt 
c       calculation is not accurate enough to give accurate results using 
c       eq. 16.16 Landau&Lifshitz 
	eos(iri(8))=eos(iri(8))+gam1t-gam1t_norad 
	eos(iri(9))=eos(iri(9))+gam2pt-gam2pt_norad 
c	eos(iri(10))=eos(iri(10))+gam3pt-gam3pt_norad 
	 
	RETURN 
	 
	END SUBROUTINE radsub 
	 
c******************************************************************** 
 
	REAL (kind=dp) FUNCTION rhoofp(x,ztab,t6,p,irad,conv) 
	 
	USE mod_kind 
	 
	IMPLICIT NONE	 
 
	INTEGER :: irad, ilo, ihi, imd, mlo, klo, icount 
	 
	REAL (kind=dp) :: rat,pr,t6,p,pnr,ztab,x,pmax,pmin, 
	1 rhog1,rhog2,rhog3,p1,p2,p3,xinit,rinit,tinit 
	 
	LOGICAL, INTENT(out) :: conv 
	LOGICAL, SAVE :: logic=.TRUE. 
	 
c-------------------------------------------------------------- 
	 
	rat=sigmacc ; pr=0.d0 
	IF(irad == 1)pr=4.d0/3.d0*rat*t6**4   ! Mb  
	pnr=p-pr	 
	IF(itime /= 12345678)THEN 
	 xinit=0.5d0 ; tinit=1.d0 ; rinit=1.d-3 
	 CALL esac(xinit,ztab,tinit,rinit,1,0,conv) ; IF(.NOT.conv)RETURN 
	ENDIF 
 
	ilo=2 ; ihi=mx 
8	IF(ihi-ilo > 1)THEN 
	 imd=(ihi+ilo)/2 
	 IF(x <= xa(imd)+1.d-7)THEN 
	  ihi=imd 
	 ELSE 
	  ilo=imd 
	 ENDIF 
	 GOTO 8 
	ENDIF 
	mlo=ilo ; ilo=nt ; ihi=2 
11	IF(ilo-ihi > 1)THEN 
	 imd=(ihi+ilo)/2 
	 IF(t6 == t6list(1,imd))THEN 
	  ilo=imd ; GOTO 14 
	 ENDIF 
	 IF(t6 <= t6list(1,imd))THEN 
	  ihi=imd 
	 ELSE 
	  ilo=imd 
	 ENDIF 
	 GOTO 11 
	ENDIF 
14	klo=ilo 
 
	pmax=xz(mlo,1,klo,nra(klo))*t6*rho(nra(klo)) 
	pmin=xz(mlo,1,klo,1)*t6*rho(1) 
	 
	IF((pnr > 1.25*pmax) .OR. (pnr < pmin))THEN 
c	 PAUSE'logic' 
	 IF(logic)THEN 
	  WRITE(*,15)pnr,pmax,pmin ; WRITE(2,15)pnr,pmax,pmin 
15	  FORMAT('The requested pressure-temperature not in table',/, 
	1 'pnr, pmax, pmin=',3es14.4) 
	  logic=.FALSE. ; conv=.FALSE. ; RETURN 
	 ENDIF  
	ENDIF 
       
	rhog1=rho(nra(klo))*pnr/pmax 
	CALL esac (x,ztab,t6,rhog1,1,0,conv) ; IF(.NOT.conv)RETURN 
	p1=eos(1) 
	IF(p1 > pnr)THEN 
	 p2=p1 ; rhog2=rhog1 ; rhog1=0.2*rhog1 
	 IF(rhog1 < 1.d-14) rhog1=1.d-14 
	 CALL esac(x,ztab,t6,rhog1,1,0,conv) ; IF(.NOT.conv)RETURN 
	 p1=eos(1) 
	ELSE 
	 rhog2=5.*rhog1 
	 IF(rhog2 > rho(klo)) rhog2=rho(klo) 
	 CALL esac(x,ztab,t6,rhog2,1,0,conv) ; IF(.NOT.conv)RETURN 
	 p2=eos(1) 
	ENDIF 
 
	icount=0 
1	CONTINUE 
	icount=icount+1 ; rhog3=rhog1+(rhog2-rhog1)*(pnr-p1)/(p2-p1) 
	CALL esac(x,ztab,t6,rhog3,1,0,conv) ; IF(.NOT.conv)RETURN 
	p3=eos(1) 
	IF(ABS((p3-pnr)/pnr) < 1.d-5)THEN 
	 rhoofp=rhog3 ; conv=.TRUE. 
	 RETURN 
	ENDIF 
	IF(p3 > pnr)THEN 
	 rhog2=rhog3 ; p2=p3 
	 IF(icount < 11)GOTO 1 
	 WRITE(*,1000)t6*1.d6,p*1.d12,x,rhog2 ; conv=.TRUE. 
1000	 FORMAT(/,'!ATTENTION pas de convergence dans rhoofp, t=',es10.3, 
	1 ' p=',es10.3,' X=',es10.3,' ro=',es10.3) 
	 RETURN 
	ELSE 
	 rhog1=rhog3 ; p1=p3 
	 IF(icount < 11)GOTO 1 
	 WRITE(*,1000)t6*1.e6,p*1.e12,x,rhog1 
	 conv=.TRUE. 
	 RETURN 
	ENDIF 
	 
	END FUNCTION rhoofp 
	
	END SUBROUTINE etat_opal5Z
