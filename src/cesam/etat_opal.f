
C***********************************************************************

	SUBROUTINE etat_opal(pp,tt,xchim,deriv,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c routine public du module mod_etat
	
c interface de l'équation d'état OPAL avec CESAM2k	
c package OPAL_EOS
c la recherche de ro par un newton, remplace la routine rhoofp 
c du package OPAL_EOS qui ne converge pas dans certains cas
c appel a EFF en cas de difficultés

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k
c modifs :
c	15 10 97 : variables en double précision
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	14 03 01 : F95 (P. Morel)
c	suppression du blockdata blk_eos pour les initialisations

c entrées :
c	p : pression
c	t : température
c	xchim : composition chimique
c	deriv=.FALSE. : évite le calcul de certaines dérivées

c sortie :
c	ro : densité et dérivées
c	u : énergie interne et dérivées

c----------------------------------------------------------------
	
	USE mod_donnees, only : aradia, f_eos, nchim, nom_chemin, z0
	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(nchim) :: xchim
	REAL (kind=dp), INTENT(in) :: pp, tt	
	LOGICAL, INTENT(in) :: deriv
		
	REAL (kind=dp), INTENT(out) :: ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1
	
c~~~~~initialisations du bloc data blk_eos~~~~~~~~   
	
c	iorder=10  gives all 1st and 2nd order data. See instructions
c	in esac.
c	irad=0 does not add radiation  corrections

c	on exploite les definitions de OPAL, initialement on a:

c       eos(1) is the pressure in megabars (10**12dyne/cm**2)
c       eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
c       eos(3) is the entropy in units of energy/T6
c       eos(4) is dE/dRho at constant T6
c       eos(5) is the specific heat, dE/dT6 at constant V.
c       eos(6) is dlogP/dlogRho at constant T6 Cox and Guili eq 9.82
c       eos(7) is dlogP/dlogT6 at conxtant Rho Cox and Guili eq 9.81
c       eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
c       eos(9) is gamma2/(gaamma2-1). Eqs. 9.88 Cox and Guili
c       eos(10) is gamma3-1. Eqs 9.88 Cox and Guili

c	pour CESAM
c	au lieu de: data(index(i),i=1,10)/1,2,3 ,4,5,6,7,8,9,10/
c	on a mis:   data(index(i),i=1,10)/1,4,10,9,8,2,3,7,5,6/
c	si bien que:
c	eos(1) is the pressure in megabars (10**12dyne/cm**2)
c	eos(2) is dlogP/dlogRho at constant T6. (Ki ro)
c	eos(3) is dlogP/dlogT6 at constant Rho. (Ki T)
c	eos(4) is energy in 10**12 ergs/gm. Zero is zero T6
c	eos(5) is gamma2/(gamma2-1). Eqs. 9.88 CoxGuili (dlnP/dlnT)ad=1/gradad
c	eos(6) is gamma3-1. Eqs 9.88 CoxGuili (dlnT/dlnro)ad
c	eos(7) is gamma1. Eqs. 9.88 CoxGuili (dlnP/dlnro)ad
c	eos(8) is the specific heat, dE/dT6 at constant V.
c	eos(9) is dE/dRho at constant T6
c	eos(10) is the entropy in units of energy/T6

c	on devrait appeler esac avec iorder=7, mais ca ne fonctionne pas ?????

	INTEGER, PARAMETER :: mx=5, mv=10, nr=121, nt=167		
	REAL (kind=dp), DIMENSION(mx,mv,nt,nr) :: xz
        REAL (kind=dp), DIMENSION(nr,nt) :: t6list
        REAL (kind=dp), DIMENSION(nt,nr) :: esk, esk2
	REAL (kind=dp), PARAMETER, DIMENSION(mx) :: xa=(/ 0.d0, 0.2d0,
	1 0.4d0, 0.6d0, 0.8d0 /)	
	REAL (kind=dp), DIMENSION(mx) :: dfsx, zz 
	REAL (kind=dp), DIMENSION(nr) :: dfsr, rho
	REAL (kind=dp), DIMENSION(nt) :: dfs, t6a
	REAL (kind=dp), DIMENSION(mv) :: eos	
	REAL (kind=dp) :: esact
		
	INTEGER	:: m, mf	
	INTEGER, PARAMETER, DIMENSION(10) :: index =(/ 1, 4, 10, 9, 8,
	1 2, 3, 7, 5, 6 /)
	INTEGER, SAVE, DIMENSION(10) :: iri
			
c~~~~~initialisations du bloc data blk_eos~~~~~~~~   
	
	REAL (kind=dp), DIMENSION(nchim) :: Lxchim
	REAL (kind=dp), DIMENSION(mx,nt,nr) :: epl
	REAL (kind=dp), DIMENSION(mx) :: xx	
	REAL (kind=dp), DIMENSION(4) :: h, q
		
	REAL (kind=dp), PARAMETER :: dx=1.d-3, unpdx=1.d0+dx
	REAL (kind=dp), SAVE:: aradias3, a_lim, b_lim, c_lim, d_lim,
	1 ga_lim, g1_lim, text, tlim, ztab
	REAL (kind=dp) :: stor, stor0, dstor, x, t6, r, p12, p, t, plim
	
	INTEGER :: l1, l2, l3, l4, k1, k2, k3, k4, ip, iq
	   
	LOGICAL :: init=.TRUE.
	LOGICAL :: ok
	
c----------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

c T6 is the temperature in units of 10**6 K 
c rho is the density in grams/cc
c R=Rho/T6**3
	p=pp ; t=tt 

	IF(init)THEN
	 init=.FALSE.
	 WRITE(2,1) ; WRITE(*,1)
1	 FORMAT(//,'-------------Equation d''etat-----------------',//,
	1 'Equation d''etat OPAL si dans [ 5010 , 1.e8 ], MHD ou EFF sinon',/,
	2 'der. num.; (p,T)-->(ro,T) calcule par ro_new et non rhoofp',/) 
	 
c appel pour l'interpolation de delta, cp, gradad, alfa, gamma1
c entre OPAL et EFF dans l'atmosphere hors des tables OPAL
	 aradias3=aradia/3.d0
	 tlim=5.d3
	 ztab=z0	!EOS OPAL n'utilise qu'une table en Z
	 t6=tlim*1.d-6 ; x=xchim(1)
	 plim=4.d4 ; p12=plim*1.d-12
	 CALL ro_new(p12,t6,x,ztab,ro,ok)
	 d_lim=eos(3)/eos(2) 				!delta= Ki T/Ki ro
	 ga_lim=1.d0/eos(5)				!gradad	 
	 delta=eos(3)/eos(2)				!Ki T/Ki ro
	 c_lim=p/ro/t*eos(3)*(1.d0/eos(6)+delta)	!cp
	 a_lim=1.d0/eos(2)				!alfa=1/Ki ro
	 b_lim=1.d0-aradias3*t**4/p			!beta
	 g1_lim=eos(7)					!gamma1
	 text=4200.d0	 
	ENDIF
	
c	WRITE(*,2000)p,t,z0

c appel à l'équation d'état
	p12=p*1.d-12 ; x=xchim(1) ; t6=t*1.d-6	
	CALL ro_new(p12,t6,x,ztab,ro,ok)
	IF(ok)THEN
	
c	 WRITE(*,2000)ro,eos ; CALL pause('appel a esac 7')

c après l'appel a ro_new, on utilise les valeurs des eos du COMMON /e/
c pour les grandeurs thermodynamiques
	 drop=ro/p/eos(2)		!1 / dp/dro
	 drot=-ro/t*eos(3)/eos(2)	!- dp/dT / dp/dro	
	 u=eos(4)*1.d12 ; eos(9)=eos(9)*1.e12
	 dup=eos(9)*drop			!du/dT dro/dp
	 dut=eos(9)*drot+eos(8)*1.e6	!du/dro  dro/dT + du/dt	
	 delta=eos(3)/eos(2) !Ki T/Ki ro
	 gradad=1.d0/eos(5) ; cp=p/ro/t*eos(3)*(1.d0/eos(6)+delta)	
	 alfa=1.d0/eos(2)		!1/Ki ro
	 beta=1.d0-aradias3*t**4/p ; gamma1=eos(7)

c dérivées	
	 IF(.NOT.deriv)RETURN
	
c dérivées / P
	 stor0=p ; stor=stor0*unpdx
	 IF(stor < dx)stor=dx
	 dstor=stor-stor0 ; p=stor ; p12=p*1.d-12	
	 CALL ro_new(p12,t6,x,ztab,r,ok)
	 IF(ok)THEN
          drop=(r-ro)/dstor ; dup=(eos(4)*1.d12-u)/dstor
          deltap=(eos(3)/eos(2)-delta)/dstor
          dgradadp=(1.d0/eos(5)-gradad)/dstor
          dcpp=(p/r/t*eos(3)*(1.d0/eos(6)+eos(3)/eos(2))-cp)/dstor
          p=stor0 ; p12=p*1.d-12
        
c dérivées / T
	  stor0=t ; stor=stor0*unpdx
	  IF(stor < dx)stor=dx
	  dstor=stor-stor0 ; t=stor ; t6=t*1.d-6	
	  CALL ro_new(p12,t6,x,ztab,r,ok)
	  IF(ok)THEN
           drot=(r-ro)/dstor ; dut=(eos(4)*1.d12-u)/dstor
           deltat=(eos(3)/eos(2)-delta)/dstor
           dgradadt=(1./eos(5)-gradad)/dstor
           dcpt=(p/r/t*eos(3)*(1.d0/eos(6)+eos(3)/eos(2))-cp)/dstor
           t=stor0 ; t6=t*1.d-6        
        
c dérivées / X
	   stor0=xchim(1) ; stor=stor0*unpdx
	   IF(stor < dx)stor=dx
	   dstor=stor-stor0 ; x=stor
	   CALL ro_new(p12,t6,x,ztab,r,ok)
	   IF(ok)THEN
            drox=(r-ro)/dstor ; dux=(eos(4)*1.d12-u)/dstor
            deltax=(eos(3)/eos(2)-delta)/dstor
            dgradadx=(1./eos(5)-gradad)/dstor
            dcpx=(p/r/t*eos(3)*(1.d0/eos(6)+eos(3)/eos(2))-cp)/dstor
            x=stor0
           ENDIF
	   RETURN
	  ENDIF
	 ENDIF
	ENDIF
	
c en cas de Pb avec opal appel a EFF	
	Lxchim=xchim
	CALL etat_eff(p,t,Lxchim,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c raccord par interpolation lineaire en T entre Tlim et Text
	IF(.FALSE.)THEN
c	IF(t < tlim .AND. t > text)THEN
	 x=(t-text)/(tlim-text)		!x : vt 
	 delta=(1.d0-x)*delta+x*d_lim ; cp=(1.d0-x)*cp+x*c_lim	 
	 gradad=(1.d0-x)*gradad+x*ga_lim ; alfa=(1.d0-x)*alfa+x*a_lim	 
	 beta=(1.d0-x)*beta+x*b_lim ; gamma1=(1.d0-x)*gamma1+x*g1_lim
c	 WRITE(*,2000)x,t,text,tlim,ga_lim,gradad,g1_lim,gamma1
	ENDIF
	
	RETURN
	
	CONTAINS

C***********************************************************************

	SUBROUTINE esac(xh,ztab,t6,r,iorder,irad,ok)
      
c adaptation a CESAM de la routine esac du package OPAL_EOS
c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c	16 10 97 : lecture des tables en binaire

c..... The purpose of this SUBROUTINE is to interpolate 
c      the equation of state and its derivatives in X, T6, density
c        izi=0 recalulate table indices to use; =1 keep previous

c        xh=hydrogen mass fraction
c        ztab is metal fraction of the EOSdata tables you are using.
c           included only for purpose of preventing mismatch
c        t6=T6=temperature in millions of degrees kelvin
c        r=rho=Rho=density(g/cm**3)
c..... to use esac insert COMMON/e/ esact,eos(10) in the CALLing routine.
c      This COMMON CONTAINS the interpolated EOS values for the EOS
c
c..... eos(i) are obtained from a quadradic interpolation at
c      fixed T6 at three values of Rho; followed by quadratic
c      interpolation along T6. Results smoothed by mixing
c      overlapping quadratics.
c         definitions:
c     
c            T6 is the temperature in units of 10**6 K
c 
c            rho is the density in grams/cc
c            R=Rho/T6**3

c            eos(1) is the pressure in megabars (10**12dyne/cm**2)
c            eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
c            eos(3) is the entropy in units of energy/T6
c            eos(4) is dE/dRho at constant T6
c            eos(5) is the specific heat, dE/dT6 at constant V.
c            eos(6) is dlogP/dlogRho at constant T6. 
c                   Cox and Guil1 eq 9.82
c            eos(7) is dlogP/dlogT6 at conxtant Rho.
c                   Cox and Guil1 eq 9.81
c            eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
c            eos(9) is gamma2/(gaamma2-1). Eqs. 9.88 Cox and Guili
c            eos(10) is gamma3-1. Eqs 9.88 Cox and Guili

c            iorder sets maximum index for eos(i);i.e., iorder=1
c                   gives just the pressure
c
c            irad  if =0 no radiation correction; if =1 adds radiation

c            index(i),i=1,10  sets order in which the equation of state
c            variables are stored in eos(i).  Above order corresponds
c            to block data statement:
c                 data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/. 
c            If you, for example, only want to RETURN gamma1: set iorder=1 
c            and set: data (index(i),i=1,10)/8,2,3,4,5,6,7,1,9,10/

c------------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), DIMENSION(7) :: frac
	REAL (kind=dp), PARAMETER :: aprop=83.1446304d0
	REAL (kind=dp) :: moles, xh, ztab, t6, r, slt, dixr, p0,
	1 slr, z, sum1, sum2, sum23, sum33, tmass, eground, fracz

	INTEGER	iorder,irad,i,j,ilo,ihi,imd,mg,mh,mi,mf2,ms,kf2,ir,it,is,iw,iv
	
	LOGICAL :: ok, logic1=.TRUE., logic2=.TRUE., init=.TRUE.
	      
	SAVE

c-------------------------------------------------------------------------

	IF(iorder > 10 )THEN
         WRITE (*,'(" esac : Iorder cannot exceed 10")') ; STOP
	ENDIF
	IF((irad /= 0) .AND. (irad /= 1))THEN
	 WRITE (*,'(" esac : Irad must be 0 or 1")')
	 STOP
	ENDIF
	slt=t6 ; slr=r

	IF(init)THEN
	 init=.FALSE.
	 DO i=1,10
	  DO j=1,10
	   IF(index(i) == j)iri(i)=j
	  ENDDO
	 ENDDO
	 DO  i=1,mx
	  xx(i)=xa(i)
	 ENDDO

c READ the data files
c	 CALL readco		!en ASCII
	 CALL r_opal_bin	!en binaire
         z=zz(1)
         IF(abs(ztab-z) > 0.03)THEN
          WRITE(*,10)ztab,z ; WRITE(2,10)ztab,z
10	  FORMAT(/,'!!ATTENTION, Z=',es10.3,
	1 ' est très different du Z de la table d''equation d''etat, Z=',
	2 es10.3,/,'!!ATTENTION',/)
	 ENDIF
	ENDIF               

	IF(z+xh-1.d-6 > 1.d0 )THEN
	 WRITE(*,'(" esac : Mass fractions exceed unity")')
	 PRINT*,'Z=',z,' X=',xh
	 STOP
	ENDIF

c Determine T6,rho grid points to use in the interpolation.
c	PRINT*,t6a(1),t6a(nt)
c	PRINT*,rho(1),rho(nr)
c	PAUSE'limites'	

	slt=t6 ; slr=r
	
	IF((slt > t6a(1)) .OR. (slt < t6a(nt)) .OR. 
	1  (slr < rho(1)) .OR. (slr > rho(nr)))THEN 
	 IF(logic1)THEN
	  WRITE(*,21)slt,slr
	  WRITE(2,21)slt,slr
21	  FORMAT(/'!!!!!ATTENTION!!!!!!!!!!',/
	1 'on a recontre au moins 1 fois le Pb. suivant dans la',/,
	2 'routine "esac" de livermore, on utilise MHD  ou EFF',/, 
	3 'T6 ou LogR outside of table range, T6=',es10.3,
	4 ' LogR=',es10.3,/,
c	5 'raccord par interp. lin. en T dans atmosphere',/,
	6 '!!!!!ATTENTION!!!!!!!!!!',/)		
	  logic1=.FALSE.
	 ENDIF
	 ok=.FALSE.
	 RETURN
	ENDIF

c	IF(izi == 0)THEN  ! freeze table indices if not 0
	ilo=2
	ihi=mx
8  	IF(ihi-ilo > 1)THEN
	 imd=(ihi+ilo)/2
	 IF(xh <= xa(imd)+1.e-7)THEN
	  ihi=imd
	 ELSE
	  ilo=imd
	 ENDIF
	 GOTO 8
	ENDIF
	i=ihi
	mf=i-2
	mg=i-1
	mh=i
	mi=i+1
	mf2=mi
	IF(xh < 1.e-6)THEN
	 mf=1
	 mg=1
	 mh=1
	 mi=2
	 mf2=1
	ENDIF
        IF((xh <= xa(2)+1.e-7) .OR. (xh >= xa(mx-2)-1.e-7))mf2=mh

	ilo=2
	ihi=nr
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
13	i=ihi
	l1=i-2
	l2=i-1
	l3=i
	l4=l3+1

	ilo=nt
	ihi=2
11	IF(ilo-ihi > 1)THEN
	 imd=(ihi+ilo)/2
	 IF(t6 == t6list(1,imd))THEN
	  ilo=imd
	  GOTO 14
	 ENDIF
	 IF(t6 <= t6list(1,imd))THEN
	  ihi=imd
	 ELSE
	  ilo=imd
	 ENDIF
	 GOTO 11
	ENDIF
14	i=ilo
	k1=i-2
	k2=i-1
	k3=i
	k4=k3+1
	kf2=k4
	IF(kf2 > nt)kf2=k3
	IF(k3 == 0)THEN
	 WRITE (*,'(" ihi,ilo,imd",3i5)')
	ENDIF
c	ENDIF		!izi

c check to determine if interpolation indices fall within
c table boundaries. choose largest allowed size.
	sum1=0.0
	sum2=0.0
	sum23=0.0
	sum33=0.0
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
	 DO ir=l1,l1+2
	  DO it=k1,kf2
	   sum23=sum23+xz(m,1,it,ir)
	  ENDDO
	 ENDDO
	 DO ir=l1,l1+3
	  DO it=k1,kf2
	   sum33=sum33+xz(m,1,it,ir)
	  ENDDO
	 ENDDO
	ENDDO
	iq=2
	ip=2
	IF(sum2 > 1.e+30)THEN
	 IF(sum1 < 1.e+25)THEN
	  k1=k3-3
	  k2=k1+1
	  k3=k2+1
	  l1=l3-3
	  l2=l1+1
	  l3=l2+1
	  GOTO 15
	 ELSE
	  IF(logic2)THEN
	   WRITE(*,22)xh,t6*1.d6,r
	   WRITE(2,22)xh,t6*1.d6,r	 
22	   FORMAT(/'!!!!!ATTENTION!!!!!!!!!!',/
	1  'on a recontre au moins 1 fois le Pb. suivant dans la',/,
	2  'routine "esac" de livermore, on utilise MHD ou EFF',/, 
	3  'T6/log rho in empty region od table (65)',/,
	4  'X=',es10.3,' T=',es10.3,' ro=',es10.3,/,
	5  '!!!!!ATTENTION!!!!!!!!!!',/)
	   logic2=.FALSE.
	  ENDIF
	  ok=.FALSE.
	  RETURN
	 ENDIF
	ENDIF
	
	IF(sum23 < 1.e+30) ip=3
	IF(sum33 < 1.e+30) iq=3
	IF(t6 >= t6list(1,2)+1.e-7) ip=2
	IF(slr <= rho(2)+1.e-15) iq=2
	IF(l3 == nr) iq=2
	
15	CONTINUE

	DO iv=1,iorder
	 DO m=mf,mf2
	  is=0
	  DO ir=l1,l1+iq
	   DO it=k1,k1+ip
	    epl(m,it,ir)=xz(m,iv,it,ir)
	    is=1
	   ENDDO
	  ENDDO
	 ENDDO
	 IF(zz(mg) /= zz(mf) .OR. zz(mh) /= zz(mf))THEN
	  WRITE(*,'("Z does not match Z in EOSdata files you are using")')
	  STOP
	 ENDIF
	 IF(z /= zz(mf)) GOTO 66
	 is=0
	 iw=1
	 DO ir=l1,l1+iq
	  DO it=k1,k1+ip
	   IF(mf2 == 1)THEN
	    esk(it,ir)=epl(mf,it,ir)
	   ELSE
	    esk(it,ir)=quad(is,iw,xh,epl(mf,it,ir),epl(mg,it,ir),
	1   epl(mh,it,ir),xx(mf),xx(mg),xx(mh))
	    IF(esk(it,ir) > 1.e+20)THEN
	     WRITE(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir,l3,k3,iq,ip
	     WRITE(*,'(3e12.4)')  (epl(ms,it,ir),ms=mf,mf+2)
	    ENDIF	!esk
	    is=1
	   ENDIF	!mf2
          ENDDO	!it
	 ENDDO	!ir
	 IF(mi == mf2)THEN  ! interpolate between quadratics
	  is=0
	  iw=1
	  dixr=(xx(mh)-xh)*dfsx(mh)
	  DO ir=l1,l1+iq
	   DO it=k1,k1+ip
	    esk2(it,ir)=quad(is,iw,xh,epl(mg,it,ir),epl(mh,it,ir),
	1   epl(mi,it,ir),xx(mg),xx(mh),xx(mi))
	    IF(esk(it,ir) > 1.e+20)THEN
	     WRITE(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir,l3,k3,iq,ip
	     WRITE(*,'(3e12.4)')  (epl(ms,it,ir),ms=mg,mg+2)
	    ENDIF
	    esk(it,ir)=esk(it,ir)*dixr+esk2(it,ir)*(1.-dixr)
	    is=1
	   ENDDO	!it
	  ENDDO		!ir
	 ENDIF		!mi
	 is=0

c	 completed X interpolation. Now interpolate T6 and rho on a
c	 4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
c	 mixes overlapping quadratics to obtain smoothed derivatives.

	 CALL t6rinterp(slr,slt)
	 eos(iv)=esact
	ENDDO		!iv

	p0=t6*r
	eos(iri(1))=eos(iri(1))*p0   ! interpolated in p/po
	eos(iri(2))=eos(iri(2))*t6   ! interpolated in E/T6
	eos(iri(4))=eos(iri(4))/sqrt(r*t6) ! interp dE/dr/sqrt(r/t6)
	tmass=gmass(xh,z,moles,eground,fracz,frac)
	IF(irad == 1)THEN
	 CALL radsub(t6,r,moles,tmass)
	ELSE
	 eos(iri(5))=eos(iri(5))*moles*aprop/tmass
	ENDIF
	ok=.true.
	
	RETURN

66	WRITE(*,'(" Z does not match Z in EOSdata* files you are",
	1 "dans routine esac, using (66)")')
	WRITE (*,'("mf,zz(mf)=",i5,e12.4)') mf,zz(mf)
	WRITE (*,'("  iq,ip,k3,l3,xh,t6,r,z= ",4i5,4e12.4)')
	1 ip,iq,k3,l3,xh,t6,r,z
	STOP

	END SUBROUTINE esac

C***********************************************************************

	REAL (kind=dp) FUNCTION gmass(x,z,amoles,eground,fracz,frac)
      
c	adaptation a CESAM de la routine gmass du package OPAL_EOS

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	15 10 97 : variables en double precision

c----------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), PARAMETER, DIMENSION(7) ::
	1 amas=(/ 0.00054858d0, 20.179d0, 15.9994d0, 14.0067d0,
	2 12.011d0, 4.0026d0, 1.0079d0 /)		
	REAL (kind=dp), PARAMETER, DIMENSION(6) ::
	1 anum=(/ 10.d0, 8.d0, 7.d0, 6.d0, 2.d0, 1.d0/),
	2 eion=(/ -3388.637d0,-1970.918d0,-1431.717d0,
	3 -993.2303d0,-76.2315d0,-15.29409d0 /)
	REAL (kind=dp), PARAMETER :: xc=0.247137766d0, xn=0.0620782d0,
	1 xo=0.52837118d0, xne=0.1624188d0,
	2 f=xc*amas(5)+xn*amas(4)+xo*amas(3)+xne*amas(2) 	
	REAL (kind=dp), DIMENSION(7) :: frac
	REAL (kind=dp) :: x,z,amoles,eground,fracz,
	1 xc2,xn2,xo2,xh,xtot,anume,xhe,xne2
	
	INTEGER	:: i
		
c------------------------------------------------------------------------

	fracz=z/f     
	xc2=fracz*xc
	xn2=fracz*xn
	xo2=fracz*xo
	xne2=fracz*xne 
	xh=x/amas(7)
	xhe=(1.d0-x -z)/amas(6)
	xtot=xh+xhe+xc2+xn2+xo2+xne2
	frac(6)=xh/xtot
	frac(5)=xhe/xtot
	frac(4)=xc2/xtot
	frac(3)=xn2/xtot
	frac(2)=xo2/xtot
	frac(1)=xne2/xtot
	eground=0.0d0
	amoles=0.0d0
	DO i=1,6
	 eground=eground+eion(i)*frac(i)
	 amoles=amoles+(1.d0+anum(i))*frac(i)
	ENDDO
	anume=amoles-1.d0
	gmass=anume*amas(1)
	DO i=2,7
	 gmass=gmass+amas(i)*frac(i-1)
	ENDDO
      
	RETURN

	END FUNCTION gmass	

C***********************************************************************

	REAL (kind=dp) FUNCTION quad(ic,i,x,y1,y2,y3,x1,x2,x3)
	
c	adaptation a CESAM de la routine quad du package OPAL_EOS
 
c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.

c	15 10 97 : variables en double precision

c..... this FUNCTION performs a quadratic interpolation.

c-------------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: ic,i
	REAL (kind=dp), INTENT(in) :: x, y1, y2, y3, x1, x2, x3
	
	REAL (kind=dp), DIMENSION(30) :: xx12, xx13, xx23, xx1sq, xx1pxx2
	REAL (kind=dp), DIMENSION(3) :: xx, yy	
	
	REAL (kind=dp) :: c1, c2, c3
	
c-----------------------------------------------------------------------    
     
	xx(1)=x1 ; xx(2)=x2 ; xx(3)=x3
	yy(1)=y1 ; yy(2)=y2 ; yy(3)=y3
	IF(ic == 0)THEN
	 xx12(i)=1./(xx(1)-xx(2))
	 xx13(i)=1./(xx(1)-xx(3))
	 xx23(i)=1./(xx(2)-xx(3))
	 xx1sq(i)=xx(1)*xx(1)
	 xx1pxx2(i)=xx(1)+xx(2)
	ENDIF
	c3=(yy(1)-yy(2))*xx12(i)
	c3=c3-(yy(2)-yy(3))*xx23(i)
	c3=c3*xx13(i)
	c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3
	c1=yy(1)-xx(1)*c2-xx1sq(i)*c3
	quad=c1+x*(c2+x*c3)
	
	RETURN
      
	END FUNCTION quad

C***********************************************************************

	SUBROUTINE r_opal_bin
	
c adaptation a CESAM de la routine readco du package OPAL_EOS
c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c	15 10 97 : variables en double precision 
    
c..... The purpose of this SUBROUTINE is to READ the data tables

c-------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), DIMENSION(mx,nr) :: rhogr	
	REAL (kind=dp), DIMENSION(nr,nt) :: alogr	
	REAL (kind=dp), DIMENSION(mx,6) :: frac		
	REAL (kind=dp), DIMENSION(mx) :: moles, tmass, xin	

	INTEGER, DIMENSION(mx,nr) :: icycuse
	INTEGER	:: i

	LOGICAL :: ok
	CHARACTER(len=80) :: chain
     
	SAVE
 	
c------------------------------------------------------------------
 	
2000	FORMAT(8es10.3)

	chain=TRIM(nom_chemin)//f_eos(1)
	INQUIRE(file=TRIM(chain)//'.gz',exist=ok)
	IF(ok)CALL SYSTEM('gunzip '//TRIM(chain)//'.gz')			
	WRITE(*,1)chain ; WRITE(2,1)chain
1	FORMAT('données prises dans le fichier binaire:',/,a80)	
	OPEN(unit=60,file=TRIM(chain),form='unformatted')
	PRINT*,'lecture, des données EOS opal, et c''est long'
	PRINT*,SIZE(xz),mx,mv,nt,nr,mx*mv*nt*nr
	READ(60)xz,t6list,alogr,rhogr,xin,zz,moles,tmass,frac,icycuse
	CLOSE(unit=60)
	WRITE(*,2)TRIM(chain)
2	FORMAT('recompression du fichier binaire : ',a)	
	CALL system('gzip '//TRIM(chain))
c	WRITE(*,2000)(zz(i),i=1,mx)
c	PAUSE'lecture du binaire'
			
	DO i=1,nt
	 IF(t6list(1,i) /= 0.d0)THEN
	  t6a(i)=t6list(1,i)
	 ENDIF	 
	ENDDO
   
	DO i=2,nt
	 dfs(i)=1.d0/(t6a(i)-t6a(i-1))
	ENDDO
	   
	rho(1)=rhogr(1,1)
	DO i=2,nr
         rho(i)=rhogr(1,i)
	 dfsr(i)=1.d0/(rho(i)-rho(i-1))
	ENDDO  
	DO i=2,mx
	 dfsx(i)=1.d0/(xx(i)-xx(i-1))
	ENDDO
	
	PRINT*,'fin de lecture des donnees EOS opal en binaire'
	
	RETURN

	END SUBROUTINE r_opal_bin
C
C***********************************************************************
C
	SUBROUTINE radsub(t6,density,moles,tmass)
      
c	adaptation a CESAM de la routine radsub du package OPAL_EOS
 
c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.

c	15 10 97 : variables en double precision

	USE mod_kind

 	IMPLICIT NONE

  	REAL (kind=dp) :: t6, density, moles, tmass, molenak,
 	1 unitf=0.9648575d0, unitfold=0.9648575d0, sigmacc=1.8914785d-3,
	2 aprop=83.14511d0	
	
	REAL (kind=dp) :: rat,pr,er,sr,revise,fixerror,st,chir,chitt,pt,et,
	1 cvtt,gam3pt,gam1t,gam2pt,dedrhoat
	 	   
	SAVE 
	     
c------------------------------------------------------------------------

cPhysical constants
c       Na=6.0221367e+23
c       k =1.380658e-16 !   erg/degree K
c       Na*k=6.0221367E+23*1.380658e-16 erg/degree K=8.314511E+7 erg/degree K
c           =8.314511e+7*11604.5 erg/eV=0.9648575E+12 erg/eV
c           Define unitf= Na*k/e+12=0.9648575
c           unitf=0.9648575  ! obtained with latest physical constants
c           unitfold=0.9652   ! old units- still used in the EOS code
c           In these units energy/density is in units of Mb-CC/gm
c           Pressure is in units of E+12 bars=Mb
c       sigma is the Stefan-Boltzmann constant
c       sigma=5.67051E-5 !   erg /(s*cm**2*K**4)
c       c=2.99792458E+10 !   cm/sec

c     rat=sigma/c    ! dyne/(cm**2*K**4)

c     rat=rat*1.e+24  !  Convert degrees K to units 10**6 K (T replaced with T6)
      rat=sigmacc

      pr=4./3.*rat*t6**4   ! Mb 
      er=3.*pr/density   ! Mb-cc/gm
      sr=4./3.*er/t6   ! Mb-cc/(gm-unit T6)
      revise=unitf/unitfold
      eos(iri(1))=eos(iri(1))*revise
      eos(iri(2))=eos(iri(2))*revise
      eos(iri(3))=eos(iri(3))*revise
      eos(iri(4))=eos(iri(4))*revise
      eos(iri(5))=eos(iri(5))*revise
      pt=eos(iri(1))+pr
      et=eos(iri(2))+er
      fixerror=1./(.0114045*.0116045) !4/14/96 Corrects for earlier
c        multiplication by .0114045 where divide by .0116045 needed
c        Converts eV to T6
      st=eos(iri(3))*fixerror+sr
      chir=eos(iri(6))*eos(iri(1))/pt
      chitt=(eos(iri(1))*eos(iri(7))+4.*pr)/pt
c     gam1t(jcs,i)=(p(jcs,i)*gam1(jcs,i)+4./3.*pr)/pt(jcs,i)
c     gam2pt(jcs,i)=(gam2p(jcs,i)*p(jcs,i)+4.*pr)/pt(jcs,i)
c     gam3pt(jcs,i)=gam1t(jcs,i)/gam2pt(jcs,i)
      molenak=moles*aprop  ! Mb-cc/unit T6
      cvtt=(eos(iri(5))*molenak/tmass+4.*er/t6)
      gam3pt=pt*chitt/(cvtt*density*t6)
      gam1t=chir+chitt*gam3pt
      gam2pt=gam1t/gam3pt

c     normalize cvt to 3/2 when gas is ideal,non-degenerate,
c     fully-ionized, and has no radiation correction
c     cvt=(eos(5)*molenak/tmass+4.*er/t6)
c    x  /molenak
      dedrhoat=eos(iri(4))-er/density
      eos(iri(1))=pt
      eos(iri(2))=et
      eos(iri(3))=st
      eos(iri(4))=dedrhoat
      eos(iri(5))=cvtt
      eos(iri(6))=chir
      eos(iri(7))=chitt
      eos(iri(8))=gam1t
      eos(iri(9))=gam2pt
      eos(iri(10))=gam3pt
      
      RETURN
      
      END SUBROUTINE radsub
C
C***********************************************************************
C
	SUBROUTINE ro_new(p,t6,x,ztab,ro,conv)

c	package OPAL_EOS
c	recherche de ro par un newton, remplace la routine rhoofp 
c	du package OPAL_EOS qui ne converge pas dans certains cas solaires

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	15 10 97 : variables en double précision

c entrées
c	p: pression
c	t6: temperature/1.d6
c	x: X hydrogene
c	ztab: Z

c sorties:
c	ro: densite

c----------------------------------------------------------------------

	USE mod_donnees, only : granr, ah, ahe4
	USE mod_kind	

	IMPLICIT NONE
	
c	INTEGER, PARAMETER :: mv=10
	
	LOGICAL, INTENT(out) :: conv    
		
	REAL (kind=dp), PARAMETER :: epsi=1.d-5

	REAL (kind=dp) :: granrtp
	REAL (kind=dp) :: p, t6, x, ro, lr, lp, mum1, ztab, rop,
	1 pp=-1.d0, tp=-1.d0, xp=-1.d0
	
	INTEGER, PARAMETER :: nmax=20
	INTEGER :: n
	
	LOGICAL, SAVE :: init=.true.
	
c	REAL (kind=dp) :: esact, eos(mv)	
c	COMMON/e/esact,eos
	
c-------------------------------------------------------------------

2000	FORMAT(1x,8es10.3)

	IF(init)THEN
	 granrtp=granr*1.d6*1.d-12	!1.e6 pour T6, 1.e-12 pour p
	 init=.FALSE.
	ENDIF
	
c	PRINT*,index

c	initialisation de ro

	IF(max( abs(pp-p)/p,abs(tp-t6)/t6,abs(xp-x) ) < 0.1d0)THEN
	 ro=rop
	ELSE
	 mum1=2.d0*x/ah+3.d0*(1.-x-ztab)/ahe4	!totalement ionise
	 ro=p/granrtp/mum1/t6	!ro estime
	ENDIF
	 
	lp=log(p)
	lr=log(ro)
	
	n=0
	conv=.FALSE.
	DO WHILE((.NOT.conv .AND. n < nmax) .OR. n < 1)
	 CALL esac(x,ztab,t6,ro,mv,1,conv)	!calcul
	 IF(.NOT.conv)RETURN	!eos(1)=P*1.d-12
	 
	 lr=lr-(log(eos(1))-lp)/eos(2)	!eos(2)=d ln P / d ln ro
	 n=n+1
	 
	 conv=abs(eos(1)-p)/p < epsi	 
	 
	 ro=exp(lr)
	 
c	 WRITE(*,2000)p,ro,(eos(1)-p)/p	 
c	 PRINT*,conv,n
c	 PAUSE	 
	ENDDO
	
	IF(.NOT.conv)THEN
	 WRITE(*,10)t6*1.e6,p*1.e12,x,ro
10	 FORMAT(/,'!ATTENTION pas de convergence dans ro_new, t=',es10.3,
	1 ' p=',es10.3,' X=',es10.3,' ro=',es10.3)
	 pp=-1.d0
	 tp=-1.d0
	 xp=-1.d0
	ELSE
	 pp=p
	 tp=t6
	 xp=x
	 rop=ro
	ENDIF	
	
c	WRITE(*,2000)lr

	RETURN

	END SUBROUTINE ro_new
C
C***********************************************************************
C
	SUBROUTINE t6rinterp(slr,slt)
      
c	adaptation de la routine t6rinterp du package OPAL_EOS

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k
     
c	The purpose of this SUBROUTINE is to interpolate in T6 and rho

c---------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	INTEGER	iu,is,kx,iw
	
	REAL (kind=dp) :: slr, slt, dix, dix2, esact2, esactq, esactq2
      
	SAVE
	
c------------------------------------------------------------------------
	
      iu=0
      is=0

      DO kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slr,esk(kx,l1),esk(kx,l2),esk(kx,l3),
     x  rho(l1),rho(l2),rho(l3))
          IF(iq == 3)THEN
            iw=2
            q(iu)=quad(is,iw,slr,esk(kx,l2),esk(kx,l3),esk(kx,l4),
     x      rho(l2),rho(l3),rho(l4))
          ENDIF
        is=1
      ENDDO
c
      is=0 ; iw=1
c..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      esact=quad(is,iw,slt,h(1),h(2),h(3),t6a(k1),t6a(k2),t6a(k3))
        IF(iq == 3)THEN
c.....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          esactq=quad(is,iw,slt,q(1),q(2),q(3),t6a(k1),t6a(k2),t6a(k3))
        ENDIF
        IF(ip == 3)THEN
c.....    eos(i) in lower-left 3x3.
          esact2=quad(is,iw,slt,h(2),h(3),h(4),t6a(k2),t6a(k3),t6a(k4))
c.....    eos(i) smoothed in left 3x4
          dix=(t6a(k3)-slt)*dfs(k3)
          esact=esact*dix+esact2*(1.-dix)
c       ENDIF   ! moved to loc a
        IF(iq == 3)THEN
 
c.....     eos(i) in upper-right 3x3.
          esactq2=quad(is,iw,slt,q(2),q(3),q(4),t6a(k2),t6a(k3),t6a(k4))
          esactq=esactq*dix+esactq2*(1.-dix)
        ENDIF
        ENDIF  ! loc a
c
        IF(iq == 3)THEN
          dix2=(rho(l3)-slr)*dfsr(l3)
            IF(ip == 3)THEN
c.....        eos(i) smoothed in both log(T6) and log(R)
              esact=esact*dix2+esactq*(1-dix2)
            ENDIF
        ENDIF
        IF(esact > 1.e+15)THEN
          WRITE(*,'("Interpolation indices out of range",
     x              ";please report conditions.")') 
          STOP
        ENDIF
	
	RETURN
      
      END SUBROUTINE t6rinterp
      
      END SUBROUTINE etat_opal
