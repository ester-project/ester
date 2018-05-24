
c**************************************************************************

	SUBROUTINE ppcno9Fe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c	routine private du module mod_nuc

c	CYCLEs PP, CNO 
c	cf. Clayton p. 380, 392 et 430

c	éléments pris en compte:
c	H1, He3, He4, C12, C13, N14, N15, O16, O17, Ex, Fe
c	Ex est l'élément fictif complément, il n'intéresse que la diffusion
c	Fe n'intéresse que la diffusion
c	H2, Li7, Be7 à l'équilibre

c	au premier appel une table ppcno est crée automatiquement par
c	rq_reac --> tabul_nuc

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées :
c	t : température cgs
c	ro : densité cgs
c	comp : abondances
c	deriv=.true. : on calcule le jacobien
c	fait=1 : initialisation de la composition chimique
c	    =2 : calcul de dcomp et jacobien si deriv
c	    =3 : énergie nucléaire et dérivées / t et ro
c	    =4 : production de neutrinos

c sorties
c	dcomp : dérivée temporelle (unité de temps : 10**6 ans)
c	jac : jacobien (unité de temps : 10**6 ans)
c	epsilon, et, ero, ex : énergie thermonucléaire (unité de temps : s)
c			   : et dérivées /t, ro ,X
c	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
c	hhe réaction : H1(p,e+ nu)H2
c	be7e réaction : Be7(e-,nu g)Li7
c	b8e réaction : B8(,e+ nu)Be8
c	n13e réaction : N13(,e+ nu)C13
c	o15e réaction : O15(e+,nu)N15 
c	f17e réaction : F17(,e+ nu)O17

c initialisation de
c	ab_min : abondances négligeables
c	ab_ini : abondances initiales

c	r(1) : réaction H1(p,e+ nu)H2			PP
c	r(2) : réaction H2(p,g)H3
c	r(3) : réaction He3(He3,2p)He4
c	r(4) : réaction He3(a,g)Be7
c	r(5) : réaction Li7(p,a)He4
c	r(6) : réaction Be7(e-,nu g)Li7
c	r(7) : réaction Be7(p,g)B8(,e+ nu)Be8(a)He4

c	r(8) : réaction C12(p,g)N13(,e+ nu)C13		CNO
c	r(9) : réaction C13(p,g)N14
c	r(10) : réaction N14(p,g)O15(e+,nu)N15
c	r(11) : réaction N15(p,g)O16
c	r(12) : réaction N15(p,a)C12
c	r(13) : réaction O16(p,g)F17(,e+ nu)O17
c	r(14) : réaction O17(p,a)N14

c	indices des éléments

c	H1 : 1
c	He3 : 2
c	He4 : 3
c	C12 : 4
c	C13 : 5
c	N14 : 6
c	N15 : 7
c	O16 : 8
c	O17 : 9
c	Ex : 10
c	Fe : 11

c----------------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon, ife56, ihe4,
	1 i_ex, langue, nchim, nom_elem, nom_xheavy,
	2 nucleo, secon6, t_inf, x0, y0, zi, z0
	USE mod_kind
	USE mod_numerique, ONLY : gauss_band
		
	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: fait
	LOGICAL, INTENT(in) :: deriv
	REAL (kind=dp), INTENT(in):: t, ro
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: comp
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: jac	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dcomp, ex, epsilon
	REAL (kind=dp), INTENT(out) :: et, ero, hhe, be7e, b8e, n13e,
	1 o15e, f17e
	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: drx, dqx
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a, b
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: anuc, comp_dex,
	1 dmuex, dh2x, denx, dbe7x, dli7x, drt, dro, r, q, dqt, dqo		
	REAL (kind=dp) :: mue, nbz, den, be7, h2, li7, dh2t, dh2ro,
	1 dent, denro, dbe7t, dbe7ro, dli7t, dli7ro,
	2 mass_ex, charge_ex, sum_a
		
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc
	INTEGER :: i, j
	
	LOGICAL :: inversible
	
	CHARACTER (len=2) :: text
		
c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(5es15.8)
2002	FORMAT(11es8.1)
	
c	initialisations

	SELECT CASE(fait)
	CASE(0)
	 
c	 définition de nchim: nombre d'éléments chimiques dont on
c	 calcule l'abondance H1, He3, He4, C13, C13, N14, N15, O16, O17,
c	 Ex, Fe

	 nchim=9+2

c	 appel d'initialisation pour tabulation des réactions nucléaires
c	 allocations fictives

	 ALLOCATE(drx(1,1),dqx(1,1),r(1),drt(1),dro(1),q(1),
	1 dqt(1),dqo(1),dmuex(1))
	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex) 
c	 PAUSE'après case0'
	 
	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 
	CASE(1)

c	 détermination des abondances initiales
c	 He3+He4=Y0
c	 Z0 = somme des éléments plus lourds que helium
c	 dans Z rapports en nombre

	 CALL abon_ini
	 
c	 on ajoute Fe56

	 nucleo(nchim)=m(26)	!nucleo de Fe
	 zi(nchim)=c(26)	!charge de Fe
	 WRITE(text,10)NINT(m(26))
	 i=NINT(c(26))
	 nom_elem(nchim)=elem(i)//text 	!nom de Fe
	  
c	 Ex : élément fictif moyenne des éléments # CNO

	 charge_ex=0.d0 ; mass_ex=0.d0 ; sum_a=0.d0
	 B1: DO i=3,nelem_ini		!à partir de Li
	  IF(elem(i) == ' C')CYCLE B1
	  IF(elem(i) == ' N')CYCLE B1
	  IF(elem(i) == ' O')CYCLE B1
	  IF(elem(i) == 'Fe')CYCLE B1	  
	  charge_ex=charge_ex+c(i)*ab(i) ; mass_ex=mass_ex+m(i)*ab(i)
	  sum_a=sum_a+ab(i)
	 ENDDO B1
	 charge_ex=NINT(charge_ex/sum_a) ; mass_ex=mass_ex/sum_a
	 WRITE(text,10)NINT(mass_ex)
10	 FORMAT(i2)

c élément fictif
	 nucleo(nchim-1)=mass_ex	!nucleo de l'élément chimique reliquat
	 zi(nchim-1)=charge_ex	!charge de l'élément chimique reliquat
	 i=NINT(charge_ex)
	 nom_elem(nchim-1)=elem(i)//text	!nom reliquat
	 nom_xheavy=nom_elem(nchim-1)
	 i_ex=nchim-1
 	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1023)TRIM(nom_elem(nchim-1)),NINT(mass_ex),NINT(charge_ex)
	  WRITE(2,1023)TRIM(nom_elem(nchim-1)),NINT(mass_ex),NINT(charge_ex)	 
1023	  FORMAT(a,': fictitious species /= CNO, of mass : ',i3,/,
	1 'and charge :',i3)	 
	 CASE DEFAULT	 
	  WRITE(*,23)TRIM(nom_elem(nchim-1)),NINT(mass_ex),NINT(charge_ex)
	  WRITE(2,23)TRIM(nom_elem(nchim-1)),NINT(mass_ex),NINT(charge_ex)	 
23	  FORMAT(a,': élément fictif /= CNO, de masse : ',i3,/,
	1 'et de charge :',i3)
	 END SELECT	 	 
c	 PRINT*,nchim ; WRITE(*,2000)nucleo	 
	 
c détermination des abondances initiales, a(équation,élément)

	 ALLOCATE(a(nchim,nchim),indpc(nchim),b(1,nchim))
	 a=0.d0 ; b=0.d0 ; indpc=1	
		
	 a(1,1)=nucleo(1)  	!H1
	 b(1,1)=x0
	
	 a(2,2)=nucleo(2)  	!He3
	 a(2,3)=nucleo(3)	!He4
	 b(1,2)=y0

	 DO j=4,nchim
	  a(3,j)=nucleo(j)	!somme i > 4 comp(i)*nucleo(i)=Z0
	  a(4,j)=-abon_rela(6)	!somme comp(i) C, C/Z
	  a(5,j)=-abon_rela(7)	!somme comp(i) N, N/Z	 
	  a(6,j)=-abon_rela(8)	!somme comp(i) O, O/Z
	  a(11,j)=-abon_rela(26)	!somme comp(i) Fe	  
	 ENDDO
	 b(1,3)=z0			!Z
	
	 a(4,4)=a(4,4)+1.d0	!C12	 		
	 a(4,5)=a(4,5)+1.d0	!C13
	
	 a(5,6)=a(5,6)+1.d0	!N14	 		
	 a(5,7)=a(5,7)+1.d0	!N15
	
	 a(6,8)=a(6,8)+1.d0	!O16	 		
	 a(6,9)=a(6,9)+1.d0	!O17
	 a(11,11)=a(11,11)+1.d0	!Fe	 
	
c	 rapports isotopiques	
	
c	 a(7,1)=1.d0		!H1		
	 a(7,2)=1.d0		!He3
	 a(7,3)=-he3she4z	!He3/He4 avec H2 dans He3		

	 a(8,5)=1.d0		!C13
	 a(8,4)=-c13sc12	!C13/C12
	
	 a(9,7)=1.d0		!N15
	 a(9,6)=-n15sn14	!N15/N14
	
	 a(10,9)=1.d0		!O17
	 a(10,8)=-o17so16	!O17/O16
	
c	 PRINT*,nchim
c	 DO i=1,nchim
c	  WRITE(*,2002)a(i,1:nchim),b(1,i)
c	 ENDDO

	 CALL gauss_band(a,b,indpc,nchim,nchim,nchim,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'ppcno9Fe, matrice du calcul des abond. non inversible'
	  PRINT*,'ARRET' ; STOP
	 ENDIF

c	 allocations diverses

	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),
	2 dmuex(nchim),dh2x(nchim),denx(nchim),dbe7x(nchim),dli7x(nchim))

c abondances initiales et abondances négligeables
	 
	 comp(1:nchim)=MAX(1.d-29,b(1,1:nchim))
	 ab_ini(1:nchim)=comp(1:nchim)*nucleo(1:nchim)
	 ab_min=ab_ini*fmin_abon	 
	 
c nombre/volume des métaux dans Z, indice de Fe56
		
	 nbz=SUM(comp(ihe4+1:nchim)) ; ife56=11	 
	 
c abondances en DeX, H=12

	 ALLOCATE(comp_dex(nchim))
	 comp_dex=12.d0+LOG10(comp/comp(1))
	 
c écritures

	 WRITE(2,2) ; WRITE(*,2) 
2	 FORMAT(/,'Réactions thermonucléaires des cycles PP, CNO',/)
	 WRITE(2,3)nreac ; WRITE(*,3)nreac 
3	 FORMAT('nombre de réactions : ',i3)
	 WRITE(2,4)nreac ; WRITE(*,4)nchim
4	 FORMAT('nombre d''éléments chimiques : ',i3)
	 WRITE(2,20)x0,y0,z0,z0/x0 ; WRITE(*,20)x0,y0,z0,z0/x0
20	 FORMAT(/,'abondances initiales/gramme déduites de:',/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,/,
	2 'H1=X0, H2+He3+He4=Y0, avec H2 dans He3',/,
	3 'Z0 = 1-X0-Y0 = C12+C13+N14+N15+O16+O17+Ex',/)	
	 WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
1	 FORMAT('H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	1 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	2 ', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,
	3 ', Ex:',es10.3,/,'Fe56:',es10.3)
	 WRITE(2,9)comp_dex ; WRITE(*,9)comp_dex
9	 FORMAT(/,'Abondances initiales en nombre: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	3 ', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,', Ex:',es10.3,/,
	4 'Fe56:',es10.3)
	 WRITE(2,21)(comp(4)+comp(5))/nbz,(comp(6)+comp(7))/nbz,
	1 (comp(8)+comp(9))/nbz,comp(10)/nbz
	 WRITE(*,21)(comp(4)+comp(5))/nbz,(comp(6)+comp(7))/nbz,
	1 (comp(8)+comp(9))/nbz,comp(10)/nbz,comp(11)/nbz
21	 FORMAT(/,'rapports en nombre dans Z:',/,'C/Z:',es10.3,', N/Z:',
	1 es10.3,', O/Z:',es10.3,', Ex/Z:',es10.3,', Fe/Z:',es10.3)
	 WRITE(2,14)he3she4z,c13sc12,n15sn14,o17so16
	 WRITE(*,14)he3she4z,c13sc12,n15sn14,o17so16
14	 FORMAT(/,'Rapports isotopiques en nombre:',/,
	1 'HE3/HE4=',es10.3,', C13/C12=',es10.3,
	2 ', N15/N14=',es10.3,', O17/O16=',es10.3)	
	 WRITE(2,5)ab_min(1:nchim) ; WRITE(*,5)ab_min(1:nchim)
5	 FORMAT(/,'abondances/gramme négligeables:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', C12:',es10.3,', C13:',es10.3,/,'N14:',es10.3,
	3 ', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,', Ex:',es10.3,/,
	4 'Fe56:',es10.3)
	 WRITE(2,6) ; WRITE(*,6)
6	 FORMAT(/,'H2, Li7, Be7 à l''équilibre')
	 WRITE(2,7) ; WRITE(*,7)
7	 FORMAT(/,'on utilise une table')
	 WRITE(2,8) ; WRITE(*,8)
8	 FORMAT(/,'évol. temporelle, test de précision sur H1 et He4')

c définitions diverses

	 DO i=1,nchim
	  ab_min(i)=ab_min(i)/nucleo(i)
	  anuc(i)=ANINT(nucleo(i))		!nombre atomique
	 ENDDO

c nettoyage

	 DEALLOCATE(a,b,comp_dex,indpc)
	 	 
c les réactions	 
	 
	CASE(2)
	 dcomp=0.d0 ; jac=0.d0

	 IF(t < t_inf)RETURN	
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)	

c	 WRITE(*,*)'comp' ; WRITE(*,2000)comp(1:nchim)
c	 WRITE(*,*)'réactions' ; WRITE(*,2000)r(1:nreac)

c	 équations d'évolution

	 dcomp(1)=-(3.*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+r(10)*comp(6)
	1 +(r(11)+r(12))*comp(7)+r(13)*comp(8)+r(14)*comp(9))*comp(1)
	2 +(2.d0*r(3)*comp(2)-r(4)*comp(3))*comp(2)		!H1
	 dcomp(2)=r(1)*comp(1)**2-(2.*r(3)*comp(2)
	1 +r(4)*comp(3))*comp(2)				!He3
	 dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
	1 +(r(12)*comp(7)+r(14)*comp(9))*comp(1)		!He4
	 dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)		!C12
	 dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)		!C13
	 dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9))*comp(1)!N14
	 dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7))*comp(1)	!N15
	 dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)		!O16
	 dcomp(9)=(r(13)*comp(8)-r(14)*comp(9))*comp(1)		!O17

c	   Pour vérifications SUM dcomp*nucleo=0

c	 PRINT*,'ppcno9Fe, vérifications SUM dcomp*nucleo=0'
c	 WRITE(*,2000)DOT_PRODUCT(dcomp,anuc) ; PAUSE'vérif' 

	 dcomp(10)=-DOT_PRODUCT(anuc,dcomp)/anuc(10)		!Ex	
	 
c	 calcul du jacobien

	 IF(deriv)THEN	!jac(i,j) : équation, i : élément j
	
c	  équation dcomp(1)
c	  dcomp(1)=-(3.*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)+
c	1 r(10)*comp(6)+(r(11)+r(12))*comp(7)+r(13)*comp(8)+
c	2 r(14)*comp(9))*comp(1)+(2.*r(3)*comp(2)
c	3 -r(4)*comp(3))*comp(2)	!H1

	  jac(1,1)=-6.*r(1)*comp(1)-r(8)*comp(4)-r(9)*comp(5)
	1 -r(10)*comp(6)-(r(11)+r(12))*comp(7)-r(13)*comp(8)
	2 -r(14)*comp(9)				!d /H1
	  jac(1,2)=4.*r(3)*comp(2)-r(4)*comp(3)		!d /He3
	  jac(1,3)=-r(4)*comp(2)			!d /He4
	  jac(1,4)=-r(8)*comp(1)			!d /C12
	  jac(1,5)=-r(9)*comp(1)			!d /C13
	  jac(1,6)=-r(10)*comp(1)			!d /N14
	  jac(1,7)=-(r(11)+r(12))*comp(1)		!d /N15
	  jac(1,8)=-r(13)*comp(1)			!d /O16
	  jac(1,9)=-r(14)*comp(1)			!d /O17
	 
	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(1,i)=jac(1,i) 
	1  -(3.*drx(1,i)*comp(1)+drx(8,i)*comp(4)
	2  +drx(9,i)*comp(5)+drx(10,i)*comp(6)
	3  +(drx(11,i)+drx(12,i))*comp(7)
	4  +drx(13,i)*comp(8)+drx(14,i)*comp(9))*comp(1)
	5  +(2.*drx(3,i)*comp(2)-drx(4,i)*comp(3))*comp(2)
	  ENDDO
			 
c	  équation dcomp(2)
c	  dcomp(2)=r(1)*comp(1)**2-(2.*r(3)*comp(2)
c	  +r(4)*comp(3))*comp(2)			!He3

	  jac(2,1)=2.*r(1)*comp(1)			!d /H1
	  jac(2,2)=-4.*r(3)*comp(2)-r(4)*comp(3)	!d /He3
	  jac(2,3)=-r(4)*comp(2)			!d /He4

	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(2,i)=jac(2,i)
	1  +drx(1,i)*comp(1)**2-(2.*drx(3,i)*comp(2)
	2  +drx(4,i)*comp(3))*comp(2)
	  ENDDO
	 
c	  équation dcomp(3)
c	  dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
c	1	+(r(12)*comp(7)+r(14)*comp(9))*comp(1)	!He4

	  jac(3,1)=r(12)*comp(7)+r(14)*comp(9)		!d /H1
	  jac(3,2)=2.*r(3)*comp(2)+r(4)*comp(3)		!d /He3
	  jac(3,3)=r(4)*comp(2)				!d /He4
	  jac(3,7)=r(12)*comp(1)			!d /N15
	  jac(3,9)=r(14)*comp(1)			!d /O17
	 
	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(3,i)=jac(3,i)
	1  +(drx(3,i)*comp(2)+drx(4,i)*comp(3))*comp(2)
	2  +(drx(12,i)*comp(7)+drx(14,i)*comp(9))*comp(1)
	  ENDDO
	 
c	  équation dcomp(4)
c	  dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)	!C12

	  jac(4,1)=-r(8)*comp(4)+r(12)*comp(7)			!d /H1
	  jac(4,4)=-r(8)*comp(1)				!d /C12
	  jac(4,7)=r(12)*comp(1)				!d /N15
	 
	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(4,i)=jac(4,i)
	1  +(-drx(8,i)*comp(4)+drx(12,i)*comp(7))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(5)
c	  dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)	!C13
	
	  jac(5,1)=r(8)*comp(4)-r(9)*comp(5)		!d /H1
	  jac(5,4)=r(8)*comp(1)				!d /C12
	  jac(5,5)=-r(9)*comp(1)				!d /C13

	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(5,i)=jac(5,i)+(drx(8,i)*comp(4)-drx(9,i)*comp(5))*comp(1)
	  ENDDO
	
c	  équation dcomp(6)
c	  dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9))*comp(1)!N14

	  jac(6,1)=r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9)	!d /H1
	  jac(6,5)=r(9)*comp(1)				!d /C13
	  jac(6,6)=-r(10)*comp(1)			!d /N14
	  jac(6,9)=r(14)*comp(1)			!d /O17

	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(6,i)=jac(6,i)
	1  +(drx(9,i)*comp(5)-drx(10,i)*comp(6)
	2  +drx(14,i)*comp(9))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(7)
c	  dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7))*comp(1)!N15
	
	  jac(7,1)=r(10)*comp(6)-(r(11)+r(12))*comp(7)	!d /H1
	  jac(7,6)=r(10)*comp(1)			            !d /N14
	  jac(7,7)=-(r(11)+r(12))*comp(1)		        !d /N15
	 
	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(7,i)=jac(7,i)
	1  +(drx(10,i)*comp(6)-(drx(11,i)+drx(12,i))*comp(7))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(8)
c	  dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)!O16

	  jac(8,1)=r(11)*comp(7)-r(13)*comp(8)		!d /H1
	  jac(8,7)=r(11)*comp(1)			!d /N15
	  jac(8,8)=-r(13)*comp(1)			!d /O16
	 
	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(8,i)=jac(8,i)+(drx(11,i)*comp(7)-drx(13,i)*comp(8))*comp(1)
	  ENDDO
	 
c	  équation dcomp(9)
c	  dcomp(9)=(r(13)*comp(8)-r(14)*comp(9))*comp(1)!O17

	  jac(9,1)=r(13)*comp(8)-r(14)*comp(9)		!d /H1
	  jac(9,8)=r(13)*comp(1)			!d /O16
	  jac(9,9)=-r(14)*comp(1)			!d /O17
	 
	  DO i=1,9		!dépendances dues à l'effet d'écran
	   jac(9,i)=jac(9,i)+(drx(13,i)*comp(8)-drx(14,i)*comp(9))*comp(1)
	  ENDDO		 

c équation dcomp(10)	 
c	  dcomp(10)=-SUM(anuc*dcomp)/anuc(10)!conservation des baryons

	  DO j=1,10
	   DO i=1,9
	    jac(10,j)=jac(10,j)+anuc(i)*jac(i,j)
	   ENDDO
	   jac(10,j)=-jac(10,j)/anuc(10)	   
	  ENDDO

c	  unités de temps pour intégration temporelle

	  jac=jac*secon6
	  
c	  PAUSE'après case2 deriv'
	 ENDIF		!deriv

	 dcomp=dcomp*secon6

	CASE(3)

c	 calcul de la production d'énergie nucléaire et dérivées
c	 pour H2(H,g)He3, q(2)H**2=q(2)*r(1)/r(2)
	 
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0
	 IF(t <= t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

c	 mue : nombre d'électrons / mole /g = 1/poids mol. moy. par e-

	 IF(comp(1) > 0.d0)THEN
	  h2=r(1)/r(2)*comp(1) ; den=r(6)*mue+r(7)*comp(1)
	  be7=r(4)*comp(2)*comp(3)/den ; li7=r(6)*be7*mue/r(5)/comp(1)
	 ELSE
	  h2=0.d0 ; be7=0.d0 ; li7=0.d0
	 ENDIF
	
c	 PRINT*,'h2,li7,be7' ; WRITE(*,2000)h2,li7,be7

	 epsilon(2)=(q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7)*comp(1)
	1 +(q(3)*comp(2)+q(4)*comp(3))*comp(2)+q(6)*mue*be7
	 epsilon(3)=(q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)
	1 +(q(11)+q(12))*comp(7)+q(13)*comp(8)+q(14)*comp(9))*comp(1)
	 DO i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 ENDDO
c	 PAUSE'apres case3 avant deriv'
	 
	 IF(deriv)THEN	
	  IF(h2 > 0.d0)THEN
	   dh2t=h2*(drt(1)/r(1)-drt(2)/r(2))
	   dh2ro=h2*(dro(1)/r(1)-dro(2)/r(2))
	   DO i=1,nchim
	    dh2x(i)=h2*(drx(1,i)/r(1)-drx(2,i)/r(2))
	   ENDDO
	   dh2x(1)=dh2x(1)+h2/comp(1)	  
	  ELSE
	   dh2t=0.d0 ; dh2ro=0.d0 ; dh2x=0.d0
	  ENDIF

	  IF(be7 > 0.d0)THEN
	   dent= drt(6)*mue+drt(7)*comp(1)
	   denro=dro(6)*mue+dro(7)*comp(1)
	   DO i=1,nchim-1	  
	    denx(i)=drx(6,i)*mue+r(6)*dmuex(i)+drx(7,i)*comp(1)
	   ENDDO
	   denx(1)=denx(1)+r(7)
	   	   
	   dbe7t= be7*(drt(4)/r(4)- dent/den)
	   dbe7ro=be7*(dro(4)/r(4)-denro/den)
	   DO i=1,nchim-1
	    dbe7x(i)=be7*(drx(4,i)/r(4)-denx(i)/den)
	   ENDDO
	   dbe7x(2)=dbe7x(2)+be7/comp(2)
	   dbe7x(3)=dbe7x(3)+be7/comp(3)	   
	  ELSE
	   dbe7t=0.d0 ; dbe7ro=0.d0 ; dbe7x=0.d0
	  ENDIF

	  IF(li7 > 0.d0)THEN
	   dli7t= li7*(drt(6)/r(6) +dbe7t/be7-drt(5)/r(5))
	   dli7ro=li7*(dro(6)/r(6)+dbe7ro/be7-dro(5)/r(5))
	   DO i=1,nchim-1	   
	    dli7x(i)=li7*(drx(6,i)/r(6)+dbe7x(i)/be7
	1	+dmuex(i)/mue-drx(5,i)/r(5))
	   ENDDO
	   dli7x(1)=dli7x(1)-li7/comp(1)	   
	  ELSE
	   dli7t=0.d0 ; dli7ro=0.d0 ; dli7x=0.d0
	  ENDIF
	  	
 	  et=(dqt(1)*comp(1)+dqt(2)*h2+dqt(5)*li7+dqt(7)*be7)*comp(1)
	1  +(dqt(3)*comp(2)+dqt(4)*comp(3))*comp(2)
	2  +dqt(6)*mue*be7+(dqt(8)*comp(4)+dqt(9)*comp(5)
	3  +dqt(10)*comp(6)+(dqt(11)+dqt(12))*comp(7)+dqt(13)*comp(8)
	4  +dqt(14)*comp(9))*comp(1)
	5  +(q(2)*dh2t+q(5)*dli7t+q(7)*dbe7t)*comp(1)+q(6)*mue*dbe7t
	
 	  ero=(dqo(1)*comp(1)+dqo(2)*h2+dqo(5)*li7+dqo(7)*be7)*comp(1)
	1  +(dqo(3)*comp(2)+dqo(4)*comp(3))*comp(2)
	2  +dqo(6)*mue*be7+(dqo(8)*comp(4)+dqo(9)*comp(5)
	3  +dqo(10)*comp(6)+(dqo(11)+dqo(12))*comp(7)+dqo(13)*comp(8)
	4  +dqo(14)*comp(9))*comp(1)
	5  +(q(2)*dh2ro+q(5)*dli7ro+q(7)*dbe7ro)*comp(1)+q(6)*mue*dbe7ro

	  ex(1)=2.d0*q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7
	1  +q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)+(q(11)+q(12))*comp(7)
	2  +q(13)*comp(8)+q(14)*comp(9)
	  ex(2)=2.d0*q(3)*comp(2)+q(4)*comp(3)
	  ex(3)=q(4)*comp(2) ; ex(4)=q(8)*comp(1) ; ex(5)=q(9)*comp(1)
	  ex(6)=q(10)*comp(1) ; ex(7)=(q(11)+q(12))*comp(1)
	  ex(8)=q(13)*comp(1) ; ex(9)=q(14)*comp(1)
	 
	  DO i=1,nchim-1	!contributions des écrans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*h2
	1   +q(2)*dh2x(i)+dqx(5,i)*li7+q(5)*dli7x(i)
	2   +dqx(7,i)*be7+q(7)*dbe7x(i))*comp(1)
	3   +(dqx(3,i)*comp(2)+dqx(4,i)*comp(3))*comp(2)
	4   +dqx(6,i)*mue*be7+q(6)*dmuex(i)*be7+q(6)*mue*dbe7x(i)
	5   +(dqx(8,i)*comp(4)+dqx(9,i)*comp(5)
	6   +dqx(10,i)*comp(6)+(dqx(11,i)+dqx(12,i))*comp(7)
	7   +dqx(13,i)*comp(8)+dqx(14,i)*comp(9))*comp(1)
	  ENDDO
c	  PAUSE'apres case3 deriv'
	 
	 ENDIF	!deriv
	 
	CASE(4)		!taux de production des neutrinos

	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)	 
	  be7=r(4)*comp(2)*comp(3)/(r(6)*mue+r(7)*comp(1))
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*be7/amu
	  b8e=r(7)*comp(1)*be7/amu ; n13e=r(8)*comp(1)*comp(4)/amu
	  o15e=r(10)*comp(1)*comp(6)/amu ; f17e=r(13)*comp(1)*comp(8)/amu
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0
	  n13e=0.d0 ; o15e=0.d0 ; f17e=0.d0
	 ENDIF
c	 PAUSE'apres case4'

	CASE DEFAULT
	 PRINT*,'ppcno0, fait ne peut prendre que les valeurs 1, 2, 3 ou 4'
	 PRINT*,'ERREUR: fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; STOP
	 
	END SELECT
	
	RETURN

	END SUBROUTINE ppcno9Fe
