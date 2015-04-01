
c**************************************************************************

	SUBROUTINE ppcno11(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c	routine private du module mod_nuc

c	cycles PP et CNO
c	cf. Clayton p. 380, 392 et 430,

c	éléments pris en compte:
c	H1, H2, He3, He4, Li7, C12, C13, N14, N15, O16, O17, Ex
c	Ex est l'élément fictif complément, il n'intéresse que la diffusion
c	Be7 est a l'équilibre


c	un premier appel a rq_reac initialise et définit le nb.
c	d'éléments chimiques pour lesquels les reac. nuc. sont tabulées
c	dans ppcno11 on ajoute Ex, soit nchim+1, puis

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	CESAM2k

c entree :
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
c	jac : jacobien (unité  de temps : 10**6 ans)
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
c	r(3) : réaction He3(He3,2H)He4
c	r(4) : réaction He4(He3,g)Be7
c	r(5) : réaction Li7(p,He4)He4
c	r(6) : réaction Be7(e-,nu g)Li7
c	r(7) : réaction Be7(p,g)B8(,e+ nu)Be8(,He4)He4

c	r(8) : réaction C12(p,g)N13(,e+ nu)C13		CNO
c	r(9) : réaction C13(p,g)N14
c	r(10) : réaction N14(p,g)O15(e+,nu)N15
c	r(11) : réaction N15(p,g)O16
c	r(12) : réaction N15(p,He4)C12
c	r(13) : réaction O16(p,g)F17(,e+ nu)O17
c	r(14) : réaction O17(p,He4)N14

c	indices des éléments
c	H1 : 1
c	H2 : 2
c	He3 : 3
c	He4 : 4
c	Li7 : 5
c	C12 : 6
c	C13 : 7
c	N14 : 8
c	N15 : 9
c	O16 : 10
c	O17 : 11
c	Ex : 12

c--------------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon, ihe4,
	1 ili7, i_ex, langue, nchim, nom_elem, nom_xheavy, nucleo,
	2 secon6, t_inf, x0, y0, zi, z0
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
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: anuc, dmuex,
	1 denx, dbe7x, drt, dro, r, q, dqt, dqo		
	REAL (kind=dp) :: mue, nbz, den, be7, dbe7he3, dbe7he4,
	1 dbe7mue, dbe7h, mass_ex, charge_ex, sum_a

	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc	
	INTEGER :: i, j
	
	LOGICAL :: inversible
	
	CHARACTER (len=2) :: text
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(5es15.8)
2002	FORMAT(11es8.1)

c	initialisations

	SELECT CASE(fait)
	CASE(0)
	 
c	 définition de nchim: nombre d'éléments chimiques dont on
c	 calcule l'abondance H1, H2, He3, He4, Li7, C13, C13, N14, N15,
c	 O16, O17, Ex

	 nchim=11+1 ; ili7=5

c	 appel d'initialisation pour tabulation des réactions nucléaires
c	 allocations fictives

	 ALLOCATE(drx(1,1),dqx(1,1),r(1),drt(1),dro(1),q(1),
	1 dqt(1),dqo(1),dmuex(1))
	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	 
	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 	 
	CASE(1)
	
c	 détermination des abondances initiales
c	 He3+He4=Y0
c	 Z0 = somme des éléments plus lourds que helium
c	 dans Z rapports en nombre
 
c	 Ex : élément fictif moyenne des éléments # Li + Be + CNO

	 CALL abon_ini

	 charge_ex=0.d0 ; mass_ex=0.d0 ; sum_a=0.d0
	 b1: DO i=4,nelem_ini		!a partir de Li+1=3+1
	  IF(elem(i) == ' C')CYCLE b1
	  IF(elem(i) == ' N')CYCLE b1
	  IF(elem(i) == ' O')CYCLE b1
	  charge_ex=charge_ex+c(i)*ab(i)	  	 
	  mass_ex=mass_ex+m(i)*ab(i)
	  sum_a=sum_a+ab(i)
	 ENDDO b1
	 charge_ex=nint(charge_ex/sum_a) ; mass_ex=nint(mass_ex/sum_a)	
	 WRITE(text,10)nint(mass_ex)
10	 FORMAT(i2)

c élément complément
	 i_ex=nchim		!indice de l'élément chimique reliquat
	 nucleo(nchim)=mass_ex	!nucleo de l'élément chimique reliquat
	 zi(nchim)=charge_ex	!charge de l'élément chimique reliquat
       i = nint(charge_ex)
	 nom_elem(nchim)=elem(i)//text !nom elem. chim. rel.
	 nom_xheavy=nom_elem(nchim)
 	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1023)TRIM(nom_elem(nchim)),nint(mass_ex),nint(charge_ex)
	  WRITE(2,1023)TRIM(nom_elem(nchim)),nint(mass_ex),nint(charge_ex)	 
1023	  FORMAT(a,': fictitious species /= CNO, of mass : ',i3,/,
	1 'and charge :',i3)	 
	 CASE DEFAULT	 
	  WRITE(*,23)TRIM(nom_elem(nchim)),nint(mass_ex),nint(charge_ex)
	  WRITE(2,23)TRIM(nom_elem(nchim)),nint(mass_ex),nint(charge_ex)	 
23	  FORMAT(a,': élément fictif /= CNO, de masse : ',i3,/,
	1 'et de charge :',i3)
	 END SELECT
c	 PRINT*,nchim,ihe4 ; WRITE(*,2000)nucleo(1:nchim)
c	 WRITE(*,2000)abon_rela
	 	 		
	 ALLOCATE(a(nchim,nchim),indpc(nchim),b(1,nchim))
	 a=0.d0 ; b=0.d0 ; indpc=1
	
	 a(1,1)=nucleo(1)  	!H1
	 a(1,2)=nucleo(2)  	!H2	
	 b(1,1)=x0
	
	 a(2,3)=nucleo(3)  	!He3
	 a(2,4)=nucleo(4)	!He4
	 b(1,2)=y0

	 DO j=5,12		!de Li7 a Ex
	  a(3,j)=nucleo(j)	!somme j >= 5 comp(j)*nucleo(j)=Z0
	  a(4,j)=-abon_rela(6)	!somme comp(i) C, C/Z
	  a(5,j)=-abon_rela(7)	!somme comp(i) N, N/Z	 
	  a(6,j)=-abon_rela(8)	!somme comp(i) O, O/Z
	  a(12,j)=-abon_rela(3)	!somme comp(i) Li, Li/Z	  
	 ENDDO		 		
	 b(1,3)=z0		!Z
	
	 a(4,6)=a(4,6)+1.d0	!C12	 		
	 a(4,7)=a(4,7)+1.d0	!C13
	
	 a(5,8)=a(5,8)+1.d0	!N14	 		
	 a(5,9)=a(5,9)+1.d0	!N15
	
	 a(6,10)=a(6,10)+1.d0	!O16	 		
	 a(6,11)=a(6,11)+1.d0	!O17
	 
	 a(12,5)=a(12,5)+1.d0	!Li7	 
	
c	 rapports isotopiques

	 a(7,1)=1.d0		!H1
	 a(7,2)=-1.d0/h2sh1	!H1/H2	
	
	 a(8,4)=1.d0		!He4
	 a(8,3)=-1.d0/he3she4	!He4/He3

	 a(9,7)=1.d0		!C13
	 a(9,6)=-c13sc12	!C13/C12
	
	 a(10,9)=1.d0		!N15
	 a(10,8)=-n15sn14	!N15/N14
	
	 a(11,11)=1.d0		!O17
	 a(11,10)=-o17so16	!O17/O16
	
c	 PRINT*,nchim
c	 DO i=1,nchim
c	  WRITE(*,2002)a(i,1:nchim),b(1,i)
c	 ENDDO

	 CALL gauss_band(a,b,indpc,nchim,nchim,nchim,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'ppcno11, matrice calcul des abondances non inversible'
	  PRINT*,'ARRET'
	  stop
	 ENDIF

c	 allocations diverses

	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),
	2 dmuex(nchim), denx(nchim),dbe7x(nchim))
	 
c	 abondances initiales et abondances negligeables
		 
	 comp(1:nchim)=max(1.d-29,b(1,1:nchim))
	 ab_ini(1:nchim)=comp(1:nchim)*nucleo(1:nchim)
	
c	 ab_min(1)=1.d-3	!H1
c	 ab_min(2)=1.d-20	!H2
c	 ab_min(3)=5.d-7	!He3
c	 ab_min(4)=1.d-3	!He4
c	 ab_min(5)=1.d-14	!Li7
c	 ab_min(6)=5.d-6	!C12
c	 ab_min(7)=1.d-7	!C13
c	 ab_min(8)=1.d-6	!N14
c	 ab_min(9)=5.d-9	!N15
c	 ab_min(10)=1.d-5	!O16
c	 ab_min(11)=5.d-9	!O17
c	 ab_min(12)=1.d-6	!Ex
	 
	 ab_min=ab_ini*fmin_abon
	 
c	 nombre/volume des métaux dans Z
		
	 nbz=sum(comp(ihe4+1:nchim))
		
	 WRITE(2,*)
	 WRITE(2,*)'Réactions thermonucléaires des cycles PP, CNO'
	 WRITE(2,*)
	 WRITE(2,"(' nombre de réactions : ',i3)")nreac
	 WRITE(2,"(' nombre d''éléments chimiques : ',i3)")nchim
	 WRITE(2,*) ; WRITE(2,20)x0,y0,z0,z0/x0	
20	 FORMAT('abondances initiales deduites de X0=',es10.3,
	1 ', Y0=',es10.3,', Z0=',es10.3,/,' Z0/X0=',es10.3,/,
	2 'H1+H2=X0, He3+He4=Y0',/,
	3 'Z0 = 1-X0-Y0 = Li7+C12+C13+N14+N15+O16+O17+Ex',/)		
	 WRITE(2,1)ab_ini(1:nchim)
1	 FORMAT(' H1:',es10.3,', H2:',es10.3,', He3:',es10.3,
	1 ', He4:',es10.3,', Li7:',es10.3,/,
	2 'C12:',es10.3,', C13:',es10.3,', N14:',es10.3,
	3 ', N15:',es10.3,', O16:',es10.3,/,'O17:',es10.3,', Ex:',es10.3)
	 WRITE(2,*)
	 WRITE(2,21)comp(5)/nbz,(comp(6)+comp(7))/nbz,	!C
	1 (comp(8)+comp(9))/nbz,	!N
	2 (comp(10)+comp(11))/nbz,	!O
	3 comp(12)/nbz	!other
21	 FORMAT(1x,'rapports en nombre, Li/Z:',es10.3,', C/Z:',es10.3,
	1	', N/Z:',es10.3,/,' O/Z:',es10.3,', Ex/Z:',es10.3,/)
	 WRITE(2,*)'abondances negligeables:' ; WRITE(2,1)ab_min(1:nchim)
	 WRITE(2,*) ; WRITE(2,*)'H2 est pris a l''équilibre'
	 WRITE(2,*) ; WRITE(2,*)'on utilise une table' ; WRITE(2,*)
	 WRITE(2,*)'évol. temporelle, test de précision sur H1 et He4'
	 WRITE(2,*)
	
	 PRINT* ; PRINT*,'Réactions thermonucléaires des cycles PP, CNO'
	 PRINT* ; WRITE(*,"(' nombre de réactions : ',i3)")nreac
	 WRITE(*,"(' nombre d''éléments chimiques : ',i3)")nchim
	 PRINT* ; WRITE(*,20)x0,y0,z0,z0/x0 ; WRITE(*,1)ab_ini(1:nchim)
	 PRINT*
	 WRITE(*,21)comp(5)/nbz,(comp(6)+comp(7))/nbz,	!C
	1 (comp(8)+comp(9))/nbz,		!N
	2 (comp(10)+comp(11))/nbz,		!O
	3 comp(12)/nbz		!other
	 PRINT*,'abondances negligeables:' ; WRITE(*,1)ab_min(1:nchim)
	 PRINT* ; PRINT*,'H2 est pris a l''équilibre' ; PRINT*
	 PRINT*,'on utilise une table' ; PRINT*
	 PRINT*,'évol. temporelle, test de précision sur H1 et He4'
	 PRINT*

	 DO i=1,nchim
	  ab_min(i)=ab_min(i)/nucleo(i)
	  anuc(i)=anint(nucleo(i))		!nombre atomique
	 ENDDO
	
	 DEALLOCATE(a,b,indpc)

c	 les réactions

	CASE(2)
	 dcomp=0.d0 ; jac=0.d0
	
	 IF(t < t_inf)return
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	
c	 PRINT*,'comp' ; WRITE(*,2000)(comp(i),i=1,nchim)
c	 PRINT*,'réactions' ; WRITE(*,2000)(r(i),i=1,nreac)

c	 Be7

	 den=(r(6)*mue+r(7)*comp(1)) ; be7=r(4)*comp(3)*comp(4)/den
	 dbe7he3=be7/comp(3) ; dbe7he4=be7/comp(4)
	 dbe7mue=-be7*r(6)/den ; dbe7h=-be7*r(7)/den
	
c	 WRITE(*,2000)be7,dbe7he3,dbe7he4,dbe7mue,dbe7h ; PAUSE

c	 équations d'evolution

	 dcomp(1)=-(2.d0*r(1)*comp(1)+r(2)*comp(2)+r(5)*comp(5)
	1 +r(7)*be7+r(8)*comp(6)+r(9)*comp(7)+r(10)*comp(8)
	2 +(r(11)+r(12))*comp(9)+r(13)*comp(10)
	3 +r(14)*comp(11))*comp(1)+2.d0*r(3)*comp(3)**2		!H1
	 dcomp(2)=(r(1)*comp(1)-r(2)*comp(2))*comp(1)		!H2
	 dcomp(3)=r(2)*comp(1)*comp(2)-(2.d0*r(3)*comp(3)
	1 +r(4)*comp(4))*comp(3)				!He3
	 dcomp(4)=(r(3)*comp(3)-r(4)*comp(4))*comp(3)
	1 +(2.d0*(r(5)*comp(5)+r(7)*be7)+r(12)*comp(9)
	2 +r(14)*comp(11))*comp(1)				!He4
	 dcomp(5)=-r(5)*comp(1)*comp(5)+r(6)*be7*mue		!Li7
	 dcomp(6)=(-r(8)*comp(6)+r(12)*comp(9))*comp(1)		!C12
	 dcomp(7)=(r(8)*comp(6)-r(9)*comp(7))*comp(1)		!C13
	 dcomp(8)=(r(9)*comp(7)-r(10)*comp(8)+r(14)*comp(11))*comp(1) !N14
	 dcomp(9)=(r(10)*comp(8)-(r(11)+r(12))*comp(9))*comp(1)	!N15
	 dcomp(10)=(r(11)*comp(9)-r(13)*comp(10))*comp(1)	!O16
	 dcomp(11)=(r(13)*comp(10)-r(14)*comp(11))*comp(1)	!O17


c	   Pour vérifications SUM dcomp*nucleo=0

c	 PRINT*,'ppcno11, vérifications SUM dcomp*nucleo=0'
c	 WRITE(*,2000)DOT_PRODUCT(dcomp,anuc) ; PAUSE'vérif' 

	 dcomp(12)=-DOT_PRODUCT(dcomp,anuc)/anuc(12) 	!cons. des baryons


c	 calcul du jacobien

	 IF(deriv)THEN	!jac(i,j) : équation, j : élément i
	
c	  équation
c	  dcomp(1)=-(2.d0*r(1)*comp(1)+r(2)*comp(2)+r(5)*comp(5)
c	1 +r(7)*be7+r(8)*comp(6)+r(9)*comp(7)+r(10)*comp(8)
c	2 +(r(11)+r(12))*comp(9)+r(13)*comp(10)
c	3 +r(14)*comp(11))*comp(1)+2.d0*r(3)*comp(3)**2		!H1

	  jac(1,1)=-4.d0*r(1)*comp(1)-r(2)*comp(2)-r(5)*comp(5)
	1 -r(7)*be7-r(8)*comp(6)-r(9)*comp(7)-r(10)*comp(8)
	2 -(r(11)+r(12))*comp(9)-r(13)*comp(10)-r(14)*comp(11)
	3 -r(7)*comp(1)*dbe7h					!d /H1
	  jac(1,2)=-r(2)*comp(1)				!d /H2
	  jac(1,3)=4.d0*r(3)*comp(3)-r(7)*comp(1)*dbe7he3	!d /He3
	  jac(1,4)=-r(7)*comp(1)*dbe7he4			!d /He3	 
	  jac(1,5)=-r(5)*comp(1)				!d /Li7
	  jac(1,6)=-r(8)*comp(1)				!d /C12
	  jac(1,7)=-r(9)*comp(1)				!d /C13
	  jac(1,8)=-r(10)*comp(1)			!d /N14
	  jac(1,9)=-(r(11)+r(12))*comp(1)		!d /N15
	  jac(1,10)=-r(13)*comp(1)			!d /O16
	  jac(1,11)=-r(14)*comp(1)			!d /O17
	 
	  DO i=1,nchim	!dependances dues a l'effet d'ecran et be7/muex
	   jac(1,i)=jac(1,i)
	1  -(2.d0*drx(1,i)*comp(1)+drx(2,i)*comp(2)
	2  +drx(5,i)*comp(5)+drx(7,i)*be7
	3  +drx(8,i)*comp(6)+drx(9,i)*comp(7)
	4  +drx(10,i)*comp(8)+(drx(11,i)
	5  +drx(12,i))*comp(9)+drx(13,i)*comp(10)
	6  +drx(14,i)*comp(11)+r(7)*dbe7mue*dmuex(i))*comp(1)
	7  +2.d0*drx(3,i)*comp(3)**2
	  ENDDO
	   	 
c	  équation dcomp(2)
c	  dcomp(2)=(r(1)*comp(1)-r(2)*comp(2))*comp(1)		!H2

	  jac(2,1)=2.d0*r(1)*comp(1)-r(2)*comp(2)		!d /H1
	  jac(2,2)=-r(2)*comp(1)				!d /H2
	 
	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(2,i)=jac(2,i)+(drx(1,i)*comp(1)-drx(2,i)*comp(2))*comp(1)
	  ENDDO
	 	 	 
c	  équation dcomp(3)
c	  dcomp(3)=r(2)*comp(1)*comp(2)-(2.d0*r(3)*comp(3)
c	1 +r(4)*comp(4))*comp(3)				!He3

	  jac(3,1)=r(2)*comp(2)				!d /H1
	  jac(3,2)=r(2)*comp(1)				!d /H2
	  jac(3,3)=-4.d0*r(3)*comp(3)-r(4)*comp(4)		!d /He3
	  jac(3,4)=-r(4)*comp(3)				!d /He4
	 
	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(3,i)=jac(3,i)
	1  +drx(2,i)*comp(1)*comp(2)-(2.d0*drx(3,i)*comp(3)
	2  +drx(4,i)*comp(4))*comp(3)
	  ENDDO

c	  équation dcomp(4)
c	  dcomp(4)=(r(3)*comp(3)-r(4)*comp(4))*comp(3)
c	1 +(2.d0*(r(5)*comp(5)+r(7)*be7)+r(12)*comp(9)
c	2 +r(14)*comp(11))*comp(1)				!He4

	  jac(4,1)=2.d0*(r(5)*comp(5)+r(7)*be7+r(7)*dbe7h*comp(1))
	1 +r(12)*comp(9)+r(14)*comp(11)				!d /H1
	  jac(4,3)=2.d0*r(3)*comp(3)-r(4)*comp(4)+2.d0*r(7)*dbe7he3 !d /He3
	  jac(4,4)=-r(4)*comp(3)+2.d0*r(7)*dbe7he4		!d /He4
	  jac(4,5)=r(5)*comp(1)*2.d0				!d /Li7
	  jac(4,9)=r(12)*comp(1)				!d /N15
	  jac(4,11)=r(14)*comp(1)				!d /O17
	 
	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(4,i)=jac(4,i)
	1  +(drx(3,i)*comp(3)-drx(4,i)*comp(4))*comp(3)
	2  +(2.d0*(drx(5,i)*comp(5)+drx(7,i)*be7)
	3  +drx(12,i)*comp(9)+2.d0*r(7)*dbe7mue*dmuex(i)
	4  +drx(14,i)*comp(11))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(5)
c	  dcomp(5)=-r(5)*comp(1)*comp(5)+r(6)*be7*mue	!Li7

	  jac(5,1)=-r(5)*comp(5)+r(6)*dbe7h*mue		!d /H1
	  jac(5,3)=r(6)*dbe7he3*mue			!d /He3	 
	  jac(5,4)=r(6)*dbe7he4*mue			!d /He4	 
	  jac(5,5)=-r(5)*comp(1)			!d /Li7
	 
	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(5,i)=jac(5,i)
	1	-drx(5,i)*comp(1)*comp(5)
	2	+drx(6,i)*be7*mue+r(6)*(dbe7mue*mue+be7)*dmuex(i)
	  ENDDO
	 	 	 
c	  équation dcomp(6)	 
c	  dcomp(6)=(-r(8)*comp(6)+r(12)*comp(9))*comp(1)	!C12

	  jac(6,1)=-r(8)*comp(6)+r(12)*comp(9)		!d /H1
	  jac(6,6)=-r(8)*comp(1)				!d /C12
	  jac(6,9)=r(12)*comp(1)				!d /N15
	 
	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(6,i)=jac(6,i)
	1	+(-drx(8,i)*comp(6)+drx(12,i)*comp(9))*comp(1)
	  ENDDO	 
	 	 
c	  équation dcomp(7)
c	  dcomp(7)=(r(8)*comp(6)-r(9)*comp(7))*comp(1)		!C13

	  jac(7,1)=r(8)*comp(6)-r(9)*comp(7)			!d /H1
	  jac(7,6)=r(8)*comp(1)				!d /C12
	  jac(7,7)=-r(9)*comp(1)				!d /C13

	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(7,i)=jac(7,i)
	1  +(drx(8,i)*comp(6)-drx(9,i)*comp(7))*comp(1)
	  ENDDO
	 
c	  équation dcomp(8)	!N14	 
c	  dcomp(8)=(r(9)*comp(7)-r(10)*comp(8)+r(14)*comp(11))*comp(1)

	  jac(8,1)=r(9)*comp(7)-r(10)*comp(8)+r(14)*comp(11)	!d /H1
	  jac(8,7)=r(9)*comp(1)				!d /C13
	  jac(8,8)=-r(10)*comp(1)				!d /N14
	  jac(8,11)=r(14)*comp(1)				!d /O17

	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(8,i)=jac(8,i)
	1  +(drx(9,i)*comp(7)-drx(10,i)*comp(8)
	2  +drx(14,i)*comp(11))*comp(1)	  
	  ENDDO
	 
c	  équation dcomp(9)	 
c	  dcomp(9)=(r(10)*comp(8)-(r(11)+r(12))*comp(9))*comp(1) !N15

	  jac(9,1)=r(10)*comp(8)-(r(11)+r(12))*comp(9)	!d /H1
	  jac(9,8)=r(10)*comp(1)				!d /N14
	  jac(9,9)=-(r(11)+r(12))*comp(1)			!d /N15

	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(9,i)=jac(9,i)
	1  +(drx(10,i)*comp(8)-(drx(11,i)
	2  +drx(12,i))*comp(9))*comp(1)
	  ENDDO	  
	 
c	  équation dcomp(10)	 
c	  dcomp(10)=(r(11)*comp(9)-r(13)*comp(10))*comp(1)	!O16

	  jac(10,1)=r(11)*comp(9)-r(13)*comp(10)		!d /H1
	  jac(10,9)=r(11)*comp(1)				!d /N15
	  jac(10,10)=-r(13)*comp(1)				!d /O16
	 
	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(10,i)=jac(10,i)	 
	1  +(drx(11,i)*comp(9)-drx(13,i)*comp(10))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(11)
c	  dcomp(11)=(r(13)*comp(10)-r(14)*comp(11))*comp(1)	!O17

	  jac(11,1)=r(13)*comp(10)-r(14)*comp(11)		!d /H1
	  jac(11,10)=r(13)*comp(1)				!d /O16
	  jac(11,11)=-r(14)*comp(1)				!d /O17

	  DO i=1,nchim		!dependances dues a l'effet d'ecran
	   jac(11,i)=jac(11,i)
	1  +(drx(13,i)*comp(10)-drx(14,i)*comp(11))*comp(1)	  
	  ENDDO		 
	 
c	  DO i=1,11					!Ex
c	   dcomp(12)=dcomp(12)+anuc(i)*dcomp(i)
c	  ENDDO
c	  dcomp(12)=-dcomp(12)/anuc(12) !conservation des baryons
	  DO j=1,12
	   DO i=1,11
	    jac(12,j)=jac(12,j)+anuc(i)*jac(i,j)
	   ENDDO
	   jac(12,j)=-jac(12,j)/anuc(12)
	  ENDDO

c	  unité s de temps pour integration temporelle

	  jac=jac*secon6
	 
	 ENDIF
	
	 dcomp=dcomp*secon6

c	 calcul de la production d'énergie nucléaire et dérivées

	CASE(3)
	
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0
	 IF(t <= t_inf)return
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	 
c	 mue : nombre d'électrons / mole /g = 1/poids mol. moy. par e-

c	 Be7

	 den=(r(6)*mue+r(7)*comp(1)) ; be7=r(4)*comp(3)*comp(4)/den
	 dbe7he3=be7/comp(3); dbe7he4=be7/comp(4)
	 dbe7mue=-r(6)*be7/den ; dbe7h=-r(7)*be7/den

	 epsilon(2)=(q(1)*comp(1)+q(2)*comp(2)+q(5)*comp(5)+q(7)*be7)
	1 *comp(1)+(q(3)*comp(3)+q(4)*comp(4))*comp(3)+q(6)*mue*be7
	 epsilon(3)=(q(8)*comp(6)+q(9)*comp(7)+q(10)*comp(8)+
	1 (q(11)+q(12))*comp(9)+q(13)*comp(10)+q(14)*comp(11))*comp(1)
	 DO i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 ENDDO

	 IF(deriv)THEN	
	  et=(dqt(1)*comp(1)+dqt(2)*comp(2)+dqt(5)*comp(5)+dqt(7)*be7)
	1 *comp(1)+(dqt(3)*comp(3)+dqt(4)*comp(4))*comp(3)
	2 +dqt(6)*mue*be7+(dqt(8)*comp(6)+dqt(9)*comp(7)
	3 +dqt(10)*comp(8)+(dqt(11)+dqt(12))*comp(9)+dqt(13)*comp(10)
	4 +dqt(14)*comp(11))*comp(1)
		
	  ero=(dqo(1)*comp(1)+dqo(2)*comp(2)+dqo(5)*comp(5)+dqo(7)*be7)
	1 *comp(1)+(dqo(3)*comp(3)+dqo(4)*comp(4))*comp(3)
	2 +dqo(6)*mue*be7+(dqo(8)*comp(6)+dqo(9)*comp(7)
	3 +dqo(10)*comp(8)+(dqo(11)+dqo(12))*comp(9)+dqo(13)*comp(10)
	4 +dqo(14)*comp(11))*comp(1)
		
	  ex(1)=2.d0*q(1)*comp(1)+q(2)*comp(2)+q(5)*comp(5)
	1 +q(8)*comp(6)+q(9)*comp(7)+q(10)*comp(8)
	2 +(q(11)+q(12))*comp(9)+q(13)*comp(10)+q(14)*comp(11)
	3 +q(7)*(be7+dbe7h*comp(1))+q(6)*mue*dbe7h
	  ex(2)=q(2)*comp(1)
	  ex(3)=2.d0*q(3)*comp(3)+q(4)*comp(4)
	1 +q(7)*dbe7he3*comp(1)+q(6)*mue*dbe7he3
	  ex(4)=q(4)*comp(3)+q(7)*dbe7he4*comp(1)+q(6)*mue*dbe7he4
	  ex(5)=q(5)*comp(1)
	  ex(6)=q(8)*comp(1)
	  ex(7)=q(9)*comp(1)
	  ex(8)=q(10)*comp(1)
	  ex(9)=(q(11)+q(12))*comp(1)
	  ex(10)=q(13)*comp(1)
	  ex(11)=q(14)*comp(1)
	 
	  DO i=1,nchim	!contributions des ecrans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*comp(2)
	1 +dqx(5,i)*comp(5)+dqx(7,i)*be7)*comp(1)
	2 +(dqx(3,i)*comp(3)+dqx(4,i)*comp(4))*comp(3)
	3 +dqx(6,i)*mue*be7+(dqx(8,i)*comp(6)
	4 +dqx(9,i)*comp(7)+dqx(10,i)*comp(8)
	5 +(dqx(11,i)+dqx(12,i))*comp(9)
	6 +dqx(13,i)*comp(10)+dqx(14,i)*comp(11))*comp(1)
	7 +(q(7)*dbe7mue*comp(1)+q(6)*(dbe7mue*mue+be7))*dmuex(i)
	  ENDDO	 
	 ENDIF	!deriv

c	 production de neutrinos

	CASE(4)
	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	  den=(r(6)*mue+r(7)*comp(1)) ; be7=r(4)*comp(3)*comp(4)/den	
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*be7/amu
	  b8e=r(7)*comp(1)*be7/amu ; n13e=r(8)*comp(1)*comp(6)/amu
	  o15e=r(10)*comp(1)*comp(8)/amu ; f17e=r(13)*comp(1)*comp(10)/amu
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0 ; n13e=0.d0 
	  o15e=0.d0 ; f17e=0.d0
	 ENDIF
	 
	CASE DEFAULT
	 PRINT*,'ppcno11, fait ne peut valoir que 1, 2, 3 ou 4'
	 PRINT*,'ERREUR fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; stop
	
	END SELECT
	
	RETURN

	END SUBROUTINE ppcno11
