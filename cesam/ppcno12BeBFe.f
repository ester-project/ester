
c**************************************************************************

	SUBROUTINE ppcno12BeBFe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c	routine private du module mod_nuc

c	cycles PP et CNO
c	cf. Clayton p. 380, 392 et 430, 

c	éléments pris en compte:
c	H1, H2, He3, He4, Li6, Li7, Be7, Be9, B11, C12, C13, N14, N15, O16, O17,
c	Fe56, Ex
c	Ex est l'élément fictif complément, il n'intéresse que la diffusion
c	aucun élément a l'équilibre

c	un premier appel à rq_reac-->tabul_nuc initialise et définit le nb.
c	d'isotopes pour lesquels les réac. nuc. sont tabulées
c	dans ppcno12BeBFe

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrée :
c	t:température cgs
c	ro:densite cgs
c	comp:abondances
c	deriv=.true.:on calcule le jacobien
c	fait=1:initialisation de la composition chimique
c	    =2:calcul de dcomp et jacobien si deriv
c	    =3:énergie nucléaire et dérivées / t et ro
c	    =4:production de neutrinos

c sorties
c	dcomp: dérivée temporelle (unité de temps:10**6 ans)
c	jac: jacobien (unité de temps:10**6 ans)
c	epsilon, et, ero, ex: énergie thermonucléaire (unité de temps:s)
c			  :et dérivées /t, ro ,X

c	hhe, be7e, b8e, n13e, o15e, f17e:nombre de neutrinos g/s
c	hhe réaction:H1(p,e+ nu)H2
c	be7e réaction:Be7(e-,nu g)Li7
c	b8e réaction:B8(,e+ nu)Be8
c	n13e réaction:N13(,e+ nu)C13
c	o15e réaction:O15(e+,nu)N15 
c	f17e réaction:F17(,e+ nu)O17

c initialisation
c	ab_min:abondances negligeables
c	ab_ini:abondances initiales

c	r(1) : réaction H1(p,e+ nu)H2			PP
c	r(2) : réaction H2(p,g)He3
c	r(3) : réaction He3(He3,2p)He4
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

c	r(15) : réaction Be9(p,d)2He4	autres réactions
c	r(16) : réaction Li6(p,He3)He4
c	r(17) : réaction Li6(p,g)Be7
c	r(18) : réaction Be9(p,a)Li6
c	r(19) : réaction B11(p,a)2He4
c	r(20) : réaction B11(p,g)C12

c	indices des éléments
c	H1 : 1
c	H2 : 2
c	He3 : 3
c	He4 : 4
c	Li7 : 5
c	Be7 : 6
c	C12 : 7
c	C13 : 8
c	N14 : 9
c	N15 : 10
c	O16 : 11
c	O17 : 12
c	Be9 : 13
c	Ex : 14
c	B11 : 15
c	Fe : 16
c	Li6 : 17

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, ife56, ihe4,
	1 langue, nchim, nom_elem, nom_xheavy, nucleo,
	2 rot_solid, secon6, t_inf, x0, y0, zi, z0
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
	1 dmuex, drt, dro, r, q, dqt, dqo		
	REAL (kind=dp) :: mue, nbz, mass_ex, charge_ex, sum_a

	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc		
	INTEGER :: i, j
	
	LOGICAL :: inversible
	
	CHARACTER (len=2) :: text
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(5es15.8)
2002	FORMAT(12es8.1)

c	initialisations

	SELECT CASE(fait)
	CASE(0)

c	 définition de nchim: nombre d'éléments chimiques dont on
c	 calcule l'abondance H1, H2, He3, He4, Li7, Be7, C12, C13,
c	 N14, N15, O16, O17, Be9, Ex, B11, Fe56, Li6

	 nchim=17

c	 appel d'initialisation pour tabulation des réactions nucléaires
c	 allocations fictives

	 ALLOCATE(drx(0,0),dqx(0,0),r(0),drt(0),dro(0),q(0),
	1 dqt(0),dqo(0),dmuex(0))

	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 
	CASE(1)
	
c	 détermination des abondances initiales
c	 He3+He4=Y0
c	 Z0 = somme des éléments plus lourds que hélium
c	 dans Z rapports en nombre

	 CALL abon_ini	 
	 
c	 Ex : élément fictif moyenne des éléments # Li + Be + CNO + B + Fe
	 
	 charge_ex=0.d0 ; mass_ex=0.d0 ; sum_a=0.d0
	 b1: DO i=3,nelem_ini		!à partir de Li
	  IF(elem(i) == 'Li')CYCLE B1
	  IF(elem(i) == 'Be')CYCLE B1
	  IF(elem(i) == ' B')CYCLE B1	 
	  IF(elem(i) == ' C')CYCLE B1
	  IF(elem(i) == ' N')CYCLE B1
	  IF(elem(i) == ' O')CYCLE B1
	  IF(elem(i) == 'Fe')CYCLE B1	  
	  charge_ex=charge_ex+c(i)*ab(i)	  	 
	  mass_ex=mass_ex+m(i)*ab(i) ; sum_a=sum_a+ab(i)
	 ENDDO B1
	 charge_ex=nint(charge_ex/sum_a) ; mass_ex=mass_ex/sum_a
	 WRITE(text,10)nint(mass_ex)
10	 FORMAT(i2)

c élément fictif
	 nucleo(14)=mass_ex	!nucleo de l'élément chimique reliquat
	 zi(14)=charge_ex	!charge de l'élément chimique reliquat
       i = nint(charge_ex)
	 nom_elem(14)=elem(i)//text !nom elem. chim. rel
 	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(2,1011)nom_elem(14) ; WRITE(*,1011)nom_elem(14)
1011	  FORMAT('fictitious element /= H, He, Li, Be, B, CNO, Fe : ',a)
	  WRITE(2,1012)nint(mass_ex) ; WRITE(*,1012)nint(mass_ex)
1012	  FORMAT('mass of the fictitious element : ',i3)	
	  WRITE(2,1013)nint(charge_ex) ; WRITE(*,1013)nint(charge_ex)
1013	  FORMAT('charge of the fictitious element : ',i3)	 
	 CASE DEFAULT	 
	  WRITE(2,11)nom_elem(14) ; WRITE(*,11)nom_elem(14)
11	  FORMAT('Ex, elem. fictif /= H, He, Li, Be, B, CNO, Fe : ',a)
	  WRITE(2,12)nint(mass_ex) ; WRITE(*,12)nint(mass_ex)
12	  FORMAT('masse de l''élément fictif : ',i3)	
	  WRITE(2,13)nint(charge_ex) ; WRITE(*,13)nint(charge_ex)
13	  FORMAT('charge de l''élément fictif : ',i3)	
	 END SELECT	 	 	 
c	 PRINT*,nchim ; WRITE(*,2000)nucleo ; PAUSE'début'
	 
c détermination des abondances initiales, a(équation,élément) 
	 ALLOCATE(a(nchim,nchim),indpc(nchim),b(1,nchim))
	 a=0.d0 ; b=0.d0 ; indpc=1
	
	 a(1,1)=nucleo(1)  	!H1  abondance de H 
	 a(1,2)=nucleo(2)  	!H2	
	 b(1,1)=x0
	
	 a(2,3)=nucleo(3)  	!He3  abondance de He
	 a(2,4)=nucleo(4)	!He4
	 b(1,2)=y0

	 DO j=5,nchim		!à partir de Li7 ie. sauf H et He
	  a(3,j)=nucleo(j)	!somme i > 5 comp(i)*nucleo(i)=Z0
	  a(4,j)=-abon_rela(6)	!somme comp(i) C, C/Z
	  a(5,j)=-abon_rela(7)	!somme comp(i) N, N/Z	 
	  a(6,j)=-abon_rela(8)	!somme comp(i) O, O/Z
	  a(12,j)=-abon_rela(3)	!somme comp(i) Li, Li/Z	  
	  a(14,j)=-abon_rela(4)	!somme comp(i) Be, Be/Z
	  a(15,j)=-abon_rela(5)	!somme comp(i) B, B/Z
	  a(16,j)=-abon_rela(26)!somme comp(i) Fe, Fe/Z	  
	 ENDDO		 		
	 b(1,3)=z0		!Z
	
	 a(4,7)=a(4,7)+1.d0	!C12	 		
	 a(4,8)=a(4,8)+1.d0	!C13
	
	 a(5,9)=a(5,9)+1.d0	!N14	 		
	 a(5,10)=a(5,10)+1.d0	!N15
	
	 a(6,11)=a(6,11)+1.d0	!O16	 		
	 a(6,12)=a(6,12)+1.d0	!O17
	 
	 a(12,5)=a(12,5)+1.d0	!Li7
	 a(12,17)=a(12,17)+1.d0	!Li6
	 	 
	 a(14,6)=a(14,6)+1.d0	!Be7
	 a(14,13)=a(14,13)+1.d0	!Be9
	 	 
	 a(15,15)=a(15,15)+1.d0	!B11
	 
	 a(16,16)=a(16,16)+1.d0	!Fe	 	 	 
	
c	 rapports isotopiques

	 a(7,1)=1.d0		!H1
	 a(7,2)=-1./h2sh1	!H1/H2	
	
	 a(8,4)=1.d0		!He4
	 a(8,3)=-1./he3she4	!He4/He3

	 a(9,8)=1.d0		!C13
	 a(9,7)=-c13sc12	!C13/C12
	
	 a(10,10)=1.d0		!N15
	 a(10,9)=-n15sn14	!N15/N14
	
	 a(11,12)=1.d0		!O17
	 a(11,11)=-o17so16	!O17/O16
	
	 a(13,6)=1.d0		!Be7
	 a(13,13)=-be7sbe9	!Be7/Be9
	 
	 a(17,17)=1.d0		!Li6
	 a(17,5)=-li6sli7	!Li6/Li7
	
c	 PRINT*,nchim
c	 DO i=1,nchim
c	  WRITE(*,2002)a(i,1:nchim),comp(i)
c	 ENDDO
	
	 CALL gauss_band(a,b,indpc,nchim,nchim,nchim,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'ppcno12BeBFe, matrice du calcul des abon. non inversible'
	  PRINT*,'ARRET'
	  STOP
	 ENDIF

c allocations diverses

	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),dmuex(nchim))
	 
c abondances initiales et abondances négligeables
		 
	 comp(1:nchim)=max(1.d-29,b(1,1:nchim))
	 ab_ini=comp*nucleo		 
	 ab_min=ab_ini*1.d-2	 

c nombre/volume des métaux dans Z, indice de Fe56
		
	 nbz=SUM(comp(ihe4+1:nchim)) ; ife56=16
	 
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
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,
	2 ', H1+H2=X0, He3+He4=Y0',/,
	3 'Z0 = 1-X0-Y0 = Li7+Be7+Be9+B11+C12+C13+N14+N15+O16+O17+Ex+Fe',/)
	 WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
1	 FORMAT('H1:',es10.3,', H2:',es10.3,', He3:',
	1 es10.3,', He4:',es10.3,', Li7:',es10.3,/,
	2 'Be7:',es10.3,', C12:',es10.3,', C13:',es10.3,
	3 ', N14:',es10.3,', N15:',es10.3,/,'O16:',es10.3,', O17:',es10.3,
	4 ', Be9:',es10.3,', Ex:',es10.3,', B11:',es10.3,/,'Fe56:',
	5 es10.3,', Li6:',es10.3)
	 WRITE(2,9)comp_dex ; WRITE(*,9)comp_dex
9	 FORMAT(/,'Abondances initiales en nombre: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', H2:',es10.3,', He3:',
	2 es10.3,', He4:',es10.3,', Li7:',es10.3,/,
	3 'Be7:',es10.3,', C12:',es10.3,', C13:',es10.3,
	4 ', N14:',es10.3,', N15:',es10.3,/,'O16:',es10.3,', O17:',es10.3,
	5 ', Be9:',es10.3,', Ex:',es10.3,', B11:',es10.3,/,'Fe56:',
	6 es10.3,', Li6:',es10.3)	
	 WRITE(2,21)(comp(5)+comp(17))/nbz,	!Li/Z
	1 (comp(6)+comp(13))/nbz,		!Be/Z
	2 comp(15)/nbz,				!B/Z	
	3 (comp(7)+comp(8))/nbz,		!C/Z
	4 (comp(9)+comp(10))/nbz,		!N/Z
	5 (comp(11)+comp(12))/nbz,		!O/Z
	6 comp(14)/nbz,comp(16)/nbz		!Ex/Z, Fe/Z
	 WRITE(*,21)(comp(5)+comp(17))/nbz,	!Li/Z
	1 (comp(6)+comp(13))/nbz,		!Be/Z
	2 comp(15)/nbz,				!B/Z	
	3 (comp(7)+comp(8))/nbz,		!C/Z
	4 (comp(9)+comp(10))/nbz,		!N/Z
	5 (comp(11)+comp(12))/nbz,		!O/Z
	6 comp(14)/nbz,comp(16)/nbz		!Ex/Z, Fe/Z
21	 FORMAT(/,'rapports en nombre dans Z, Li/Z:',es10.3,
	1 ', Be/Z:',es10.3, ', B/Z:',es10.3,/,
	1 'C/Z:',es10.3,', N/Z:',es10.3,', O/Z:',es10.3,', Ex/Z:',es10.3,
	2 ', Fe/Z:',es10.3)
	 WRITE(2,14)h2sh1,he3she4,c13sc12,n15sn14,o17so16,be7sbe9,li6sli7
	 WRITE(*,14)h2sh1,he3she4,c13sc12,n15sn14,o17so16,be7sbe9,li6sli7	  
14	 FORMAT(/,'Rapports isotopiques en nombre:',/,
	1 'H2/H1=',es10.3,', HE3/HE4=',es10.3,', C13/C12=',es10.3,
	2 ', N15/N14=',es10.3,/,'O17/O16=',es10.3,', BE7/BE9=',es10.3,
	3 ', LI6/LI7=',es10.3)			
	 WRITE(2,5)ab_min(1:nchim) ; WRITE(*,5)ab_min(1:nchim)
5	 FORMAT(/,'abondances/gramme négligeables:',/,
	1'H1:',es10.3,', H2:',es10.3,', He3:',
	2 es10.3,', He4:',es10.3,', Li7:',es10.3,/,
	3 'Be7:',es10.3,', C12:',es10.3,', C13:',es10.3,
	4 ', N14:',es10.3,', N15:',es10.3,/,'O16:',es10.3,', O17:',es10.3,
	5 ', Be9:',es10.3,', Ex:',es10.3,', B11:',es10.3,/,'Fe56:',
	6 es10.3,', Li6:',es10.3)
	 WRITE(2,6) ; WRITE(*,6)
6	 FORMAT(/,'aucun élément à l''équilibre')
	 WRITE(2,7) ; WRITE(*,7)
7	 FORMAT(/,'on utilise une table')	 
	 WRITE(2,8) ; WRITE(*,8)
8	 FORMAT(/,'évol. temporelle, test de précision sur H1 et He4')

	 ab_min=ab_min/nucleo ; anuc=anint(nucleo)	!nombre atomique
	
c	 stor=0.d0
c	 DO i=1,4
c	  stor=stor+ab_ini(i)
c	 ENDDO
c	 sum=0.d0
c	 DO i=1,nchim
c	  sum=sum+ab_ini(i)
c	 ENDDO
c	 WRITE(*,2000)(comp(i),i=1,nchim); PRINT*
c	 WRITE(*,2000)(ab_ini(i),i=1,nchim) ; PRINT*
c	 WRITE(*,2001)stor,1.d0-stor,z0,sum ; PAUSE

c nettoyage

	 DEALLOCATE(a,b,comp_dex,indpc)
	
c les réactions

	CASE(2)
	 dcomp=0.d0 ; jac=0.d0

	 IF(t < t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	
c	 PRINT*,'comp' ; WRITE(*,2000)comp(1:nchim)
c	 PRINT*,'réactions' ; WRITE(*,2000)r(1:nreac)

c	 équations d'évolution

	 dcomp(1)=-(2.d0*r(1)*comp(1)+r(2)*comp(2)+r(5)*comp(5)
	1 +r(7)*comp(6)+r(8)*comp(7)+r(9)*comp(8)+r(10)*comp(9)
	2 +(r(11)+r(12))*comp(10)+r(13)*comp(11)
	3 +r(14)*comp(12)+(r(15)+r(18))*comp(13)+(r(16)+r(17))*comp(17)
	4 +(r(19)+r(20))*comp(15))*comp(1)
	5 +2.d0*r(3)*comp(3)**2					!H1

	 dcomp(2)=(r(1)*comp(1)-r(2)*comp(2))*comp(1)	
	1 +r(15)*comp(1)*comp(13)				!H2

	 dcomp(3)=r(2)*comp(1)*comp(2)-(2.d0*r(3)*comp(3)
	1 +r(4)*comp(4))*comp(3)+r(16)*comp(17)*comp(1)		!He3

	 dcomp(4)=(r(3)*comp(3)-r(4)*comp(4))*comp(3)
	1 +(2.d0*(r(5)*comp(5)+r(7)*comp(6))+r(12)*comp(10)
	2 +r(14)*comp(12)+2.d0*r(15)*comp(13)
	3 +r(16)*comp(17)+3.d0*r(19)*comp(15))*comp(1)		!He4

	 dcomp(5)=-r(5)*comp(1)*comp(5)+r(6)*comp(6)*mue	!Li7

	 dcomp(6)=r(4)*comp(3)*comp(4)-(r(6)*mue+r(7)*comp(1))*comp(6)
	1 +r(17)*comp(17)*comp(1) 				!Be7

	 dcomp(7)=(-r(8)*comp(7)+r(12)*comp(10)
	1 +r(20)*comp(15))*comp(1)				!C12

	 dcomp(8)=(r(8)*comp(7)-r(9)*comp(8))*comp(1)		!C13

	 dcomp(9)=(r(9)*comp(8)-r(10)*comp(9)+r(14)*comp(12))*comp(1) !N14

	 dcomp(10)=(r(10)*comp(9)-(r(11)+r(12))*comp(10))*comp(1) !N15

	 dcomp(11)=(r(11)*comp(10)-r(13)*comp(11))*comp(1)	!O16

	 dcomp(12)=(r(13)*comp(11)-r(14)*comp(12))*comp(1)	!O17

	 dcomp(13)=-(r(15)+r(18))*comp(1)*comp(13)		!Be9

	 dcomp(14)=0.d0 					!Ex

	 dcomp(15)=-(r(19)+r(20))*comp(1)*comp(15)		!B11
	 
c	 dcomp(16)=0.d0 	!Fe56	 

	 dcomp(17)=(-(r(16)+r(17))*comp(17)+r(18)*comp(13))*comp(1)  !Li6
	 
c	   Pour vérifications SUM dcomp*nucleo=0

c	 PRINT*,'ppcno12BeBFe, vérifications SUM dcomp*nucleo=0'
c	 WRITE(*,2000)DOT_PRODUCT(dcomp,anuc) ; PAUSE'vérif' 

	 dcomp(14)=-DOT_PRODUCT(anuc,dcomp)/anuc(14) !Ex, cons baryons

c	   calcul du jacobien

	 IF(deriv)THEN	!jac(i,j):équation, i:élément j
	
c	  équation
c	  dcomp(1)=-(2.d0*r(1)*comp(1)+r(2)*comp(2)+r(5)*comp(5)
c	1 +r(7)*comp(6)+r(8)*comp(7)+r(9)*comp(8)+r(10)*comp(9)
c	2 +(r(11)+r(12))*comp(10)+r(13)*comp(11)
c	3 +r(14)*comp(12)+(r(15)+r(18))*comp(13)+(r(16)+r(17))*comp(17)
c	4 +(r(19)+r(20))*comp(15))*comp(1)
c	5 +2.d0*r(3)*comp(3)**2					!H1

	  jac(1,1)=-4.d0*r(1)*comp(1)-r(2)*comp(2)-r(5)*comp(5)
	1   -r(7)*comp(6)-r(8)*comp(7)-r(9)*comp(8)-r(10)*comp(9)
	2   -(r(11)+r(12))*comp(10)-r(13)*comp(11)-r(14)*comp(12)
	3   -(r(15)+r(18))*comp(13)-(r(16)+r(17))*comp(17)
	4   -(r(19)+r(20))*comp(15)		!d /H1
	  jac(1,2)=-r(2)*comp(1)				!d /H2
	  jac(1,3)=4.d0*r(3)*comp(3)				!d /He3
	  jac(1,5)=-r(5)*comp(1)				!d /Li7
	  jac(1,6)=-r(7)*comp(1)				!d /Be7
	  jac(1,7)=-r(8)*comp(1)				!d /C12
	  jac(1,8)=-r(9)*comp(1)				!d /C13
	  jac(1,9)=-r(10)*comp(1)				!d /N14
	  jac(1,10)=-(r(11)+r(12))*comp(1)			!d /N15
	  jac(1,11)=-r(13)*comp(1)				!d /O16
	  jac(1,12)=-r(14)*comp(1)				!d /O17
	  jac(1,13)=-(r(15)+r(18))*comp(1)			!d /Be9
	  jac(1,15)=-(r(19)+r(20))*comp(1)			!d /B11
	  jac(1,17)=-(r(16)+r(17))*comp(1)			!d /Li6
	  	  	 
	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(1,i)=jac(1,i)
	1  -(2.d0*drx(1,i)*comp(1)+drx(2,i)*comp(2)
	2  +drx(5,i)*comp(5)+drx(7,i)*comp(6)
	3  +drx(8,i)*comp(7)+drx(9,i)*comp(8)+drx(10,i)*comp(9)+(drx(11,i)
	4  +drx(12,i))*comp(10)+drx(13,i)*comp(11)+drx(14,i)*comp(12)
	5  +(drx(15,i)+drx(18,i))*comp(13)+(drx(16,i)+drx(17,i))*comp(17)
	6  +(drx(19,i)+drx(20,i))*comp(15))*comp(1)
	7  +2.d0*drx(3,i)*comp(3)**2				!H1	
	  ENDDO
	   	 
c	  équation dcomp(2)
c	  dcomp(2)=(r(1)*comp(1)-r(2)*comp(2))*comp(1)
c	1 +r(15)*comp(1)*comp(13)				!H2

	  jac(2,1)=2.d0*r(1)*comp(1)-r(2)*comp(2)+r(15)*comp(13)	!d /H1
	  jac(2,2)=-r(2)*comp(1)				!d /H2
	  jac(2,13)=r(15)*comp(1)				!d /Be9
	 
	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(2,i)=jac(2,i)
	1  +(drx(1,i)*comp(1)-drx(2,i)*comp(2))*comp(1)
	2  +drx(15,i)*comp(1)*comp(13)
	  ENDDO
	 	 	 
c	  équation dcomp(3)
c	  dcomp(3)=r(2)*comp(1)*comp(2)-(2.d0*r(3)*comp(3)
c	1 +r(4)*comp(4))*comp(3)+r(16)*comp(17)*comp(1)		!He3

	  jac(3,1)=r(2)*comp(2)+r(16)*comp(17)			!d /H1
	  jac(3,2)=r(2)*comp(1)					!d /H2
	  jac(3,3)=-2.d0*r(3)*comp(3)-r(4)*comp(4)		!d /He3
	  jac(3,4)=-r(4)*comp(3)				!d /He4
	  jac(3,17)=r(16)*comp(1)				!d /Li6	  
	 
	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(3,i)=jac(3,i)
	1  +drx(2,i)*comp(1)*comp(2)-(2.d0*drx(3,i)*comp(3)
	2  +drx(4,i)*comp(4))*comp(3)+drx(16,i)*comp(17)*comp(1)
	  ENDDO

c	  équation dcomp(4)
c	  dcomp(4)=(r(3)*comp(3)-r(4)*comp(4))*comp(3)
c	1 +(2.d0*(r(5)*comp(5)+r(7)*comp(6))+r(12)*comp(10)
c	2 +r(14)*comp(12)+2.d0*r(15)*comp(13)
c	3 +r(16)*comp(17)+3.d0*r(19)*comp(15))*comp(1)	!He4

	  jac(4,1)=2.d0*(r(5)*comp(5)+r(7)*comp(6))
	1 +r(12)*comp(10)+r(14)*comp(12)+2.d0*r(15)*comp(13)
	2 +r(16)*comp(17)+3.d0*r(19)*comp(15)		!d /H1
	  jac(4,3)=2.d0*r(3)*comp(3)-r(4)*comp(4)	!d /He3
	  jac(4,4)=-r(4)*comp(3)			!d /He4
	  jac(4,5)=2.d0*r(5)*comp(1)			!d /Li7
	  jac(4,6)=2.d0*r(7)*comp(1)			!d /Be7
	  jac(4,10)=r(12)*comp(1)			!d /N15
	  jac(4,12)=r(14)*comp(1)			!d /O17
	  jac(4,13)=2.d0*r(15)*comp(1)			!d /Be9
	  jac(4,15)=3.d0*r(19)*comp(1)			!d /B11
	  jac(4,17)=r(16)*comp(1)			!d /Li6

	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(4,i)=jac(4,i)
	1  +(drx(3,i)*comp(3)-drx(4,i)*comp(4))*comp(3)
	2  +(2.d0*(drx(5,i)*comp(5)+drx(7,i)*comp(6))
	3  +drx(12,i)*comp(10)+drx(14,i)*comp(12)+2.d0*drx(15,i)*comp(13)
	3  +drx(16,i)*comp(17)+3.d0*drx(19,i)*comp(15))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(5)
c	  dcomp(5)=-r(5)*comp(1)*comp(5)+r(6)*comp(6)*mue	!Li7

	  jac(5,1)=-r(5)*comp(5)				!d /H1
	  jac(5,5)=-r(5)*comp(1)				!d /Li7
	  jac(5,6)=r(6)*mue					!d /Be7
	 
	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(5,i)=jac(5,i)
	1    -drx(5,i)*comp(1)*comp(5)
	2    +drx(6,i)*comp(6)*mue+r(6)*comp(6)*dmuex(i)
	  ENDDO
	 	 
c	  équation dcomp(6)	!Be7
c	  dcomp(6)=r(4)*comp(3)*comp(4)-(r(6)*mue+r(7)*comp(1))*comp(6)
c	1 +r(17)*comp(17)*comp(1) 			!Be7

	  jac(6,1)=-r(7)*comp(6)+r(17)*comp(17)		!d /H1
	  jac(6,3)=r(4)*comp(4)				!d /He3
	  jac(6,4)=r(4)*comp(3)				!d /He4
	  jac(6,6)=-r(6)*mue-r(7)*comp(1)		!d /Be7
	  jac(6,17)=r(17)*comp(1)			!d /Li6

	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(6,i)=jac(6,i)
	1    +drx(4,i)*comp(3)*comp(4)-(drx(6,i)*mue
	2    +r(6)*dmuex(i)+drx(7,i)*comp(1))*comp(6)
	3    +drx(17,i)*comp(17)*comp(1)
	  ENDDO
	 
c	  équation dcomp(7)	 
c	  dcomp(7)=(-r(8)*comp(7)+r(12)*comp(10)	!C12
c	1 +r(20)*comp(15))*comp(1)

	  jac(7,1)=-r(8)*comp(7)+r(12)*comp(10)+r(20)*comp(15)	!d /H1
	  jac(7,7)=-r(8)*comp(1)			!d /C12
	  jac(7,10)=r(12)*comp(1)			!d /N15
	  jac(7,15)=r(20)*comp(1)			!d /B11
	 
	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(7,i)=jac(7,i)
	1  +(-drx(8,i)*comp(7)+drx(12,i)*comp(10)
	2  +drx(20,i)*comp(15))*comp(1)
	  ENDDO	 
	 	 
c	  équation dcomp(8)
c	  dcomp(8)=(r(8)*comp(7)-r(9)*comp(8))*comp(1)	!C13

	  jac(8,1)=r(8)*comp(7)-r(9)*comp(8)		!d /H1
	  jac(8,7)=r(8)*comp(1)				!d /C12
	  jac(8,8)=-r(9)*comp(1)			!d /C13

	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(8,i)=jac(8,i)
	1    +(drx(8,i)*comp(7)-drx(9,i)*comp(8))*comp(1)
	  ENDDO
	 
c	  équation dcomp(9)	!N14	 
c	  dcomp(9)=(r(9)*comp(8)-r(10)*comp(9)+r(14)*comp(12))*comp(1)

	  jac(9,1)=r(9)*comp(8)-r(10)*comp(9)+r(14)*comp(12)	!d /H1
	  jac(9,8)=r(9)*comp(1)					!d /C13
	  jac(9,9)=-r(10)*comp(1)				!d /N14
	  jac(9,12)=r(14)*comp(1)				!d /O17

	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(9,i)=jac(9,i)
	1    +(drx(9,i)*comp(8)-drx(10,i)*comp(9)
	2    +drx(14,i)*comp(12))*comp(1)
	  ENDDO
	 
c	  équation dcomp(10)
c	  dcomp(10)=(r(10)*comp(9)-(r(11)+r(12))*comp(10))*comp(1) !N15

	  jac(10,1)=r(10)*comp(9)-(r(11)+r(12))*comp(10)	!d /H1
	  jac(10,9)=r(10)*comp(1)				!d /N14
	  jac(10,10)=-(r(11)+r(12))*comp(1)			!d /N15

	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(10,i)=jac(10,i)
	1  +(drx(10,i)*comp(9)-(drx(11,i)+drx(12,i))*comp(10))*comp(1)
	  ENDDO	  
	 
c	  équation dcomp(11)	 
c	  dcomp(11)=(r(11)*comp(10)-r(13)*comp(11))*comp(1)	!O16

	  jac(11,1)=r(11)*comp(10)-r(13)*comp(11)		!d /H1
	  jac(11,10)=r(11)*comp(1)				!d /N15
	  jac(11,11)=-r(13)*comp(1)				!d /O16
	 
	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(11,i)=jac(11,i)	 
	1  +(drx(11,i)*comp(10)-drx(13,i)*comp(11))*comp(1)
	  ENDDO
	 	 
c	  équation dcomp(12)
c	  dcomp(12)=(r(13)*comp(11)-r(14)*comp(12))*comp(1)	!O17

	  jac(12,1)=r(13)*comp(11)-r(14)*comp(12)		!d /H1
	  jac(12,11)=r(13)*comp(1)				!d /O16
	  jac(12,12)=-r(14)*comp(1)				!d /O17

	  DO i=1,17		!dépendances dues à l'effet d'écran
	   jac(12,i)=jac(12,i)
	1    +(drx(13,i)*comp(11)-drx(14,i)*comp(12))*comp(1)	  
	  ENDDO		 
	 
c	  équation dcomp(13)	 
c	  dcomp(13)=-(r(15)+r(18))*comp(1)*comp(13)	!Be9	 

	  jac(13,1)=-(r(15)+r(18))*comp(13)		!d /H1
	  jac(13,13)=-(r(15)+r(18))*comp(1)		!d /Be9

	  DO i=1,17	!dépendances dues à l'effet d'écran
	   jac(13,i)=jac(13,i)
	1   -(drx(15,i)+drx(18,i))*comp(1)*comp(13)
	  ENDDO
	  
c	  équation dcomp(14)	  
c	  dcomp(14)=0.d0 					!Ex
c	  jac(14,14)=1.d0
 
c	  équation dcomp(15)	  
c	  dcomp(15)=-(r(19)+r(20))*comp(1)*comp(15)		!B11

	  jac(15,1)=-(r(19)+r(20))*comp(15)			!d /H1
	  jac(15,15)=-(r(19)+r(20))*comp(1)			!d /B11

c	  équation dcomp(16)	  
c	  dcomp(16)=0.dO
c	  jac(16,16)=0.d0				!Fe56
	  
	  DO i=1,17	!dépendances dues à l'effet d'écran
	   jac(15,i)=jac(15,i)
	1    -(drx(19,i)+drx(20,i))*comp(1)*comp(15)
	  ENDDO
	  	  	  	  	 
c	  équation dcomp(17)
c	  dcomp(17)=(-(r(16)+r(17))*comp(17)+r(18)*comp(13))*comp(1) !Li6
	  
	  jac(17,1)=-(r(16)+r(17))*comp(17)+r(18)*comp(13)	!d /H1
	  jac(17,13)=r(18)*comp(1)				!d /Be9	  
	  jac(17,17)=-(r(16)+r(17))*comp(1)			!d /Li6

	  DO i=1,17	!dépendances dues à l'effet d'écran
	   jac(17,i)=jac(17,i)
	1    +(-(drx(16,i)+drx(17,i))*comp(17)+drx(18,i)*comp(13))*comp(1)
	  ENDDO

c conservation des baryons

c	  dcomp(14)=-SUM(anuc*dcomp)/anuc(14)

	  DO j=1,nchim
	   DO i=1,nchim
	    IF(i == 14)CYCLE
	    jac(14,j)=jac(14,j)+anuc(i)*jac(i,j)
	   ENDDO
	   jac(14,j)=-jac(14,j)/anuc(14)	   
	  ENDDO

c	  unités de temps pour intégration temporelle

	  jac=jac*secon6

	 ENDIF
	
	 dcomp=dcomp*secon6

c	 calcul de la production d'énergie nucléaire et dérivées

	CASE(3)
	
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0	 
	 IF(t <= t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	 
c	 mue:nombre d'electrons / mole /g = 1/poids mol. moy. par e-

	 epsilon(2)=(q(1)*comp(1)+q(2)*comp(2)+q(5)*comp(5)+q(7)*comp(6)
	1 +(q(15)+q(18))*comp(13)+(q(16)+q(17))*comp(17)
	2 +(q(19)+q(20))*comp(15))*comp(1)
	3 +(q(3)*comp(3)+q(4)*comp(4))*comp(3)+q(6)*mue*comp(6)
	 epsilon(3)=(q(8)*comp(7)+q(9)*comp(8)+q(10)*comp(9)
	1 +(q(11)+q(12))*comp(10)+q(13)*comp(11)+q(14)*comp(12))*comp(1)
	 epsilon(1)=SUM(epsilon(2:3))

	 IF(deriv)THEN
	
	  et=(dqt(1)*comp(1)+dqt(2)*comp(2)+dqt(5)*comp(5)+dqt(7)*comp(6))
	1 *comp(1)+(dqt(3)*comp(3)+dqt(4)*comp(4))*comp(3)
	2 +dqt(6)*mue*comp(6)+(dqt(8)*comp(7)+dqt(9)*comp(8)
	3 +dqt(10)*comp(9)+(dqt(11)+dqt(12))*comp(10)+dqt(13)*comp(11)
	4 +dqt(14)*comp(12)+(dqt(15)+dqt(18))*comp(13)
	5 +(dqt(16)+dqt(17))*comp(17)+(dqt(19)+dqt(20))*comp(15))*comp(1)
		
	  ero=(dqo(1)*comp(1)+dqo(2)*comp(2)+dqo(5)*comp(5)+dqo(7)*comp(6))
	1 *comp(1)+(dqo(3)*comp(3)+dqo(4)*comp(4))*comp(3)
	2 +dqo(6)*mue*comp(6)+(dqo(8)*comp(7)+dqo(9)*comp(8)
	3 +dqo(10)*comp(9)+(dqo(11)+dqo(12))*comp(10)+dqo(13)*comp(11)
	4 +dqo(14)*comp(12)+(dqo(15)+dqo(18))*comp(13)
	5 +(dqo(16)+dqo(17))*comp(17)+(dqo(19)+dqo(20))*comp(15))*comp(1)
		
	  ex(1)=2.*q(1)*comp(1)+q(2)*comp(2)+q(5)*comp(5)+q(7)*comp(6)
	1 +q(8)*comp(7)+q(9)*comp(8)+q(10)*comp(9)
	2 +(q(11)+q(12))*comp(10)+q(13)*comp(11)+q(14)*comp(12)
	3 +(q(15)+q(18))*comp(13)+(q(16)+q(17))*comp(17)
	4 +(q(19)+q(20))*comp(15)	
	  ex(2)=q(2)*comp(1)
	  ex(3)=2.*q(3)*comp(3)+q(4)*comp(4)
	  ex(4)=q(4)*comp(3)+q(16)*comp(13)
	  ex(5)=q(5)*comp(1)
	  ex(6)=q(6)*mue+q(7)*comp(1)
	  ex(7)=q(8)*comp(1)
	  ex(8)=q(9)*comp(1)
	  ex(9)=q(10)*comp(1)
	  ex(10)=(q(11)+q(12))*comp(1)
	  ex(11)=q(13)*comp(1)
	  ex(12)=q(14)*comp(1)
	  ex(13)=(q(15)+q(18))*comp(1)
	  ex(15)=(q(19)+q(20))*comp(1)
	  ex(17)=(q(16)+q(17))*comp(1)
	  DO i=1,nchim	!contributions des écrans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*comp(2)
	1  +dqx(5,i)*comp(5)+dqx(7,i)*comp(6))*comp(1)
	2  +(dqx(3,i)*comp(3)+dqx(4,i)*comp(4))*comp(3)
	3  +dqx(6,i)*mue*comp(6)+(dqx(8,i)*comp(7)
	4  +dqx(9,i)*comp(8)+dqx(10,i)*comp(9)
	5  +(dqx(11,i)+dqx(12,i))*comp(10)
	6  +dqx(13,i)*comp(11)+dqx(14,i)*comp(12)
	7  +(dqx(15,i)+dqx(18,i))*comp(13)+(dqx(16,i)+dqx(17,i))*comp(17)
	8  +(dqx(19,i)+dqx(20,i))*comp(15))*comp(1)

	  ENDDO
	 
	 ENDIF	!deriv

c	 production de neutrinos

	 CASE(4)
	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*comp(6)/amu
	  b8e=r(7)*comp(1)*comp(6)/amu ; n13e=r(8)*comp(1)*comp(7)/amu
	  o15e=r(10)*comp(1)*comp(9)/amu ; f17e=r(13)*comp(1)*comp(11)/amu
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0 ; n13e=0.d0
	  o15e=0.d0 ; f17e=0.d0
	 ENDIF
	 
	CASE DEFAULT
	 PRINT*,'ppcno12Be, fait ne peut valoir que 1, 2, 3 ou 4'
	 PRINT*,'ERREUR fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; STOP
	 
	END SELECT
	
	RETURN

	END SUBROUTINE ppcno12BeBFe
