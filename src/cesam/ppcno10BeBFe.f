
c**************************************************************************

	SUBROUTINE ppcno10BeBFe(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)
	
c       routine private du module mod_nuc

c	cycles PP et CNO
c	cf. Clayton p. 380, 392 et 430,

c	éléments pris en compte:
c	H1, He3, He4, Li7, C12, C13, N14, N15, O16, O17, Fe56, Ex,
c	Li6, Be9, B11
c	H2 et Be7 à l'équilibre, Ex est l'élément moyen de complément
c	Fe56 et Ex n'intéressent que la diffusion

c	un premier appel a rq_reac initialise et définit le nb.
c	d'éléments chimiques pour lesquels les reac. nuc. sont tabulées
c	dans ppcno10BeBFe on ajoute Ex, soit nchim+1, puis

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
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
c	jac : jacobien (unité de temps : 10**6 ans)
c	epsilon, et, ero, ex : énergie thermonucléaire (unite de temps : s)
c			   : et dérivées /t, ro ,X

c	Neutrinos
c	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
c	hhe réaction : H1(p,e+ nu)H2
c	be7e réaction : Be7(e-,nu g)Li7
c	b8e réaction : B8(,e+ nu)Be8
c	n13e réaction : N13(,e+ nu)C13
c	o15e réaction : O15(e+,nu)N15 
c	f17e réaction : F17(,e+ nu)O17

c initialisations
c	ab_min : abondances négligeables
c	ab_ini : abondances initiales

c	r(1) : réaction H1(p,e+ nu)H2			PP
c	r(2) : réaction H2(p,g)He3			H2 équilibre
c	r(3) : réaction He3(He3,2H)He4
c	r(4) : réaction He4(He3,g)Be7			Be7 équilibre
c	r(5) : réaction Li7(p,He4)He4
c	r(6) : réaction Be7(e-,nu g)Li7			Be7 équilibre
c	r(7) : réaction Be7(p,g)B8(,e+ nu)Be8(,He4)He4	Be7 équilibre

c	r(8) : réaction C12(p,g)N13(,e+ nu)C13		CNO
c	r(9) : réaction C13(p,g)N14
c	r(10) : réaction N14(p,g)O15(e+,nu)N15
c	r(11) : réaction N15(p,g)O16
c	r(12) : réaction N15(p,He4)C12
c	r(13) : réaction O16(p,g)F17(,e+ nu)O17
c	r(14) : réaction O17(p,He4)N14

c	r(15) : réaction Be9(p,d)2He4	autres réactions H2 équilibre
c	r(16) : réaction Li6(p,He3)He4
c	r(17) : réaction Li6(p,g)Be7		Be7 équilibre
c	r(18) : réaction Be9(p,a)Li6
c	r(19) : réaction B11(p,a)2He4
c	r(20) : réaction B11(p,g)C12

c indices des éléments
c	H1 : 1
c	He3 : 2
c	He4 : 3
c	Li7 : 4
c	C12 : 5
c	C13 : 6
c	N14 : 7
c	N15 : 8
c	O16 : 9
c	O17 : 10
c	Fe56 : 11
c	Ex  : 12
c	Li6 : 13
c	Be9 : 14
c	B11 : 15

c----------------------------------------------------------------------
	
	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon,
	1 ife56, ihe4, ili7, i_ex, langue, nchim, nom_elem, nom_xheavy, nucleo,
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
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: anuc, comp_dex,
	1 dmuex, dh2x, denx, dbe7x, drt, dro, r, q, dqt, dqo
	REAL (kind=dp) :: mue, nbz, h2, dh2be9, dh2h, den, be7, dbe7he3,
	1 dbe7he4, dbe7li6, dbe7mue, dbe7h, charge_ex , mass_ex , sum_a
		
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
	 
c définition de nchim : nombre d'isotopes chimiques dont on
c	 calcule l'abondance

	 nchim=15 ; ili7=4
	 
c ordre des isotopes :	 
c	 H1(1), He3(2), He4,(3) Li7(4), C12(5), C13(6), N14(7), N15(8),
c	 O16(9), O17(10), Fe(11), Ex(12), Li6(13), Be9(14), B11(15)

c	 appel d'initialisation pour tabulation des réactions nucléaires
c	 allocations fictives

	 ALLOCATE(drx(0,0),dqx(0,0),r(0),drt(0),dro(0),q(0),
	1 dqt(0),dqo(0),dmuex(0))
	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 	 	 
	CASE(1)

c détermination des abondances initiales
c	 He3+He4=Y0
c	 Z0 = somme des éléments plus lourds que hélium
c	 dans Z rapports en nombre

	 CALL abon_ini
	 	 
c	 Ex : élément fictif moyenne des éléments # Li, Be, B, CNO, Fe

	 charge_ex=0.d0 ; mass_ex=0.d0 ; sum_a=0.d0
	 B1: DO i=3,nelem_ini		!à partir de Li
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


c élément chimique reliquat
	 i_ex=12		!indice de l'élément chimique reliquat
	 nucleo(12)=mass_ex	!nucleo de l'élément chimique reliquat
	 zi(12)=charge_ex	!charge de l'élément chimique reliquat
	 i=NINT(charge_ex)
	 nom_elem(12)=elem(i)//text !nom elem. chim. rel.
	 nom_xheavy=nom_elem(12)
 	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1002)TRIM(nom_elem(12)),nint(mass_ex),nint(charge_ex)
	  WRITE(2,1002)TRIM(nom_elem(12)),nint(mass_ex),nint(charge_ex)	 
1002	  FORMAT(a,': fictitious species /= CNO, of mass : ',i3,/,
	1 'and charge :',i3)	 
	 CASE DEFAULT	 
	  WRITE(*,2)TRIM(nom_elem(12)),nint(mass_ex),nint(charge_ex)
	  WRITE(2,2)TRIM(nom_elem(12)),nint(mass_ex),nint(charge_ex)	 
2	  FORMAT(a,': élément fictif /= CNO, de masse : ',i3,/,
	1 'et de charge :',i3)
	 END SELECT	 	 
c	 PRINT*,nchim ; WRITE(*,2000)nucleo ; PAUSE'nchim' 
	
	 ALLOCATE(a(nchim,nchim),indpc(nchim),b(1,nchim))
	 a=0.d0 ; b=0.d0 ; indpc=1	
			
	 a(1,1)=nucleo(1)  	!H1 abondance de H
	 b(1,1)=x0
	
	 a(2,2)=nucleo(2)  	!He3	abondance de He
	 a(2,3)=nucleo(3)	!He4
	 b(1,2)=y0

	 DO j=4,nchim		!abondances dans Z
	  a(3,j)=nucleo(j)	!somme j >= 5 comp(j)*nucleo(j)=Z0
	  a(4,j)=-abon_rela(6)		!somme comp(i) C, C/Z
	  a(5,j)=-abon_rela(7)		!somme comp(i) N, N/Z	 
	  a(6,j)=-abon_rela(8)		!somme comp(i) O, O/Z
	  a(11,j)=-abon_rela(3)		!somme comp(i) Li, Li/Z
	  a(12,j)=-abon_rela(26)	!somme comp(i) Fe, Fe/Z
	  a(13,j)=-abon_rela(4) 	!somme comp(i) Be, Be/Z
	  a(14,j)=-abon_rela(5)		!somme comp(i) B, B/Z	  	  
	 ENDDO
		 		
	 b(1,3)=z0		!Z
	
	 a(4,5)=a(4,5)+1.d0	!C12	coefficients des isotopes	 		
	 a(4,6)=a(4,6)+1.d0	!C13	équation 4 pour C
	
	 a(5,7)=a(5,7)+1.d0	!N14	équation 5 pour N 		
	 a(5,8)=a(5,8)+1.d0	!N15
	
	 a(6,9)=a(6,9)+1.d0	!O16	équation 6 pour O 		
	 a(6,10)=a(6,10)+1.d0	!O17
	 
	 a(11,4)=a(11,4)+1.d0	!Li7	équation 11 pour Li 
	 a(11,13)=a(11,13)+1.d0	!Li6	 	 
	 a(12,11)=a(12,11)+1.d0	!Fe	équation 12 pour Fe
	 a(13,14)=a(13,14)+1.d0	!Be9	équation 13 pour Be
	 a(14,15)=a(14,15)+1.d0	!B11	équation 14 pour B	  	 
	
c	 rapports isotopiques équations 7 pour He, 8 pour C...15 pour Li
	
	 a(7,2)=1.d0		!He3
	 a(7,3)=-he3she4z	!He3/He4, H2 est dans He3

	 a(8,6)=1.d0		!C13
	 a(8,5)=-c13sc12	!C13/C12
	
	 a(9,8)=1.d0		!N15
	 a(9,7)=-n15sn14	!N15/N14
	
	 a(10,10)=1.d0		!O17
	 a(10,9)=-o17so16	!O17/O16
	 
	 a(15,13)=1.d0		!Li6
	 a(15,4)=-li6sli7	!Li6/Li7
			
c	 PRINT*,nchim
c	 DO i=1,nchim
c	  WRITE(*,2002)a(i,1:nchim),b(1,i)
c	 ENDDO

	 CALL gauss_band(a,b,indpc,nchim,nchim,nchim,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'ppcno10BeBFe, matrice calcul abondances non inversible'
	  PRINT*,'ARRET'
	  stop
	 ENDIF	

c allocations diverses

	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),
	2 dmuex(nchim),dh2x(nchim),denx(nchim),dbe7x(nchim))

c abondances initiales et abondances négligeables

	 comp(1:nchim)=max(1.d-29,b(1,1:nchim))
	 ab_ini(1:nchim)=comp(1:nchim)*nucleo(1:nchim)		 
	 ab_min=ab_ini*fmin_abon
	
c nombre/volume des métaux dans Z, indice de Fe56

	 nbz=sum(comp(ihe4+1:nchim)) ; ife56=11
	 
c abondances en DeX, H=12

	 ALLOCATE(comp_dex(nchim))
	 comp_dex=12.d0+LOG10(comp/comp(1))
	 		
c écritures

	 WRITE(2,6) ; WRITE(*,6) 
6	 FORMAT(/,'Réactions thermonucléaires des cycles PP, CNO',/)
	 WRITE(2,7)nreac ; WRITE(*,7)nreac 
7	 FORMAT('nombre de réactions : ',i3)
	 WRITE(2,8)nreac ; WRITE(*,8)nchim
8	 FORMAT('nombre d''éléments chimiques : ',i3)
	 WRITE(2,20)x0,y0,z0,z0/x0 ; WRITE(*,20)x0,y0,z0,z0/x0	
20	 FORMAT(/,'abondances initiales déduites de X0=',es10.3,
	1 ', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,
	2 ', H1=X0, He3+He4=Y0',/
	3 'Z0 = 1-X0-Y0 = Li7+Be9+B11+C12+C13+N14+N15+O16+O17+Ex+Fe',/)
	 WRITE(2,1)ab_ini(1:nchim) ; WRITE(*,1)ab_ini(1:nchim)
1	 FORMAT('abondances initiales / gramme:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', Li7:',es10.3,', C12:',es10.3,/,'C13:',es10.3,
	3 ', N14:',es10.3,', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,/,
	4 'Fe56:',es10.3,', Ex:',es10.3,', Li6:',es10.3,', Be9:',
	5 es10.3,', B11:',es10.3)
	 WRITE(2,9)comp_dex ; WRITE(*,9)comp_dex	
9	 FORMAT(/,'Abondances initiales en nombre: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', Li7:',es10.3,', C12:',es10.3,/,'C13:',es10.3,
	3 ', N14:',es10.3,', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,/,
	4 'Fe56:',es10.3,', Ex:',es10.3,', Li6:',es10.3,', Be9:',
	5 es10.3,', B11:',es10.3)
	 WRITE(2,21)(comp(4)+comp(13))/nbz,	!Li/Z
	1 (comp(5)+comp(6))/nbz,		!C/Z
	2 (comp(7)+comp(8))/nbz,		!N/Z
	3 (comp(9)+comp(10))/nbz,		!O/Z
	4 comp(12)/nbz,comp(11)/nbz,		!Ex/Z, Fe/Z
	5 comp(14)/nbz,				!Be/Z
	6 comp(15)/nbz				!B/Z		
	 WRITE(*,21)(comp(4)+comp(13))/nbz,	!Li/Z
	1 (comp(5)+comp(6))/nbz,		!C/Z
	2 (comp(7)+comp(8))/nbz,		!N/Z
	3 (comp(9)+comp(10))/nbz,		!O/Z
	4 comp(12)/nbz,comp(11)/nbz,		!Ex/Z, Fe/Z
	5 comp(14)/nbz,				!Be/Z
	6 comp(15)/nbz				!B/Z		
21	 FORMAT(/,'rapports en nombre dans Z:',/,'Li/Z:',es10.3,
	1 ', C/Z:',es10.3,', N/Z:',es10.3,', O/Z:',es10.3,/,'Fe/Z:',es10.3,
	2 ', Ex/Z:',es10.3,', Be/Z:',es10.3,', B/Z:',es10.3,/)
	 WRITE(2,14)ab_min(1:nchim) ; WRITE(*,14)ab_min(1:nchim)
14	 FORMAT('abondances / gramme négligeables:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,
	2 ', Li7:',es10.3,', C12:',es10.3,/,'C13:',es10.3,
	3 ', N14:',es10.3,', N15:',es10.3,', O16:',es10.3,', O17:',es10.3,/,
	4 'Fe56:',es10.3,', Ex:',es10.3,', Li6:',es10.3,', Be9:',
	5 es10.3,', B11:',es10.3)
	 WRITE(2,11) ; WRITE(*,11)
11	 FORMAT(/,'H2, Be7 à l''équilibre')
	 WRITE(2,12) ; WRITE(*,12)
12	 FORMAT(/,'on utilise une table')	 
	 WRITE(2,13) ; WRITE(*,13)
13	 FORMAT(/,'évol. temporelle, test de précision sur H1 et He4')

	 DO i=1,nchim
	  ab_min(i)=ab_min(i)/nucleo(i)
	  anuc(i)=anint(nucleo(i))		!nombre atomique
	 ENDDO

c nettoyage

	 DEALLOCATE(a,b,comp_dex,indpc)

c réactions

	CASE(2)
	 dcomp=0.d0 ; jac=0.d0
	
	 IF(t < t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	
c	 PRINT*,'comp' ; WRITE(*,2000)comp(1:nchim)
c	 PRINT*,'réactions' ; WRITE(*,2000)r(1:nreac)

c H2 et Be7 à l'équilibre
c	 H2, H1(p,e+ nu)H2, H2(p,g)He3, Be9(p,d)2He4

	 dh2h=r(1)/r(2) ; dh2be9=r(15)/r(2)
	 h2=dh2h*comp(1)+dh2be9*comp(14)

	 den=r(6)*mue+r(7)*comp(1)
	 IF(den /= 0.d0)THEN
	  be7=(r(4)*comp(2)*comp(3)+r(17)*comp(1)*comp(13))/den 
	  dbe7he3=r(4)*comp(3)/den ; dbe7he4=r(4)*comp(2)/den
	  dbe7mue=-be7*r(6)/den ; dbe7h=(-be7*r(7)+r(17)*comp(13))/den
	  dbe7li6=r(17)*comp(1)/den
	 ELSE
	  be7=0.d0 ; dbe7he3=0.d0 ; dbe7mue=0.d0 ; dbe7he4=0.d0
	  dbe7h= 0.d0 ; dbe7li6=0.d0
	 ENDIF
	
c	 WRITE(*,2000)dh2h,h2,be7,dbe7he3,dbe7he4,dbe7mue,dbe7h ; PAUSE

c	 équations d'évolution

	 dcomp(1)=-(2.d0*r(1)*comp(1)+r(2)*h2+r(5)*comp(4)	!H1
	1 +r(7)*be7+r(8)*comp(5)+r(9)*comp(6)+r(10)*comp(7)
	2 +(r(11)+r(12))*comp(8)+r(13)*comp(9)+r(14)*comp(10)
	3 +(r(15)+r(18))*comp(14)+(r(16)+r(17))*comp(13)
	4 +(r(19)+r(20))*comp(15))*comp(1)+2.d0*r(3)*comp(2)**2

	 dcomp(2)=r(2)*comp(1)*h2-(2.d0*r(3)*comp(2)
	1 +r(4)*comp(3))*comp(2)+r(16)*comp(13)*comp(1)		!He3
		
	 dcomp(3)=(r(3)*comp(2)-r(4)*comp(3))*comp(2)
	1 +(2.d0*(r(5)*comp(4)+r(7)*be7)+r(12)*comp(8)
	2 +r(14)*comp(10)+2.d0*r(15)*comp(14)+r(16)*comp(13)
	3 +3.d0*r(19)*comp(15))*comp(1)				!He4
		
	 dcomp(4)=-r(5)*comp(1)*comp(4)+r(6)*be7*mue		!Li7
	 
	 dcomp(5)=(-r(8)*comp(5)+r(12)*comp(8)
	1 +r(20)*comp(15))*comp(1)				!C12
	 
	 dcomp(6)=(r(8)*comp(5)-r(9)*comp(6))*comp(1)		!C13
	 
	 dcomp(7)=(r(9)*comp(6)-r(10)*comp(7)+r(14)*comp(10))*comp(1) !N14
	 
	 dcomp(8)=(r(10)*comp(7)-(r(11)+r(12))*comp(8))*comp(1)	!N15
	 
	 dcomp(9)=(r(11)*comp(8)-r(13)*comp(9))*comp(1)		!O16
	 
	 dcomp(10)=(r(13)*comp(9)-r(14)*comp(10))*comp(1)	!O17
	 
c	 dcomp(11)=0.d0 					!Fe56	 
c	 dcomp(12)=0.d0 					!Ex

	 dcomp(13)=(-(r(16)+r(17))*comp(13)+r(18)*comp(14))*comp(1) !Li6

	 dcomp(14)=-(r(15)+r(18))*comp(14)*comp(1)	!Be9

	 dcomp(15)=-(r(19)+r(20))*comp(15)*comp(1)	!B11

c	   Pour vérifications SUM dcomp*nucleo=0

c	 PRINT*,'ppcno10BeBFe, vérifications SUM dcomp*nucleo=0'
c	 WRITE(*,2000)DOT_PRODUCT(dcomp,anuc) ; PAUSE'vérif' 

c	 conservation des baryons

	 dcomp(12)=-DOT_PRODUCT(anuc,dcomp)/anuc(12) !Ex, cons baryons	

c calcul du jacobien
 
	 IF(deriv)THEN	!jac(i,j) : équation, j : élément i
	
c équation dcomp(1)
c	 dcomp(1)=-(2.d0*r(1)*comp(1)+r(2)*h2+r(5)*comp(4)	!H1
c	 1 +r(7)*be7+r(8)*comp(5)+r(9)*comp(6)+r(10)*comp(7)
c	 2 +(r(11)+r(12))*comp(8)+r(13)*comp(9)+r(14)*comp(10)
c	 3 +(r(15)+r(18))*comp(14)+(r(16)+r(17))*comp(13)
c	 4 +(r(19)+r(20))*comp(15))*comp(1)+2.d0*r(3)*comp(2)**2

	  jac(1,1)=-(4.d0*r(1)+r(2)*dh2h)*comp(1)-r(2)*h2-r(5)*comp(4)
	1 -r(7)*be7-r(7)*comp(1)*dbe7h-r(8)*comp(5)-r(9)*comp(6)
	2 -r(10)*comp(7)-(r(11)+r(12))*comp(8)-r(13)*comp(9)
	3 -r(14)*comp(10)-(r(15)+r(18))*comp(14)-(r(16)+r(17))*comp(13)
	4 -(r(19)+r(20))*comp(15)	!d /H1	
	  jac(1,2)=-r(7)*comp(1)*dbe7he3+4.d0*r(3)*comp(2)	!d /He3	  
	  jac(1,3)=-r(7)*comp(1)*dbe7he4			!d /He4
	  jac(1,4)=-r(5)*comp(1)				!d /Li7
	  jac(1,5)=-r(8)*comp(1)				!d /C12
	  jac(1,6)=-r(9)*comp(1)				!d /C13
	  jac(1,7)=-r(10)*comp(1)				!d /N14
	  jac(1,8)=-(r(11)+r(12))*comp(1)			!d /N15
	  jac(1,9)=-r(13)*comp(1)				!d /O16
	  jac(1,10)=-r(14)*comp(1)				!d /O17
	  jac(1,13)=-(r(7)*dbe7li6+r(16)+r(17))*comp(1)		!d /Li6	  
	  jac(1,14)=-(r(2)*dh2be9+r(15)+r(18))*comp(1)		!d /Be9
	  jac(1,15)=-(r(19)+r(20))*comp(1)			!d /B11
	  	 
	  DO i=1,nchim		!dépendances effet d'écran et be7/muex
	   jac(1,i)=jac(1,i)-r(7)*dbe7mue*dmuex(i)*comp(1)
	1  -(2.d0*drx(1,i)*comp(1)+drx(2,i)*h2
	2  +drx(5,i)*comp(4)+drx(7,i)*(be7+dbe7mue*dmuex(i))
	3  +drx(8,i)*comp(5)+drx(9,i)*comp(6)
	4  +drx(10,i)*comp(7)+(drx(11,i)+drx(12,i))*comp(8)
	5  +drx(13,i)*comp(9)+drx(14,i)*comp(10)
	6  +(drx(15,i)+drx(18,i))*comp(14)+(drx(16,i)+drx(17,i))*comp(13)
	7  +(drx(19,i)+drx(20,i))*comp(15))*comp(1)
	8  +2.d0*drx(3,i)*comp(2)**2
	  ENDDO
	   	 	 	 	 
c	  équation dcomp(2)
c	  dcomp(2)=r(2)*comp(1)*h2-(2.d0*r(3)*comp(2)
c	1 +r(4)*comp(3))*comp(2)+r(16)*comp(13)*comp(1)		!He3

	  jac(2,1)=r(2)*(h2+comp(1)*dh2h)+r(16)*comp(13)	!d /H1
	  jac(2,2)=-4.d0*r(3)*comp(2)-r(4)*comp(3)		!d /He3
	  jac(2,3)=-r(4)*comp(2)				!d /He4
	  jac(2,13)=r(16)*comp(1)				!d /Li6
	  jac(2,14)=r(2)*comp(1)*dh2be9				!d /Be9
	  	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(2,i)=jac(2,i)
	1  +drx(2,i)*comp(1)*h2-(2.d0*drx(3,i)*comp(2)
	2  +drx(4,i)*comp(3))*comp(2)+drx(16,i)*comp(13)*comp(1)
	  ENDDO

c équation dcomp(3)
c	  dcomp(3)=(r(3)*comp(2)-r(4)*comp(3))*comp(2)
c	1 +(2.d0*(r(5)*comp(4)+r(7)*be7)+r(12)*comp(8)
c	2 +r(14)*comp(10)+2.d0*r(15)*comp(14)+r(16)*comp(13)
c	3 +3.d0*r(19)*comp(15))*comp(1)				!He4

	  jac(3,1)=2.d0*(r(5)*comp(4)+r(7)*be7+r(7)*dbe7h*comp(1))
	1 +r(12)*comp(8)+r(14)*comp(10)+2.d0*r(15)*comp(14)
	2 +r(16)*comp(13)+3.d0*r(19)*comp(15)			!d /H1	
	  jac(3,2)=2.d0*r(3)*comp(2)-r(4)*comp(3)
	1 +2.d0*r(7)*dbe7he3*comp(1) 				!d /He3
	  jac(3,3)=-r(4)*comp(2)+2.d0*r(7)*dbe7he4*comp(1)	!d /He4
	  jac(3,4)=2.d0*r(5)*comp(1)				!d /Li7
	  jac(3,8)=r(12)*comp(1)				!d /N15
	  jac(3,10)=r(14)*comp(1)				!d /O17
	  jac(3,13)=(r(7)*dbe7li6+r(16))*comp(1)		!d /Li6	  
	  jac(3,14)=2.d0*r(15)*comp(1)				!d /Be9	  
	  jac(3,15)=3.d0*r(19)*comp(1)				!d /B11	
	   
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(3,i)=jac(3,i)   
	1  +(drx(3,i)*comp(2)-drx(4,i)*comp(3))*comp(2)
	2  +(2.d0*(drx(5,i)*comp(4)+drx(7,i)*be7
	3  +2.d0*r(7)*dbe7mue*dmuex(i))+drx(12,i)*comp(8)
	4  +drx(14,i)*comp(10)+drx(16,i)*comp(13)
	5  +3.d0*drx(19,i)*comp(15))*comp(1)	
	  ENDDO
	 	 
c équation dcomp(4)
c	  dcomp(4)=-r(5)*comp(1)*comp(4)+r(6)*be7*mue	!Li7

	  jac(4,1)=-r(5)*comp(4)+r(6)*dbe7h*mue   	!d /H1
	  jac(4,2)=r(6)*dbe7he3*mue			!d /He3	 
	  jac(4,3)=r(6)*dbe7he4*mue			!d /He4	 
	  jac(4,4)=-r(5)*comp(1)			!d /Li7
	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(4,i)=jac(4,i)
	1  -drx(5,i)*comp(1)*comp(4)
	2  +drx(6,i)*be7*mue+r(6)*(dbe7mue*mue+be7)*dmuex(i)
	  ENDDO
	 	 	 
c équation dcomp(5)	 
c	  dcomp(5)=(-r(8)*comp(5)+r(12)*comp(8)+r(20)*comp(15))*comp(1)!C12

	  jac(5,1)=-r(8)*comp(5)+r(12)*comp(8)+r(20)*comp(15)	!d /H1
	  jac(5,5)=-r(8)*comp(1)			!d /C12
	  jac(5,8)=r(12)*comp(1)			!d /N15
	  jac(5,15)=r(20)*comp(1)			!d /B11
	  	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(5,i)=jac(5,i)
	1  +(-drx(8,i)*comp(5)+drx(12,i)*comp(8)
	2  +drx(20,i)*comp(15))*comp(1)
	  ENDDO
	 	 
c équation dcomp(6)
c	  dcomp(6)=(r(8)*comp(5)-r(9)*comp(6))*comp(1)		!C13

	  jac(6,1)=r(8)*comp(5)-r(9)*comp(6)			!d /H1
	  jac(6,5)=r(8)*comp(1)					!d /C12
	  jac(6,6)=-r(9)*comp(1)				!d /C13

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(6,i)=jac(6,i)
	1  +(drx(8,i)*comp(5)-drx(9,i)*comp(6))*comp(1)
	  ENDDO
	 
c équation dcomp(7)	!N14	 
c	  dcomp(7)=(r(9)*comp(6)-r(10)*comp(7)+r(14)*comp(10))*comp(1)

	  jac(7,1)=r(9)*comp(6)-r(10)*comp(7)+r(14)*comp(10)	!d /H1
	  jac(7,6)=r(9)*comp(1)				!d /C13
	  jac(7,7)=-r(10)*comp(1)				!d /N14
	  jac(7,10)=r(14)*comp(1)				!d /O17

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(7,i)=jac(7,i)
	1  +(drx(9,i)*comp(6)-drx(10,i)*comp(7)+drx(14,i)*comp(10))*comp(1)
	  ENDDO
	 
c équation dcomp(8)	 
c	  dcomp(8)=(r(10)*comp(7)-(r(11)+r(12))*comp(8))*comp(1)  !N15

	  jac(8,1)=r(10)*comp(7)-(r(11)+r(12))*comp(8)		!d /H1
	  jac(8,7)=r(10)*comp(1)				!d /N14
	  jac(8,8)=-(r(11)+r(12))*comp(1)			!d /N15

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(8,i)=jac(8,i)
	1  +(drx(10,i)*comp(7)-(drx(11,i)+drx(12,i))*comp(8))*comp(1)
	  ENDDO	  
	 
c équation dcomp(9)	 
c	  dcomp(9)=(r(11)*comp(8)-r(13)*comp(9))*comp(1)	!O16

	  jac(9,1)=r(11)*comp(8)-r(13)*comp(9)			!d /H1
	  jac(9,8)=r(11)*comp(1)				!d /N15
	  jac(9,9)=-r(13)*comp(1)				!d /O16
	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(9,i)=jac(9,i)	 
	1  +(drx(11,i)*comp(8)-drx(13,i)*comp(9))*comp(1)
	  ENDDO
	 	 
c équation dcomp(10)
c	  dcomp(10)=(r(13)*comp(9)-r(14)*comp(10))*comp(1)	!O17

	  jac(10,1)=r(13)*comp(9)-r(14)*comp(10)		!d /H1
	  jac(10,9)=r(13)*comp(1)				!d /O16
	  jac(10,10)=-r(14)*comp(1)				!d /O17

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(10,i)=jac(10,i)
	1  +(drx(13,i)*comp(9)-drx(14,i)*comp(10))*comp(1)	  
	  ENDDO	
	  
c équation dcomp(11)=0.d0	Fe56
c	    jac(11,:)=0.d0
	  
c équation dcomp(12)=0.d0	Ex

c équation dcomp(13)
c	  dcomp(13)=(-(r(16)+r(17))*comp(13)+r(18)*comp(14))*comp(1) !Li6

	  jac(13,1)=-(r(16)+r(17))*comp(13)+r(18)*comp(14)	!d /H1
	  jac(13,13)=-(r(16)+r(17))*comp(1)			!d /Li6
	  jac(13,14)=r(18)*comp(1)				!d /Be9

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(13,i)=jac(13,i)
	1  +(-(drx(16,i)+drx(17,i))*comp(13)+drx(18,i)*comp(14))*comp(1)
	  ENDDO
	  
c équation dcomp(14)
c 	  dcomp(14)=-(r(15)+r(18))*comp(14)*comp(1)	!Be9

	  jac(14,1)=-(r(15)+r(18))*comp(14)		!d /H1
	  jac(14,14)=-(r(15)+r(18))*comp(1)		!d /Be9

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(14,i)=jac(14,i)
	1  -(drx(15,i)+drx(18,i))*comp(14)*comp(1)	  
	  ENDDO	

c	  dcomp(15)=-(r(19)+r(20))*comp(15)*comp(1)	!B11
  
	  jac(15,1)=-(r(19)+r(20))*comp(15)		!d /H1
	  jac(15,15)=-(r(19)+r(20))*comp(1)		!d /B11

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(15,i)=jac(15,i)
	1  -(drx(19,i)+drx(20,i))*comp(15)*comp(1)	  
	  ENDDO	

c conservation des baryons, 12 est l'indice de Ex	 

c	  dcomp(12)=-SUM(anuc*dcomp)/anuc(12)

	  DO j=1,nchim
	   DO i=1,nchim
	    IF(i == 12)CYCLE
	    jac(12,j)=jac(12,j)+anuc(i)*jac(i,j)
	   ENDDO
	   jac(12,j)=-jac(12,j)/anuc(12)	   
	  ENDDO

c unités de temps pour intégration temporelle

	  jac=jac*secon6

	 ENDIF

	 dcomp=dcomp*secon6

c	 calcul de la production d'énergie nucléaire et dérivées

	CASE(3)
	
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0
	 
	 IF(t <= t_inf)RETURN
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	 
c	 mue : nombre d'electrons / mole /g = 1/poids mol. moy. par e-

c H2 et Be7 à l'équilibre
c	 H2, H1(p,e+ nu)H2, H2(p,g)He3, Be9(p,d)2He4

	 dh2h=r(1)/r(2) ; dh2be9=r(15)/r(2)
	 h2=dh2h*comp(1)+dh2be9*comp(14)

	 den=r(6)*mue+r(7)*comp(1)
	 be7=(r(4)*comp(2)*comp(3)+r(17)*comp(1)*comp(13))/den 
	 dbe7he3=r(4)*comp(3)/den ; dbe7he4=r(4)*comp(2)/den
	 dbe7mue=-be7*r(6)/den ; dbe7h=(-be7*r(7)+r(17)*comp(13))/den
	 dbe7li6=r(17)*comp(1)/den

c l'énergie

	 epsilon(2)=(q(1)*comp(1)+q(2)*h2+q(5)*comp(4)+q(7)*be7
	1 +(q(15)+q(18))*comp(14)+(q(16)+q(17))*comp(13)
	2 +(q(19)+q(20))*comp(15))*comp(1)
	3 +(q(3)*comp(2)+q(4)*comp(3))*comp(2)+q(6)*mue*be7
	 epsilon(3)=(q(8)*comp(5)+q(9)*comp(6)+q(10)*comp(7)+
	1 (q(11)+q(12))*comp(8)+q(13)*comp(9)+q(14)*comp(10))*comp(1)
	 DO i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 ENDDO
	 epsilon(1)=SUM(epsilon(2:3))

	 IF(deriv)THEN
	
	  et=(dqt(1)*comp(1)+dqt(2)*h2+dqt(5)*comp(4)+dqt(7)*be7
	1 +(dqt(15)+dqt(18))*comp(14)+(dqt(16)+dqt(17))*comp(13)
	2 +(dqt(19)+dqt(20))*comp(15)+dqt(8)*comp(5)+dqt(9)*comp(6)
	3 +dqt(10)*comp(7)+(dqt(11)+dqt(12))*comp(8)+dqt(13)*comp(9)
	4 +dqt(14)*comp(10))*comp(1)
	5 +(dqt(3)*comp(2)+dqt(4)*comp(3))*comp(2)+dqt(6)*mue*be7
		
	  ero=(dqo(1)*comp(1)+dqo(2)*h2+dqo(5)*comp(4)+dqo(7)*be7
	1 +(dqo(15)+dqo(18))*comp(14)+(dqo(16)+dqo(17))*comp(13)
	2 +(dqo(19)+dqo(20))*comp(15)+dqo(8)*comp(5)+dqo(9)*comp(6)
	3 +dqo(10)*comp(7)+(dqo(11)+dqo(12))*comp(8)+dqo(13)*comp(9)
	4 +dqo(14)*comp(10))*comp(1)
	5 +(dqo(3)*comp(2)+dqo(4)*comp(3))*comp(2)+dqo(6)*mue*be7
		
	  ex(1)=2.d0*q(1)*comp(1)+q(2)*(dh2h+h2*comp(1))+q(5)*comp(4)
	1 +q(8)*comp(5)+q(9)*comp(6)+q(10)*comp(7)
	2 +(q(11)+q(12))*comp(8)+q(13)*comp(9)+q(14)*comp(10)
	3 +q(7)*(be7+dbe7h*comp(1))+q(6)*mue*dbe7h
	4 +(q(15)+q(18))*comp(14)+(q(16)+q(17))*comp(13)
	5 +(q(19)+q(20))*comp(15)	
	  ex(2)=2.d0*q(3)*comp(2)+q(4)*comp(3)+q(7)*dbe7he3*comp(1)
	1 +q(6)*mue*dbe7he3
	  ex(3)=q(4)*comp(2)+q(7)*dbe7he4*comp(1)+q(6)*mue*dbe7he4
	  ex(4)=q(5)*comp(1) ; ex(5)=q(8)*comp(1) ; ex(6)=q(9)*comp(1)
	  ex(7)=q(10)*comp(1) ; ex(8)=(q(11)+q(12))*comp(1)
	  ex(9)=q(13)*comp(1) ; ex(10)=q(14)*comp(1)
	  ex(13)=(q(16)+q(17))*comp(1) ; ex(14)=(q(15)+q(18))*comp(1)
	  ex(15)=(q(19)+q(20))*comp(1)
	  DO i=1,nchim	!contributions des écrans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*h2
	1  +dqx(7,i)*be7+dqx(5,i)*comp(4)+(dqx(15,i)+dqx(18,i))*comp(14)
	2  +(dqx(16,i)+dqx(17,i))*comp(13)+(dqx(19,i)
	3  +dqx(20,i))*comp(15)+dqx(8,i)*comp(5)
	4  +dqx(9,i)*comp(6)+dqx(10,i)*comp(7)
	5  +(dqx(11,i)+dqx(12,i))*comp(8)+dqx(13,i)*comp(9)
	6  +dqx(14,i)*comp(10))*comp(1)
	7  +(dqx(3,i)*comp(2)+dqx(4,i)*comp(3))*comp(2)+dqx(6,i)*mue*be7
	  ENDDO
	 
	 ENDIF	!deriv
	   
c	 production de neutrinos

	CASE(4)

	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

c	  Be7

	  den=(r(6)*mue+r(7)*comp(1)) ; be7=r(4)*comp(2)*comp(3)/den
	
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*be7/amu
	  b8e=r(7)*comp(1)*be7/amu ; n13e=r(8)*comp(1)*comp(5)/amu
	  o15e=r(10)*comp(1)*comp(7)/amu ; f17e=r(13)*comp(1)*comp(9)/amu
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0 ; n13e=0.d0
	  o15e=0.d0 ; f17e=0.d0
	 ENDIF
	 
	CASE default
	 PRINT*,'ppcno10Fe, fait ne peut valoir que 1, 2, 3, 4'
	 PRINT*,'ERREUR fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; STOP
	 
	END SELECT
	
	RETURN

	END SUBROUTINE ppcno10BeBFe
