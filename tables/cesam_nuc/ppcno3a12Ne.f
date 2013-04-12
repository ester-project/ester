
c**************************************************************************

	SUBROUTINE ppcno3a12Ne(t,ro,comp,dcomp,jac,deriv,fait,
	1 epsilon,et,ero,ex,hhe,be7e,b8e,n13e,o15e,f17e)

c routine private du module mod_nuc

c cycles PP, CNO, 3 alpha, Carbone, cf. Clayton p. 380, 392 et 430

c éléments pris en compte:
c H1, He3, He4, C12, C13, N14, N15, O16, O17, O18, Ne20, Ex
c H2, Li7, Be7 à l'équilibre

c un premier appel a rq_reac initialise et définit le nb.
c d'éléments chimiques pour lesquels les reac. nuc. sont tabulées
c dans ppcno9a12Ne on ajoute Ex, soit nchim+1

c Auteurs: P.Morel & B.Pichon Département Cassiopée, O.C.A., CESAM2k

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
c	r(2) : réaction H2(p,g)He3
c	r(3) : réaction He3(He3,2p)He4
c	r(4) : réaction He3(a,g)Be7
c	r(5) : réaction Li7(p,a)He4
c	r(6) : réaction Be7(e-,nu g)Li7
c	r(7) : réaction Be7(p,g)B8(e+ nu)Be8(a)He4

c	r(8) : réaction C12(p,g)N13(e+ nu)C13		CNO
c	r(9) : réaction C13(p,g)N14
c	r(10) : réaction N14(p,g)O15(e+,nu)N15
c	r(11) : réaction N15(p,g)O16
c	r(12) : réaction N15(p,a)C12
c	r(13) : réaction O16(p,g)F17(e+ nu)O17
c	r(14) : réaction O17(p,a)N14

c	r(15) : réaction He4(2a,g)C12		3 alpha
c	r(16) : réaction C12(a,g)O16
c	r(17) : réaction O16(a,g)Ne20

c	r(18) : réaction N14(a,g)F18(e+ nu)O18	(32)
c	r(19) : réaction O18(a,g)Ne22		(33)
c	r(20) : réaction O17(p,g)F18(e+ nu)O18	(24)
c	r(21) : réaction O18(p,a)N15		(28)
c	r(22) : réaction Ne20(a,g)Mg24		(29)

c	r(23) : C12(C12,g)Mg24	(18)		Carbone
c	r(24) : C12(C12,p)Na23	(20)
c	r(25) : C12(C12,a)Ne20	(21)

c indices des isotopes
c	H1 : 1
c	He3 : 2
c	He4 : 3
c	C12 : 4
c	C13 : 5
c	N14 : 6
c	N15 : 7
c	O16 : 8
c	O17 : 9
c	Ne20 :10
c	Ex : 11
c	O18 : 12

c----------------------------------------------------------------------

	USE mod_donnees, ONLY : ab_ini, ab_min, ah, amu, fmin_abon, ihe4, i_ex,
	1 langue, nchim, nom_elem, nom_xheavy, nucleo,
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

c initialisations
	SELECT CASE(fait)
	CASE(0)
	 
c définition de nchim: nombre d'éléments chimiques dont on
c calcule l'abondance H1, He3, He4, C13, C13, N14, N15, O16, O17, O18, Ne20, Ex
	 nchim=12

c appel d'initialisation pour tabulation des réactions nucléaires
	 ALLOCATE(drx(1,1),dqx(1,1),r(1),drt(1),dro(1),q(1),
	1 dqt(1),dqo(1),dmuex(1))
	 CALL rq_reac(comp,1.d7,1.d0,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

	 DEALLOCATE(dqx,drx) ; ALLOCATE(dqx(nreac,nchim),drx(nreac,nchim))
	 
	CASE(1)

c détermination des abondances initiales
c He3+He4=Y0
c Z0 = somme des éléments plus lourds que hélium
c dans Z rapports en nombre
	 CALL abon_ini
	 
c Ex : élément fictif moyenne des éléments utilisés
	 charge_ex=0.d0 ; mass_ex=0.d0 ; sum_a=0.d0
	 B1: DO i=3,nelem_ini		!à partir de Li=3
	  IF(elem(i) == ' C')CYCLE B1
	  IF(elem(i) == ' N')CYCLE B1
	  IF(elem(i) == ' O')CYCLE B1
	  IF(elem(i) == 'Ne')CYCLE B1
	  charge_ex=charge_ex+c(i)*ab(i) ; mass_ex=mass_ex+m(i)*ab(i)
	  sum_a=sum_a+ab(i)
	 ENDDO B1
	 charge_ex=NINT(charge_ex/sum_a) ; mass_ex=NINT(mass_ex/sum_a)	
	 WRITE(text,10)NINT(mass_ex)
10	 FORMAT(i2)

c élément fictif, indice 11
	 i_ex=11		!indice de l'élément chimique reliquat 
	 nucleo(11)=mass_ex	!nucleo de l'élément chimique reliquat
	 zi(11)=charge_ex	!charge de l'élément chimique reliquat
	 i=NINT(charge_ex)
	 nom_elem(11)=elem(i)//text !nom elem. chim. rel.
	 nom_xheavy=nom_elem(11)
 	 SELECT CASE(langue)	  
	 CASE('english')	
	  WRITE(*,1023)TRIM(nom_elem(11)),NINT(mass_ex),NINT(charge_ex)
	  WRITE(2,1023)TRIM(nom_elem(11)),NINT(mass_ex),NINT(charge_ex)	 
1023	  FORMAT(a,': fictitious species /= CNO, of mass : ',i3,/,
	1 'and charge :',i3)	 
	 CASE DEFAULT	 
	  WRITE(*,23)TRIM(nom_elem(11)),NINT(mass_ex),NINT(charge_ex)
	  WRITE(2,23)TRIM(nom_elem(11)),NINT(mass_ex),NINT(charge_ex)	 
23	  FORMAT(a,': élément fictif /= CNO, de masse : ',i3,/,
	1 'et de charge :',i3)
	 END SELECT	 	 
c	 PRINT*,nchim ; WRITE(*,2000)nucleo
	 	 
c détermination des abondances initiales, a(équation,élément)
c cf. Eq. (7,97) p. 135 de la NOTICE
	 ALLOCATE(a(nchim,nchim),indpc(nchim),b(1,nchim))
	 a=0.d0 ; b=0.d0 ; indpc=1	
		
	 a(1,1)=nucleo(1)  	!H1
	 b(1,1)=x0
	
	 a(2,2)=nucleo(2)  	!He3
	 a(2,3)=nucleo(3)	!He4
	 b(1,2)=y0

c abon_rela : abondance relative des métaux dans Z, défini dans abon_ini
	 DO j=4,nchim			!après He4
	  a(3,j)=nucleo(j)		!somme i > 3 comp(i)*nucleo(i)=Z0
	  a(4,j)=-abon_rela(6)		!somme comp(i) C, C/Z
	  a(5,j)=-abon_rela(7)		!somme comp(i) N, N/Z	 
	  a(6,j)=-abon_rela(8)		!somme comp(i) O, O/Z
	  a(11,j)=-abon_rela(10)	!somme comp(i) Ne, Ne/Z	  
	 ENDDO
	 b(1,3)=z0		!Z
	
	 a(4,4)=a(4,4)+1.d0	!C12	 		
	 a(4,5)=a(4,5)+1.d0	!C13
	
	 a(5,6)=a(5,6)+1.d0	!N14	 		
	 a(5,7)=a(5,7)+1.d0	!N15
	
	 a(6,8)=a(6,8)+1.d0	!O16	 		
	 a(6,9)=a(6,9)+1.d0	!O17
	 a(6,12)=a(6,12)+1.d0	!O18	 
	 
	 a(11,10)=a(11,10)+1.d0	!Ne20
	
c rapports isotopiques
	 a(7,2)=1.d0		!He3
	 a(7,3)=-he3she4z	!He3/He4, H2 est dans He3		 		

	 a(8,5)=1.d0		!C13
	 a(8,4)=-c13sc12	!C13/C12
	
	 a(9,7)=1.d0		!N15
	 a(9,6)=-n15sn14	!N15/N14
	
	 a(10,9)=1.d0		!O17
	 a(10,8)=-o17so16	!O17/O16
	 
	 a(12,12)=1.d0		!O18
	 a(12,8)=-o18so16	!O18/O16

c	 PRINT*,nchim
c	 DO i=1,nchim
c	  WRITE(*,2002)a(i,1:nchim),b(1,i)
c	 ENDDO

	 CALL gauss_band(a,b,indpc,nchim,nchim,nchim,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'ppcno3a12Ne, matrice du calcul des abon. non inversible'
	  PRINT*,'ARRET' ; STOP
	 ENDIF
c	 WRITE(*,2002)b(1,:) ; PAUSE'ab.ini'

c allocations diverses
	 DEALLOCATE(drt,dro,r,q,dqt,dqo,dmuex)
	 ALLOCATE(ab_ini(nchim),ab_min(nchim),drt(nreac),dro(nreac),
	1 r(nreac),q(nreac),dqt(nreac),dqo(nreac),anuc(nchim),
	2 dmuex(nchim),dh2x(nchim),denx(nchim),dbe7x(nchim),dli7x(nchim))

c abondances initiales et abondances négligeables
	 comp(1:nchim)=MAX(1.d-29,b(1,1:nchim))
	 ab_ini=comp*nucleo	 	
	 ab_min=ab_ini*fmin_abon ; nbz=SUM(comp(ihe4+1:nchim))
	 	 
c abondances en DeX, H=12
	 ALLOCATE(comp_dex(nchim))
	 comp_dex=12.d0+LOG10(comp/comp(1))	 

c écritures		
	 WRITE(2,2) ; WRITE(*,2) 
2	 FORMAT(/,'Réac. thermonucléaires cycles PP, CNO, 3 alpha, Carbone, T < 1GK',/)
	 WRITE(2,3)nreac ; WRITE(*,3)nreac 
3	 FORMAT('nombre de réactions : ',i3)
	 WRITE(2,4)nreac ; WRITE(*,4)nchim
4	 FORMAT('nombre d''éléments chimiques : ',i3)
	 WRITE(2,20)x0,y0,z0,z0/x0 ; WRITE(*,20)x0,y0,z0,z0/x0	
20	 FORMAT(/,'abondances initiales déduites de:',/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3,/,'Z0/X0=',es10.3,
	2 ', H1=X0, He3+He4=Y0',/,
	3 'Z0 = 1-X0-Y0 = C12+C13+N14+N15+O16+O17+O18+Ne20+Ex')
	 WRITE(2,1)ab_ini ; WRITE(*,1)ab_ini
1	 FORMAT(/,'abondances initiales/gramme:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', C12:',es10.3,
	2 ', C13:',es10.3,/'N14:',es10.3,', N15:',es10.3,', O16:',es10.3,
	3 ', O17:',es10.3,', Ne20:',es10.3/,'Ex:',es10.3,', O18:',es10.3)
	 WRITE(2,9)comp_dex ; WRITE(*,9)comp_dex
9	 FORMAT(/,'Abondances initiales en nombre: 12+Log10(Ni/Nh)',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', C12:',es10.3,
	2 ', C13:',es10.3,/'N14:',es10.3,', N15:',es10.3,', O16:',es10.3,
	3 ', O17:',es10.3,', Ne20:',es10.3/,'Ex:',es10.3,', O18:',es10.3)		
	 WRITE(2,21)(comp(4)+comp(5))/nbz,	!C/Z
	4 (comp(6)+comp(7))/nbz,		!N/Z
	5 (comp(8)+comp(9)+comp(12))/nbz,	!O/Z
	5 comp(10)/nbz,				!Ne/Z	
	6 comp(11)/nbz		!Ex/Z
	 WRITE(*,21)(comp(4)+comp(5))/nbz,	!C/Z
	4 (comp(6)+comp(7))/nbz,		!N/Z
	5 (comp(8)+comp(9)+comp(12))/nbz,	!O/Z
	5 comp(10)/nbz,				!Ne/Z	
	6 comp(11)/nbz				!Ex/Z
21	 FORMAT(/,'Rapports en nombre dans Z:',/,'C/Z:',es10.3,', N/Z:',es10.3,
	1 ', O/Z:',es10.3,', Ne/Z:',es10.3,', Ex/Z:',es10.3)
	 WRITE(2,14)he3she4,c13sc12,n15sn14,o17so16,o18so16
	 WRITE(*,14)he3she4,c13sc12,n15sn14,o17so16,o18so16	  
14	 FORMAT(/,'Rapports isotopiques en nombre:',/,
	1 'HE3/HE4=',es10.3,', C13/C12=',es10.3,
	2 ', N15/N14=',es10.3,/,'O17/O16=',es10.3,', O18/O16=',es10.3)		
	 WRITE(2,5)ab_min ; WRITE(*,5)ab_min
5	 FORMAT(/,'abondances/gramme négligeables:',/,
	1 'H1:',es10.3,', He3:',es10.3,', He4:',es10.3,', C12:',es10.3,
	2 ', C13:',es10.3,/'N14:',es10.3,', N15:',es10.3,', O16:',es10.3,
	3 ', O17:',es10.3,', Ne20:',es10.3/,'Ex:',es10.3,', O18:',es10.3)	
	 WRITE(2,6) ; WRITE(*,6)
6	 FORMAT(/,'H2, Li7, Be7 à l''équilibre')
	 WRITE(2,7) ; WRITE(*,7)
7	 FORMAT(/,'on utilise une table')	 
	 WRITE(2,8) ; WRITE(*,8)
8	 FORMAT(/,'évol. temporelle, test de précision sur H1 et He4')

c définitions diverses
	  ab_min=ab_min/nucleo ; anuc=ANINT(nucleo)	!nombre atomique
	 
c nettoyage	
	 DEALLOCATE(a,b,comp_dex,indpc)	 
	 
c les réactions	  
	CASE(2)
	 dcomp=0.d0 ; jac=0.d0
	
	 IF(t < t_inf)RETURN	 
c	 WRITE(*,2000)t,ro ; PRINT*,'entrée'	
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
c	 WRITE(*,*)'comp' ; WRITE(*,2000)comp(1:nchim)
c	 WRITE(*,*)'réactions' ; WRITE(*,2000)r(1:nreac)

c équations d'évolution
	 dcomp(1)=-(3.d0*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)
	1 +r(10)*comp(6)+(r(11)+r(12))*comp(7)+r(13)*comp(8)
	2 +(r(14)+r(20))*comp(9)+r(21)*comp(12))*comp(1)
	3 +(2.d0*r(3)*comp(2)-r(4)*comp(3))*comp(2)+r(24)*comp(4)**2	!H1

	 dcomp(2)=r(1)*comp(1)**2-(2.d0*r(3)*comp(2)
	1 +r(4)*comp(3))*comp(2)				!He3

	 dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
	1 +(r(12)*comp(7)+r(14)*comp(9)+r(21)*comp(12))*comp(1)
	2 -(3.d0*r(15)*comp(3)**2+r(16)*comp(4)+r(17)*comp(8)
	3 +r(18)*comp(6)+r(19)*comp(12)+r(22)*comp(10))*comp(3)
	4 +r(25)*comp(4)**2					!He4

	 dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)
	1 +(r(15)*comp(3)**2-r(16)*comp(4))*comp(3)
	2 -(r(23)+r(24)+r(25))*comp(4)**2			!C12

	 dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)		!C13

	 dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)
	1 +r(14)*comp(9))*comp(1)-r(18)*comp(3)*comp(6)		!N14

	 dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7)
	1 +r(21)*comp(12))*comp(1)				!N15

	 dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)
	1 +(r(16)*comp(4)-r(17)*comp(8))*comp(3)		!O16

	 dcomp(9)=(r(13)*comp(8)-(r(14)+r(20))*comp(9))*comp(1)	!O17

	 dcomp(10)=(r(17)*comp(8)-r(22)*comp(10))*comp(3)
	1 +r(25)*comp(4)**2					!Ne20

	 dcomp(12)=r(18)*comp(3)*comp(6)+r(20)*comp(1)*comp(9)	!O18
	1 -(r(19)*comp(3)+r(21)*comp(1))*comp(12)	
	
c conservation des baryons
	 dcomp(11)=-DOT_PRODUCT(dcomp,anuc)/anuc(11)

c calcul du jacobien
	 IF(deriv)THEN	!jac(i,j) : équation, j : élément i
	
c équation dcomp(1), H1
c	  dcomp(1)=-(3.d0*r(1)*comp(1)+r(8)*comp(4)+r(9)*comp(5)
c	1 +r(10)*comp(6)+(r(11)+r(12))*comp(7)+r(13)*comp(8)
c	2 +(r(14)+r(20))*comp(9)+r(21)*comp(12))*comp(1)
c	3 +(2.d0*r(3)*comp(2)-r(4)*comp(3))*comp(2)+r(24)*comp(4)**2

	  jac(1,1)=-6.d0*r(1)*comp(1)-r(8)*comp(4)-r(9)*comp(5)
	1 -r(10)*comp(6)-(r(11)+r(12))*comp(7)-r(13)*comp(8)
	2 -(r(14)+r(20))*comp(9)-r(21)*comp(12)			!d /H1
	  jac(1,2)=4.d0*r(3)*comp(2)-r(4)*comp(3)		!d /He3
	  jac(1,3)=-r(4)*comp(2)				!d /He4
	  jac(1,4)=-r(8)*comp(1)+2.d0*r(24)*comp(4)		!d /C12
	  jac(1,5)=-r(9)*comp(1)				!d /C13
	  jac(1,6)=-r(10)*comp(1)				!d /N14
	  jac(1,7)=-(r(11)+r(12))*comp(1)			!d /N15
	  jac(1,8)=-r(13)*comp(1)				!d /O16
	  jac(1,9)=-(r(14)+r(20))*comp(1)			!d /O17
	  jac(1,12)=-r(21)*comp(1)				!d /O18
	  
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(1,i)=jac(1,i) 
	1  -(3.d0*drx(1,i)*comp(1)+drx(8,i)*comp(4)
	2  +drx(9,i)*comp(5)+drx(10,i)*comp(6)
	3  +(drx(11,i)+drx(12,i))*comp(7)+drx(13,i)*comp(8)
	4  +(drx(14,i)+drx(20,i))*comp(9)+drx(21,i)*comp(12))*comp(1)
	5  +(2.d0*drx(3,i)*comp(2)-drx(4,i)*comp(3))*comp(2)
	6  +drx(24,i)*comp(4)**2
	  ENDDO
			 
c équation dcomp(2), He3
c	  dcomp(2)=r(1)*comp(1)**2-(2.d0*r(3)*comp(2)+r(4)*comp(3))*comp(2)

	  jac(2,1)=2.d0*r(1)*comp(1)				!d /H1
	  jac(2,2)=-4.d0*r(3)*comp(2)-r(4)*comp(3)		!d /He3
	  jac(2,3)=-r(4)*comp(2)				!d /He4

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(2,i)=jac(2,i)
	1  +drx(1,i)*comp(1)**2-(2.d0*drx(3,i)*comp(2)+drx(4,i)*comp(3))*comp(2)
	  ENDDO
	 
c équation dcomp(3) He4
c	  dcomp(3)=(r(3)*comp(2)+r(4)*comp(3))*comp(2)
c	1 +(r(12)*comp(7)+r(14)*comp(9)+r(21)*comp(12))*comp(1)
c	2 -(3.d0*r(15)*comp(3)**2+r(16)*comp(4)+r(17)*comp(8)
c	3 +r(18)*comp(6)+r(19)*comp(12)+r(22)*comp(10))*comp(3)
c	4 +r(25)*comp(4)**2
	  jac(3,1)=r(12)*comp(7)+r(14)*comp(9)+r(21)*comp(12)	!d /H1
	  jac(3,2)=2.d0*r(3)*comp(2)+r(4)*comp(3)		!d /He3
	  jac(3,3)=r(4)*comp(2)-9.d0*r(15)*comp(3)**2-r(16)*comp(4)
	1 -r(17)*comp(8)-r(18)*comp(6)
	2 -r(19)*comp(12)-r(22)*comp(10)			!d /He4
	  jac(3,4)=-r(16)*comp(3)+2.d0*r(25)*comp(4)		!d /C12
	  jac(3,6)=-r(18)*comp(3)				!d /N14	  
	  jac(3,7)=r(12)*comp(1)				!d /N15
	  jac(3,8)=-r(17)*comp(3)				!d /O16	 
	  jac(3,9)=r(14)*comp(1)				!d /O17
	  jac(3,10)=-r(22)*comp(3)				!d /Ne20
	  jac(3,12)=r(21)*comp(1)-r(19)*comp(3)			!d /O18

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(3,i)=jac(3,i)
	1  +(drx(3,i)*comp(2)+drx(4,i)*comp(3))*comp(2)
	2  +(drx(12,i)*comp(7)+drx(14,i)*comp(9)+drx(21,i)*comp(12))*comp(1)
	3  -(3.d0*drx(15,i)*comp(3)**2
	4  +drx(16,i)*comp(4)+drx(17,i)*comp(8)+drx(18,i)*comp(6)
	5  +drx(19,i)*comp(12)+drx(22,i)*comp(10))*comp(3)+drx(25,i)*comp(4)**2
	  ENDDO
	 
c équation dcomp(4), C12
c	  dcomp(4)=(-r(8)*comp(4)+r(12)*comp(7))*comp(1)
c	1 +(r(15)*comp(3)**2-r(16)*comp(4))*comp(3)
c	2 -(r(23)+r(24)+r(25))*comp(4)**2
	  jac(4,1)=-r(8)*comp(4)+r(12)*comp(7)			!d /H1
	  jac(4,3)=3.d0*r(15)*comp(3)**2-r(16)*comp(4)		!d /He4	 
	  jac(4,4)=-r(8)*comp(1)-r(16)*comp(3)
	1 -2.d0*(r(23)+r(24)+r(25))*comp(4)			!d /C12
	  jac(4,7)=r(12)*comp(1)				!d /N15
	 		 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(4,i)=jac(4,i)
	1  +(-drx(8,i)*comp(4)+drx(12,i)*comp(7))*comp(1)
	2  +(drx(15,i)*comp(3)**2-drx(16,i)*comp(4))*comp(3)
	3 -(drx(23,i)+drx(24,i)+drx(25,i))*comp(4)**2
	  ENDDO
	 	 
c équation dcomp(5), C13
c	  dcomp(5)=(r(8)*comp(4)-r(9)*comp(5))*comp(1)
	  	
	  jac(5,1)=r(8)*comp(4)-r(9)*comp(5)			!d /H1
	  jac(5,4)=r(8)*comp(1)					!d /C12
	  jac(5,5)=-r(9)*comp(1)				!d /C13

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(5,i)=jac(5,i)
	1  +(drx(8,i)*comp(4)-drx(9,i)*comp(5))*comp(1)	 
	  ENDDO
	
c équation dcomp(6), N14
c	  dcomp(6)=(r(9)*comp(5)-r(10)*comp(6)
c	1 +r(14)*comp(9))*comp(1)-r(18)*comp(3)*comp(6)

	  jac(6,1)=r(9)*comp(5)-r(10)*comp(6)+r(14)*comp(9)	!d /H1
	  jac(6,3)=-r(18)*comp(6)				!d /He4
	  jac(6,5)=r(9)*comp(1)					!d /C13
	  jac(6,6)=-r(10)*comp(1)-r(18)*comp(3)			!d /N14
	  jac(6,9)=r(14)*comp(1)				!d /O17

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(6,i)=jac(6,i)
	1  +(drx(9,i)*comp(5)-drx(10,i)*comp(6)
	2  +drx(14,i)*comp(9))*comp(1)-drx(18,i)*comp(3)*comp(6)
	  ENDDO
	 	 
c équation dcomp(7), N15
c	  dcomp(7)=(r(10)*comp(6)-(r(11)+r(12))*comp(7)
c	1 +r(21)*comp(12))*comp(1)
	
	  jac(7,1)=r(10)*comp(6)-(r(11)+r(12))*comp(7)
	1 +r(21)*comp(12)					!d /H1
	  jac(7,6)=r(10)*comp(1)				!d /N14
	  jac(7,7)=-(r(11)+r(12))*comp(1)			!d /N15
	  jac(7,12)=r(21)*comp(1)				!d /O18
	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(7,i)=jac(7,i)
	1  +(drx(10,i)*comp(6)-(drx(11,i)+drx(12,i))*comp(7)
	2  +drx(21,i)*comp(12))*comp(1)
	  ENDDO
	 	 
c équation dcomp(8), O16
c	  dcomp(8)=(r(11)*comp(7)-r(13)*comp(8))*comp(1)
c	1 +(r(16)*comp(4)-r(17)*comp(8))*comp(3)

	  jac(8,1)=r(11)*comp(7)-r(13)*comp(8)		!d /H1
	  jac(8,3)=r(16)*comp(4)-r(17)*comp(8)		!d /He4
	  jac(8,4)=r(16)*comp(3)			!d /C12	 	 
	  jac(8,7)=r(11)*comp(1)			!d /N15
	  jac(8,8)=-r(13)*comp(1)-r(17)*comp(3)		!d /O16
	 	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(8,i)=jac(8,i)
	1  +(drx(11,i)*comp(7)-drx(13,i)*comp(8))*comp(1)
	2  +(drx(16,i)*comp(4)-drx(17,i)*comp(8))*comp(3)
	  ENDDO
	 
c équation dcomp(9), O17
c	  dcomp(9)=(r(13)*comp(8)-(r(14)+r(20))*comp(9))*comp(1)

	  jac(9,1)=r(13)*comp(8)-(r(14)+r(20))*comp(9)	!d /H1
	  jac(9,8)=r(13)*comp(1)			!d /O16
	  jac(9,9)=-(r(14)+r(20))*comp(1)		!d /O17
	 
	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(9,i)=jac(9,i)
	1  +(drx(13,i)*comp(8)-(drx(14,i)+drx(20,i))*comp(9))*comp(1)
	  ENDDO
	 
c équation dcomp(10), Ne20	
c	  dcomp(10)=(r(17)*comp(8)-r(22)*comp(10))*comp(3)
c	1 +r(25)*comp(4)**2

	  jac(10,3)=r(17)*comp(8)-r(22)*comp(10)	!d /He4
	  jac(10,4)=2.d0*r(25)*comp(4)			!d C12	  
	  jac(10,8)=r(17)*comp(3)			!d /O16
	  jac(10,10)=-r(22)*comp(3)			!d /Ne20 

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(10,i)=jac(10,i)
	1  +(drx(17,i)*comp(8)-drx(22,i)*comp(10))*comp(3)+drx(25,i)*comp(4)**2
	  ENDDO

c équation dcomp(12), O18	  
c	  dcomp(12)=r(18)*comp(3)*comp(6)+r(20)*comp(1)*comp(9)
c	1 -(r(19)*comp(3)+r(21)*comp(1))*comp(12)

	  jac(12,1)=r(20)*comp(9)-r(21)*comp(12)	!d /H1	  
	  jac(12,3)=r(18)*comp(6)-r(19)*comp(12)	!d /He4
	  jac(12,6)=r(18)*comp(3)			!d /N14
	  jac(12,9)=r(20)*comp(1)			!d /O17
	  jac(12,12)=-r(19)*comp(3)-r(21)*comp(1)	!d /O18 

	  DO i=1,nchim		!dépendances dues à l'effet d'écran
	   jac(12,i)=jac(12,i)+drx(18,i)*comp(3)*comp(6)
	1  +drx(20,i)*comp(1)*comp(9)-(drx(19,i)*comp(3)
	2  +drx(21,i)*comp(1))*comp(12)
	  ENDDO
	  
c équation dcomp(11), conservation des baryons	 
c 	  dcomp(11)=-SUM(dcomp*anuc)/anuc(11) 	conservation des baryons	
	  DO j=1,nchim
	   B2: DO i=1,nchim
	    IF(i == 11)CYCLE B2
	    jac(11,j)=jac(11,j)+anuc(i)*jac(i,j)
	   ENDDO B2
	   jac(11,j)=-jac(11,j)/anuc(11)
	  ENDDO	  
	 
c unités de temps pour intégration temporelle
	  jac=jac*secon6

	 ENDIF		!deriv
	 
	 dcomp=dcomp*secon6

	CASE(3)
	
c	 calcul de la production d'énergie nucléaire et dérivées
c	 pour H2(H,g)He3, q(2)H**2=q(2)*r(1)/r(2)
	 
	 epsilon(1:4)=0.d0 ; et=0.d0 ; ero=0.d0 ; ex=0.d0
	 IF(t <= t_inf)return
	
	 CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

c	 mue : nombre d'electrons / mole /g = 1/poids mol. moy. par e-

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
	 epsilon(4)=(q(15)*comp(3)**2+q(16)*comp(4)+q(17)*comp(8)
	1 +q(18)*comp(6)+q(19)*comp(12)+q(22)*comp(10))*comp(3)
	2 +(q(20)*comp(9)+q(21)*comp(12))*comp(1)
	3 +(q(23)+q(24)+q(25))*comp(4)**2
	 epsilon(1)=SUM(epsilon(2:4))
	 	
	 IF(deriv)THEN
	
	  IF(h2 > 0.d0)THEN
	   dh2t=h2*(drt(1)/r(1)-drt(2)/r(2)) ; dh2ro=h2*(dro(1)/r(1)-dro(2)/r(2))
	   DO i=1,nchim
	    dh2x(i)=h2*(drx(1,i)/r(1)-drx(2,i)/r(2))
	   ENDDO
	   dh2x(1)=dh2x(1)+h2/comp(1)	  
	  ELSE
	   dh2t=0.d0 ; dh2ro=0.d0
	   DO i=1,nchim	  
	    dh2x(i)=0.d0
	   ENDDO
	  ENDIF

	  IF(be7 > 0.d0)THEN
	   dent= drt(6)*mue+drt(7)*comp(1) 
	   denro=dro(6)*mue+dro(7)*comp(1)
	   DO i=1,nchim	  
	    denx(i)=drx(6,i)*mue+r(6)*dmuex(i)+drx(7,i)*comp(1)
	   ENDDO
	   denx(1)=denx(1)+r(7)
	   
	   dbe7t= be7*(drt(4)/r(4)- dent/den)
	   dbe7ro=be7*(dro(4)/r(4)-denro/den)
	   DO i=1,nchim
	    dbe7x(i)=be7*(drx(4,i)/r(4)-denx(i)/den)
	   ENDDO
	   dbe7x(2)=dbe7x(2)+be7/comp(2); dbe7x(3)=dbe7x(3)+be7/comp(3)	   
	  ELSE
	   dbe7t=0.d0 ; dbe7ro=0.d0
	   DO i=1,nchim	   
	    dbe7x(i)=0.d0
	   ENDDO
	  ENDIF

	  IF(li7 > 0.d0)THEN
	   dli7t= li7*(drt(6)/r(6) +dbe7t/be7-drt(5)/r(5))
	   dli7ro=li7*(dro(6)/r(6)+dbe7ro/be7-dro(5)/r(5))
	   DO i=1,nchim	   
	    dli7x(i)=li7*(drx(6,i)/r(6)+dbe7x(i)/be7
	1   +dmuex(i)/mue-drx(5,i)/r(5))
	   ENDDO
	   dli7x(1)=dli7x(1)-li7/comp(1)	   
	  ELSE
	   dli7t=0.d0 ; dli7ro=0.d0
	   DO i=1,nchim
	    dli7x(i)=0.d0
	   ENDDO
	  ENDIF
	  	
 	  et=(dqt(1)*comp(1)+dqt(2)*h2+q(2)*dh2t+dqt(5)*li7+q(5)*dli7t
	1 +dqt(7)*be7+q(7)*dbe7t+dqt(8)*comp(4)+dqt(9)*comp(5)+dqt(10)*comp(6)
	2 +(dqt(11)+dqt(12))*comp(7)+dqt(13)*comp(8)+dqt(14)*comp(9)
	3 + dqt(20)*comp(9)+dqt(21)*comp(12))*comp(1)
	5 +(dqt(3)*comp(2)+dqt(4)*comp(3))*comp(2)
	6 +dqt(6)*mue*be7+q(6)*mue*dbe7t
	7 +(dqt(15)*comp(3)**2+dqt(16)*comp(4)
	8 +dqt(17)*comp(8)+dqt(18)*comp(6)+dqt(19)*comp(12)
	9 +dqt(22)*comp(10))*comp(3)+(dqt(23)+dqt(24)+dqt(25))*comp(4)**2

 	  ero=(dqo(1)*comp(1)+dqo(2)*h2+q(2)*dh2ro+dqo(5)*li7+q(5)*dli7ro
	1 +dqo(7)*be7+q(7)*dbe7ro+dqo(8)*comp(4)+dqo(9)*comp(5)+dqo(10)*comp(6)
	2 +(dqo(11)+dqo(12))*comp(7)+dqo(13)*comp(8)+dqo(14)*comp(9)
	3 +dqo(20)*comp(9)+dqo(21)*comp(12))*comp(1)
	5 +(dqo(3)*comp(2)+dqo(4)*comp(3))*comp(2)
	6 +dqo(6)*mue*be7+q(6)*mue*dbe7ro
	7 +(dqo(15)*comp(3)**2+dqo(16)*comp(4)
	8 +dqo(17)*comp(8)+dqo(18)*comp(6)+dqo(19)*comp(12)
	9 +dqo(22)*comp(10))*comp(3)+(dqo(23)+dqo(24)+dqo(25))*comp(4)**2
	
	  ex(1)=2.d0*q(1)*comp(1)+q(2)*h2+q(5)*li7+q(7)*be7
	1 +q(8)*comp(4)+q(9)*comp(5)+q(10)*comp(6)+(q(11)+q(12))*comp(7)
	2 +q(13)*comp(8)+q(14)*comp(9)+q(20)*comp(9)+q(21)*comp(12)
	  ex(2)=2.d0*q(3)*comp(2)+q(4)*comp(3)
	  ex(3)=q(4)*comp(2)+3.d0*q(15)*comp(3)**2+q(16)*comp(4)
	1 +q(17)*comp(8)+q(18)*comp(6)+q(19)*comp(12)+q(22)*comp(10)
	  ex(4)=q(8)*comp(1)+q(16)*comp(3)+2.d0*(q(23)+q(24)+q(25))*comp(4)
	  ex(5)=q(9)*comp(1)
	  ex(6)=q(10)*comp(1)+q(18)*comp(3)
	  ex(7)=(q(11)+q(12))*comp(1)
	  ex(8)=q(13)*comp(1)+q(17)*comp(3)
	  ex(9)=(q(14)+q(20))*comp(1)
	  ex(10)=q(22)*comp(3)
	  ex(12)=q(21)*comp(1)+q(19)*comp(3)
	 
	  DO i=1,nchim	!contributions des écrans
	   ex(i)=ex(i)+(dqx(1,i)*comp(1)+dqx(2,i)*h2
	1  +q(2)*dh2x(i)+dqx(5,i)*li7+q(5)*dli7x(i)+dqx(7,i)*be7+q(7)*dbe7x(i)
	2  +dqx(8,i)*comp(4)+dqx(9,i)*comp(5)+dqx(10,i)*comp(6)
	3  +(dqx(11,i)+dqx(12,i))*comp(7)+dqx(13,i)*comp(8)+dqx(14,i)*comp(9)
	4  +dqx(20,i)*comp(9)+dqx(21,i)*comp(12))*comp(1)
	5  +(dqx(3,i)*comp(2)+dqx(4,i)*comp(3))*comp(2)
	6  +dqx(6,i)*mue*be7+q(6)*dmuex(i)*be7+q(6)*mue*dbe7x(i)
	7  +(dqx(15,i)*comp(3)**2+dqx(16,i)*comp(4)+dqx(17,i)*comp(8)
	8  +dqx(18,i)*comp(6)+dqx(19,i)*comp(12)+dqx(22,i)*comp(10))*comp(3)
	9 +(dqx(23,i)+dqx(24,i)+dqx(25,i))*comp(4)**2

	 ENDDO
	 
	 ENDIF	!deriv

	CASE(4)		!taux de production des neutrinos

	 IF(t >= t_inf)THEN
	  CALL rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)
	  be7=r(4)*comp(2)*comp(3)/(r(6)*mue+r(7)*comp(1))
	  hhe=r(1)*comp(1)**2/amu ; be7e=r(6)*mue*be7/amu
	  b8e=r(7)*comp(1)*be7/amu ; n13e=r(8)*comp(1)*comp(4)/amu
	  o15e=r(10)*comp(1)*comp(6)/amu
	  f17e=r(13)*comp(1)*comp(8)/amu
	 ELSE
	  hhe=0.d0 ; be7e=0.d0 ; b8e=0.d0 ; n13e=0.d0
	  o15e=0.d0 ; f17e=0.d0
	 ENDIF

	CASE default
	 PRINT*,'ppcno3a12Ne, fait ne peut valoir que 1, 2, 3 ou 4'
	 PRINT*,'ERREUR fait a la valeur:',fait
	 PRINT*,'ARRET' ; PRINT* ; STOP
	 
	END SELECT
	 	
	RETURN

	END SUBROUTINE ppcno3a12Ne
