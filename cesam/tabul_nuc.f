
c*************************************************************************

	SUBROUTINE tabul_nuc(ar,i3al,knot_temp,m_temp,
	1 nom_reac,n_temp,q0,taux_reac,temp,ttemp)
	
c	routine private du module mod_nuc

c	tabulation des r�actions thermonucl�aires (d�riv� du programme
c	tab_reac) le nom de la routine de r�actions nucl�aires
c	(par Ex. pour ppcno9) d�termine les r�actions, les �l�ments �
c	utiliser l'intervalle de temp�ratures pour la tabulation

c	dans un premier appel � taux_nuc on initialise les �nergies,
c	les masses r�duites etc... de toutes les nreac_tot r�ac. nucl�aires
c	consid�r�es, puis on ne garde que celles � utiliser

c	des appels � taux_reac_nuc calculent les taux, on tabule en f(ln T)

c	les variables avec terminaison `t' telles que nucleot, zit etc..
c	sont dimensionn�es � nreac_tot >= nreac car toutes les r�actions
c	ne sont pas toujours toutes utilis�es, idem pour les �l�ments
c	chimiques

c	voir dans taux_reac la proc�dure pour ajouter des �l�ments, des
c	isotopes, des r�actions

c	Auteur: P.Morel, D�partement J.D. Cassini, O.C.A.
c	CESAM2k

c--------------------------------------------------------------------

	USE mod_donnees, ONLY : amu, eve, ihe4, langue, nchim,
	1 nom_elem, nom_nuc, nom_nuc_cpl, nucleo, t_inf, t_stop, zi
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(out), ALLOCATABLE, DIMENSION(:,:) ::
	1 taux_reac
	REAL (kind=dp), INTENT(out), ALLOCATABLE, DIMENSION(:) :: ar,
	1 q0, temp, ttemp
	INTEGER, INTENT(out) :: i3al, knot_temp, m_temp, n_temp
	CHARACTER (len=20), INTENT(out), ALLOCATABLE, DIMENSION(:) ::
	1 nom_reac
	
	REAL (kind=dp), DIMENSION(nreac_tot) :: rt, qt, at
	REAL (kind=dp), DIMENSION(-1:niso_tot) :: nucleot, zit	
	REAL (kind=dp) :: ti, t0, t1, pas

	INTEGER, DIMENSION(nreac_tot,2) :: izzt
	INTEGER, DIMENSION(nreac_tot) :: nb
	INTEGER, DIMENSION(-1:niso_tot) :: ind
	INTEGER :: total, i, j, nchim_reac
					
	CHARACTER (len=4), DIMENSION(-1:niso_tot) :: nom_elemt 
	CHARACTER (len=20), DIMENSION(nreac_tot) ::  nom_react

c-----------------------------------------------------------------	 

2000	FORMAT(8es10.3)
2002	FORMAT(20i4)

c	nchim_reac est le nombre d'isotopes utilis�s dans les r�actions
c	nchim=nchim_reac + Ex, nchim est d�fini dans la routine de
c	r�actions nucl�aires, Ex: ppcno9
c	PRINT*,nchim ; PAUSE'entr�e tabul'

	SELECT CASE (langue)
	CASE('english')
	 WRITE(2,1001)nom_nuc ; WRITE(*,1001)nom_nuc
1001	 FORMAT('Thermonuclear reactions tabulated for : ',a20)	
	CASE DEFAULT
	 WRITE(2,1)nom_nuc ; WRITE(*,1)nom_nuc
1	 FORMAT('R�actions nucl�aires tabul�es pour: ',a20)
	END SELECT 
	 
	SELECT CASE(nom_nuc)
	CASE ('pp1')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1002) ; WRITE(*,1002) 	
1002	  FORMAT('No tabulation for the thermonuclear reaction PP1',/)
	 CASE DEFAULT		
	  WRITE(2,2) ; WRITE(*,2) 	
2	  FORMAT('Pas de tabulation de la r�action PP1 simplifi�',/)
	 END SELECT
	 
	CASE ('ppcno9','ppcno9Fe')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1003) ; WRITE(*,1003) 
1003	  FORMAT('PP+CNO, 9 species, H2, Li7, Be7 at equilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,3) ; WRITE(*,3) 
3	  FORMAT('PP+CNO, 9 �l�ments, H2, Li7, Be7 � l''�quilibre',/)
	 END SELECT
	 	  
c nombre de r�actions et correspondance des indices
	 nreac=14	!nombre de r�actions	 
	 nchim_reac=9
	 ihe4=3	!indice de He4
	 i3al=0	!indice de la r�action 3 alpha

c indices des isotopes utilis�s dans les r�actions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
 	 
c ordre des r�actions: on garde l'ordre de reac_c 
 	 nb(1:nreac)=(/ (i, i=1,nreac) /)
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201
	 
	CASE ('pp3')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1021) ; WRITE(*,1021) 
1021	  FORMAT('PP, 3 species, H2, Li7, Be7 at equilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,21) ; WRITE(*,21) 
21	  FORMAT('PP, 3 �l�ments, H2, Li7, Be7 � l''�quilibre',/)
	 END SELECT
	 	  
c nombre de r�actions et correspondance des indices
	 nreac=7	!nombre de r�actions	 
	 nchim_reac=3
	 ihe4=3	!indice de He4
	 i3al=0	!indice de la r�action 3 alpha

c indices des isotopes utilis�s dans les r�actions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 
c ordre des r�actions: on garde l'ordre de reac_c 
 	 nb(1:nreac)=(/ (i, i=1,nreac) /)
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=80.d6 ; n_temp=201
	  
	CASE('ppcno10','ppcno10Fe','ppcno10K')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1004) ; WRITE(*,1004) 	
1004	  FORMAT('PP+CNO, 10 species, H2, Be7 at equilibrium',/)
	 CASE DEFAULT		
	  WRITE(2,4) ; WRITE(*,4) 	
4	  FORMAT('PP+CNO, 10 �l�ments, H2, Be7 � l''�quilibre',/)
	 END SELECT
	 	 
c nombre de r�actions et correspondance des indices 
	 nreac=14	!nombre de r�actions
	 nchim_reac=10	 
	 ihe4=3		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
	 	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=3	!He3 
 	 ind(3)=4	!He4
	 ind(4)=5	!Li7 
	 ind(5)=7	!C12
	 ind(6)=8	!C13	 
	 ind(7)=10	!N14	 
 	 ind(8)=11	!N15
 	 ind(9)=12	!O16
 	 ind(10)=13	!O17

c ordre des r�actions: on garde l'ordre de reac_c
	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201
	 
	CASE('ppcno10BeBFe')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1016) ; WRITE(*,1016) 	
1016	  FORMAT('PP + CNO + Be + B + Fe with H2 and Be7 at equilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,16) ; WRITE(*,16) 	
16	  FORMAT('PP + CNO + Be + B + Fe avec H2, Be7 � l''�quilibre',/)
	 END SELECT
	 
c nombre de r�actions et correspondance des indices 
	 nreac=20	!nombre de r�actions
	 nchim_reac=15	 
	 ihe4=3		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=3	!He3 
 	 ind(3)=4	!He4
	 ind(4)=5	!Li7
	 ind(5)=7	!C12
	 ind(6)=8	!C13	 
	 ind(7)=10	!N14	 
 	 ind(8)=11	!N15
 	 ind(9)=12	!O16
 	 ind(10)=13  	!O17
 	 ind(11)=27	!Fe56
 	 ind(12)=-1	!Ex
 	 ind(13)=26	!Li6	 	 	 
 	 ind(14)=25	!Be9	 
 	 ind(15)=28	!B11
	  
c ordre des r�actions
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=39	!Be9(p,d)2He4
	 nb(16)=41	!Li6(p,He3)He4
	 nb(17)=42	!Li6(p,g)Be7
	 nb(18)=43	!Be9(p,a)Li6
 	 nb(19)=44	!B11(p,a)2He4
	 nb(20)=45	!B11(p,g)C12
  	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201

	CASE('ppcno11')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1005) ; WRITE(*,1005) 	
1005	  FORMAT('PP+CNO, 11 species, Be7 at equilibrium',/)	
	 CASE DEFAULT	
	  WRITE(2,5) ; WRITE(*,5) 	
5	  FORMAT('PP+CNO, 11 �l�ments, Be7 � l''�quilibre',/)	
	 END SELECT
	 
c nombre de r�actions et correspondance des indices 
	 nreac=14	!nombre de r�actions
	 nchim_reac=11	 
	 ihe4=4		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
	 	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7 
	 ind(6)=7	!C12
	 ind(7)=8	!C13	 
	 ind(8)=10	!N14	 
 	 ind(9)=11	!N15
 	 ind(10)=12	!O16
 	 ind(11)=13	!O17
 
c ordre des r�actions: on garde l'ordre de reac_c
	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201	 
		  
	CASE('ppcno12')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1006) ; WRITE(*,1006) 	
1006	  FORMAT('PP+CNO, 12 species',/)	
	 CASE DEFAULT	
	  WRITE(2,6) ; WRITE(*,6) 	
6	  FORMAT('PP+CNO, 12 �l�ments',/)
	 END SELECT
	 
c nombre de r�actions et correspondance des indices 
	 nreac=14	!nombre de r�actions
	 nchim_reac=12	 
	 ihe4=4		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
	 	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 
c ordre des r�actions: on garde l'ordre de reac_c
	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201	 
	
	CASE('ppcno12Be')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1007) ; WRITE(*,1007) 	
1007	  FORMAT('PP+CNO+Be9, 12 species+Be9',/)
	 CASE DEFAULT	
	  WRITE(2,7) ; WRITE(*,7) 	
7	  FORMAT('PP+CNO+Be9, 12 �l�ments+Be9',/)
	 END SELECT
	 
c nombre de r�actions et correspondance des indices 
	 nreac=16	!nombre de r�actions
	 nchim_reac=13	 
	 ihe4=4		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
 	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 	 ind(13)=25	!Be9
 
c ordre des r�actions: on garde l'ordre de reac_c
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=39
	 nb(16)=40
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201
	 
	CASE('ppcno12Li')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1008) ; WRITE(*,1008) 	
1008	  FORMAT('PP+CNO+Li6, 12 species',/)	
	 CASE DEFAULT	
	  WRITE(2,8) ; WRITE(*,8) 	
8	  FORMAT('PP+CNO+Li6, 12 �l�ments',/)
	 END SELECT
	 
c nombre de r�actions et correspondance des indices
	 nreac=15	!nombre de r�actions
	 nchim_reac=13	 
	 ihe4=4		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
	 	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 	 ind(13)=26	!Li6
 
c ordre des r�actions: on garde l'ordre de reac_c
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=41	!r�action 41 Li6(p,He3)He4
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201	 

	CASE('ppcno12BeBFe')
	 SELECT CASE (langue)
	 CASE('english')
	 WRITE(2,1013) ; WRITE(*,1013) 	
1013	 FORMAT('PP + CNO + Be + B + Fe, none species at equilibrium',/)	
	 CASE DEFAULT	
	  WRITE(2,13) ; WRITE(*,13) 	
13	  FORMAT('PP + CNO + Be + B + Fe, aucun �l�ment � l''�quilibre',/)
	 END SELECT
	 
c nombre de r�actions et correspondance des indices 
	 nreac=20	!nombre de r�actions
	 nchim_reac=17	 
	 ihe4=4		!indice de He4
	 i3al=0		!indice de la r�action 3 alpha
	 	 
c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1	!H1
	 ind(2)=2	!H2
	 ind(3)=3	!He3 
 	 ind(4)=4	!He4
	 ind(5)=5	!Li7
	 ind(6)=6	!Be7	  
	 ind(7)=7	!C12
	 ind(8)=8	!C13	 
	 ind(9)=10	!N14	 
 	 ind(10)=11	!N15
 	 ind(11)=12	!O16
 	 ind(12)=13  	!O17
 	 ind(13)=25	!Be9
 	 ind(14)=-1	!Ex	 
 	 ind(15)=28	!B11
 	 ind(16)=27	!Fe56	 
 	 ind(17)=26	!Li6
	 
c ordre des r�actions
	 nb(1:14)=(/ (i, i=1,14) /)	 
	 nb(15)=39	!Be9(p,d)2He4
	 nb(16)=41	!Li6(p,He3)He4
	 nb(17)=42	!Li6(p,g)Be7
	 nb(18)=43	!Be9(p,a)Li6
 	 nb(19)=44	!B11(p,a)2He4
	 nb(20)=45	!B11(p,g)C12
  	 
c temp�ratures max/min de tabulation

	 t_inf=0.5d6 ; t_sup=100.d6 ; n_temp=201
	CASE('ppcno3a9')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1009) ; WRITE(*,1009) 	
1009	  FORMAT('PP+CNO+3alpha, 9 species, H2, Li7, Be7 at quilibrium',/)
	 CASE DEFAULT	
	  WRITE(2,9) ; WRITE(*,9) 	
9	  FORMAT('PP+CNO+3al, 9 �l�ments, H2, Li7, Be7 � l''�quilibre',/)
	 END SELECT
	 
c nombre de r�actions et correspondance des indices
	 nreac=17	!nombre de r�actions
	 nchim_reac=9	 
	 ihe4=3		!indice de He4
	 i3al=15	!indice de la r�action 3 alpha

c indices des isotopes utilis�s dans les r�actions
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
 	 
c ordre des r�actions: on garde l'ordre de reac_c 
 	 nb(1:nreac)=(/ (i, i=1,nreac) /)	 
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=200.d6 ; n_temp=301	 
	  
	CASE('ppcno3ac10')
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1010) ; WRITE(*,1010) 	
1010	  FORMAT('PP+CNO+3al+C, 10 species, H2, Li7, Be7 at equilibrium',
	1 /,'CAUTION : C12(C12,a)O16 ignored')	
	 CASE DEFAULT	
	  WRITE(2,10) ; WRITE(*,10) 	
10	  FORMAT('PP+CNO+3al+C, 10 �l�ments, H2, Li7, Be7 en �quilibre',
	1 /,'ATTENTION : C12(C12,a)O16 pas prise en compte dans taux_nuc')
	 END SELECT
	 
c nombre de r�actions et correspondance des indices
	 nreac=23	!nombre de r�actions
	 nchim_reac=10	 
	 ihe4=3		!indice de He4
	 i3al=15	!indice de la r�action 3 alpha

c indices des isotopes utilis�s dans les r�actions 
	 ind(1)=1		!H1
	 ind(2)=3		!He3	 
	 ind(3)=4		!He4	 
 	 ind(4)=7		!C12
	 ind(5)=8		!C13	 
	 ind(6)=10		!N14
	 ind(7)=11		!N15	 
	 ind(8)=12		!O16	 
 	 ind(9)=13		!O17
 	 ind(10)=16		!Ne20
 	 
c ordre des r�actions: on garde l'ordre de reac_c
 	 nb(1:22)=(/ (i, i=1,22) /)	 
	 nb(23)=29
 	 
c temp�ratures max/min de tabulation
	 t_inf=0.5d6 ; t_sup=300.d6 ; n_temp=401	 
	  
	CASE DEFAULT
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1017)nom_nuc ; WRITE(*,1017)nom_nuc
1017	  FORMAT('STOP, unknown routine of thermonuclear reactions : ',a,/,
	1 'known routines :')
	 CASE DEFAULT
	  WRITE(2,17)nom_nuc ; WRITE(*,17)nom_nuc
17	  FORMAT('ARRET, routine de r�actions nucl�aires inconnue : ',a,/,
	1 'routines connues :')
	 END SELECT
	 WRITE(2,18) ; WRITE(*,18) ; STOP	
18	 FORMAT('pp1, pp3, ppcno9, ppcno10, ppcno10Fe, ppcno10K, ppcno11',/,
	1 'ppcno12, ppcno12Be, ppcno12Li, ppcno3a9, ppcno3ac10')
	END SELECT
	
c on limite t_stop � 0.95 t_sup pour �viter des d�passements de tabulation
	t_stop=MIN(0.95d0*t_sup,t_stop)


c s�lection du type de compilation � utiliser	
	SELECT CASE(nom_nuc_cpl)
	CASE('Cau-Fow')
	 total=0
	CASE('Adelb')
	 total=1	
	CASE('NACRE')
	 total=2	
	CASE('NACRE_LUNA')
	 total=2	
	CASE DEFAULT
	 SELECT CASE (langue)
	 CASE('english')
	  WRITE(2,1019)nom_nuc_cpl ; WRITE(*,1019)nom_nuc_cpl	
1019	  FORMAT('STOP, unknown compilation of thermonuclear reactions : ',
	1 a,/,'known compilations :')
	 CASE DEFAULT	
	  WRITE(2,19)nom_nuc_cpl ; WRITE(*,19)nom_nuc_cpl	
19	  FORMAT('ARRET, compilation de r�actions nucl�aires inconnue : ',
	1 a,/,'compilations connues :')
	END SELECT	 
	WRITE(2,20) ; WRITE(*,20) ; STOP
20	FORMAT('Cau-Fow, Adelb, NACRE')	
	END SELECT

c	initiations (ti=20.d6) des �nergies: q0, masses r�duites: ar,
c	masses des noyaux: nucleo, charges des noyaux: zi
c	charges des noyaux Z1, Z2 des noyaux de la r�action : izz
c	noms des r�action: nom_reac 	
c	appel pour remplir le common/reac_nuc/
c	identification des q, nom_reac, z1, z2

	ALLOCATE(q0(nreac),ar(nreac),nom_reac(nreac),izz(nreac,2),
	1 nom_elem(nchim),nucleo(nchim),zi(nchim))
	
	ti=20.d6 ; CALL taux_nuc(ti,total,rt,zit,izzt,qt,nom_react,
	1 nucleot,at,nom_elemt)	
	DO i=1,nreac
	 q0(i)=qt(nb(i))*eve*1.d6/amu	!�nergie des r�ac. en erg/r�action
	 ar(i)=at(nb(i)) ; nom_reac(i)=nom_react(nb(i))
	 izz(i,:)=izzt(nb(i),:)
	ENDDO
	
	DO i=1,nchim_reac	 
	 nom_elem(i)=nom_elemt(ind(i)) ; nucleo(i)=nucleot(ind(i))
	 zi(i)=zit(ind(i))
	ENDDO

c tabulation des points en ln(temp�rature)

	t0=LOG(t_inf) ; t1=LOG(t_sup)
	
	ALLOCATE(temp(n_temp))
	pas=(t1-t0)/DFLOAT(n_temp-1)
	DO i=1,n_temp
	 temp(i)=t0+pas*REAL(i-1)	 
c	 WRITE(*,2000)temp(i),EXP(t(i))
	ENDDO

c tabulation des r�actions

	ALLOCATE(taux_reac(nreac,n_temp))
	DO i=1,n_temp
	 CALL taux_nuc(EXP(temp(i)),total,rt,zit,izzt,qt,nom_react,
	1 nucleot,at,nom_elemt)
	 DO j=1,nreac
	  taux_reac(j,i)=rt(nb(j))
	 ENDDO	!j
c	 WRITE(*,2000)temp(i),EXP(temp(i)),rt(1) ; WRITE(*,2000)rt(2:9)
c	 WRITE(*,2000)rt(10:17)
	ENDDO	!i
	
	m_temp=4 ; ALLOCATE(ttemp(n_temp+m_temp))
	CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .FALSE.,temp(1),j,rt,qt)		!rt, qt: Vt
	
c �critures diverses	 
	SELECT CASE (langue)
	CASE('english')
	 WRITE(2,1011) ; WRITE(*,1011)
1011	 FORMAT(/,'isotopes employed :')	
	 WRITE(2,14)nom_elem(1:nchim_reac)
	 WRITE(*,14)nom_elem(1:nchim_reac)
14	 FORMAT(a4)  
	 WRITE(2,1012) ; WRITE(*,1012)
1012	 FORMAT(/,'thermonuclear reactions at work :')
	CASE DEFAULT
	 WRITE(2,11) ; WRITE(*,11)
11	 FORMAT(/,'isotopes utilis�s dans les r�actions :')	
	 WRITE(2,14)nom_elem(1:nchim_reac)
	 WRITE(*,14)nom_elem(1:nchim_reac) 
	 WRITE(2,12) ; WRITE(*,12)
12	 FORMAT(/,'r�actions nucl�aires utilis�es :')
	END SELECT	
	DO i=1,nreac
	 WRITE(2,15)nom_reac(i) ; WRITE(*,15)nom_reac(i)
15	 FORMAT(a20)	  	  
	ENDDO
	WRITE(2,*) ; WRITE(*,*)
		
	RETURN
	 
	END SUBROUTINE tabul_nuc
