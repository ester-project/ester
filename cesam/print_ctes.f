
c***************************************************************

	SUBROUTINE print_ctes(i)

c	subroutine public du module mod_donnees	
c	écriture sur l'unité i des constantes utilisées
c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c---------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	INTEGER, INTENT(in) :: i

c---------------------------------------------------------------------

	WRITE(i,1)source
1	FORMAT(//,t5,'Constantes fondamentales selon ',a50,/)
	WRITE(i,2)msol,lsol,rsol,echarg,kbol,hpl,g,aradia,clight,
	1 eve,amu,me,lbol0,gmsol
2	FORMAT('masse solaire=',es15.8,/,
	1 'luminosité solaire=',es15.8,/,
	2 'rayon solaire=',es15.8,/,
	3 'charge de l''électron=',es15.8,/,
	4 'constante de Boltzmann=',es15.8,/,
	5 'constante de Planck=',es15.8,/,
	6 'constante gravitation=',es15.8,/,
	7 'constante radiation=',es15.8,/,
	8 'célérité de la lumière=',es15.8,/,
	9 'électron volt=',es15.8,/,
	1 'masse atomique unité=',es15.8,/,
	2 'masse de l''électron=',es15.8,/,
	3 'zéro des Mbol L=',es15.8,/,
	4 'constante de Newton G Msol=',es15.8)
 	WRITE(i,6)
6	FORMAT(/,' masses des noyaux, données de NACRE 1999',/)
	WRITE(i,4)ah, ah2, ahe3, ahe4, ali7, ali6, abe7, abe9, ac12, ac13,
	1 an13, an14, an15, ao16, ao17, ao18, af18, af19, ane20, ane22,
	2 ana23, amg24, amg25, amg26, afe56, an, ap
4	FORMAT('    H=',es15.8,
	1  t25,'  H2=',es15.8,
	2  t50,' He3=',es15.8,/,
	3      '  He4=',es15.8,
	4  t25,' Li6=',es15.8,
	5  t50,' Li7=',es15.8,/,
	6      '  Be7=',es15.8,
	7  t25,' Be9=',es15.8,
	8  t50,' C12=',es15.8,/,
	9      '  C13=',es15.8,
	1  t25,' N13=',es15.8,
	2  t50,' N14=',es15.8,/,
	3      '  N15=',es15.8,
	4  t25,' O16=',es15.8,
	5  t50,' O17=',es15.8,/,
	6      '  O18=',es15.8,
	7  t25,' F18=',es15.8,	
	8  t50,' F19=',es15.8,/,
	9      ' Ne20=',es15.8,
	1  t25,'Ne22=',es15.8,
	2  t50,'Na23=',es15.8,/,
	3      ' Mg24=',es15.8,
	4  t25,'Mg25=',es15.8,
	5  t50,'Mg26=',es15.8,/,
	6      ' Fe56=',es15.8,
	7  t25,'neut=',es15.8,
	8  t50,'prot=',es15.8,/)
	
	RETURN

	END SUBROUTINE print_ctes
