
c****************************************************************

	MODULE mod_cesam

c module contenant la routine cesam et les routines de gestion

c fonctions private:
c	add_ascii : formation des variables pour les sorties ASCII
c	ascii : création d'un fichier de sortie ASCII personalisé	
c	des : routine générique de dessin on line
c	des_m : dessin en fonction de la masse
c	des_r : dessin en fonction du rayon
c	dnunl : calcul approché de la petite et de la grande différences 	
c	list : formation du listing
c	output : routine générique des sorties en ASCII
c	osc_adia : sortie ASCII pour oscillations adiabatiques
c	osc_invers : sortie ASCII pour inversions
c	osc_nadia : sortie ASCII pour oscillations non-adiabatiques

c fonction public:
c	cesam : gestion générale du calcul
	
c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c    	CESAM2k

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	PRIVATE
	PUBLIC :: cesam

	CONTAINS

c*********************************************************************

	INCLUDE 'add_ascii.f'
	INCLUDE 'ascii.f'	
	INCLUDE 'cesam.f'
	INCLUDE 'des.f'
	INCLUDE 'des_m.f'
	INCLUDE 'des_r.f'
	INCLUDE 'dnunl.f'
	INCLUDE 'list.f'
	INCLUDE 'output.f'
	INCLUDE 'osc_adia.f'
	INCLUDE 'osc_invers.f'
	INCLUDE 'osc_nadia.f'

	END MODULE mod_cesam
	
	INCLUDE 'pgplot_factice.f'	
