
c*******************************************************

	SUBROUTINE base_chim

c routine subordonnée de diffus	
c Formation de la base pour la diffusion des éléments chimiques

c la fonction d'interpolation est d'ordre m
c le vecteur nodal a m point en x(1), x(n) et m-i aux singularités
c s'il n'y a pas de singularité une base continue est crée
c s'il n'y a des ZR & ZC base continue non dérivable aux limites
c avec diffusion du moment cinétique, base continue non dérivable en tout point

c entrées
c	0 <= i < m : ordre de continuité aux points de singularité
c	i = 0 discontinue
c	i = 1 dérivée première discontinue
c	i = 2 dérivée seconde discontinue
c	.................
c	i = m-1 fonction continue partout

c	is(0:ns+1) : indices des abscisses des singularités
c	ns : nombre de singularités

c Auteur: P. Morel, Département Cassiopée, O.C.A.
c CESAM2k

c------------------------------------------------------------------

	USE mod_donnees, ONLY : Krot
	USE mod_variables, ONLY : tot_conv
		
	IMPLICIT NONE

	INTEGER, ALLOCATABLE, DIMENSION(:) :: mult
	
	LOGICAL, SAVE :: init=.TRUE.

c------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 SELECT CASE(langue)
	 CASE('english')
          WRITE(*,1001) ; WRITE(2,1001)
1001      FORMAT(/,'Use of a continuous basis')	 
	 CASE DEFAULT 
          WRITE(*,1) ; WRITE(2,1)
1         FORMAT(/,'Utilisation de la base continue')
	 END SELECT	
	ENDIF
		
c multiplicités
	ALLOCATE(mult(n_ch))

c continuité
	IF(tot_conv)THEN
	 mult=1
	 	 
c discontinuité de la dérivée 1-ière 
	ELSE
	 SELECT CASE(Krot)
	 CASE(3,4)		!partout	
	  mult=MAX(1,m_ch-1)
	 CASE DEFAULT		!seulement aux limites ZR/ZC
	  mult=1 ; mult(idis(1:ndis))=MAX(1,m_ch-1)
	 END SELECT	 
	ENDIF	
		
c construction du vecteur nodal				
	mult(1)=m_ch ; mult(n_ch)=m_ch ; knotc=SUM(mult)
	DEALLOCATE(mct) ; ALLOCATE(mct(knotc))
	CALL noeud(mc,mct,mult,n_ch,knotc)
	DEALLOCATE(mult)
	
	RETURN
	
	END SUBROUTINE base_chim
