
c******************************************************************

	SUBROUTINE initialise_rota4
	
c routine PRIVATE du module mod_evol

c initialisation de rota4, seul Omega et Poisson sont pris en compte

c sortie :
c	rota : du module mod_evol
	
c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : m_rot, nrot, ord_rot, r_qs
	USE mod_kind
	USE mod_numerique, ONLY : bvald, coll, gauss_band, linf, noeud
	USE mod_variables, ONLY : dim_rot, knotr, mc, mrot, mrott, n_ch, n_rot, rota

	IMPLICIT NONE

	REAL(kind=dp), DIMENSION(nrot,nrot,0:1) :: ae
	REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a, b
	REAL(kind=dp), DIMENSION(0:1,ord_rot) :: d	
	REAL(kind=dp), DIMENSION(nrot) :: be
	
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc, mult
	
	INTEGER, SAVE ::  bloc
	INTEGER ::  col, i, id, ie, ind, indice, iv, j, k, l=1, ligne, nl, rang

	LOGICAL :: inversible
	
c----------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(10es10.3)

c allocations, et vecteur nodal, on ne tient pas compte de discontinuités
	n_rot=n_ch ; ALLOCATE(mrot(n_rot),mult(n_rot)) ; mrot=mc
	mult=m_rot ; mult(1)=m_rot+r_qs ; mult(n_rot)=m_rot+r_qs
	knotr=SUM(mult) ; ord_rot=m_rot+r_qs ; dim_rot=knotr-ord_rot
	ALLOCATE(mrott(knotr))
	CALL noeud(mrot,mrott,mult,n_rot,knotr) ; DEALLOCATE(mult)
	
c les points de collocation
	ncoll_rot=(n_rot-1)*m_rot ; ALLOCATE(xcoll_rot(ncoll_rot))
	CALL coll(mrot,n_rot,m_rot,xcoll_rot)
		 
c longueur du bloc des coeff. non idt. nuls	
	bloc=nrot*ord_rot
	 
c nl = rang : nombre de lignes du système linéaire
	rang=nrot*(ncoll_rot+r_qs) ; nl=rang	
	
c allocations et initialisations 
	ALLOCATE(a(nl,bloc),b(1,nl),indpc(nl))
	a=0.d0 ; b=0.d0 ; ligne=0	
	
c pour chaque point de collocation qui sont ceux de rota	
	DO k=1,ncoll_rot

c formation des coefficients des équations
	 CALL eq_ini_rota4(0,xcoll_rot(k),ae,be)

c les splines	 	
	 CALL linf(xcoll_rot(k),mrott,knotr,l)	 
	 CALL bvald(xcoll_rot(k),mrott,ord_rot,l,1,d)

c contributions au système linéaire
	 ind=nrot*(l-ord_rot)+1	!ind. prem. col. for ne next lignes
	 DO i=1,nrot
	  ligne=ligne+1 ; b(1,ligne)=be(i) ; indpc(ligne)=ind
	  DO j=1,ord_rot		!j-ième spline
	   DO iv=1,nrot
	    indice=nrot*(j-1)+iv
	    DO id=0,1 			!id-ème dérivée
	     a(ligne,indice)=a(ligne,indice)+ae(i,iv,id)*d(id,j)
	    ENDDO	!id
	   ENDDO	!iv
	  ENDDO		!j
	 ENDDO		!i
	ENDDO		!k

c les conditions limites au centre	
	CALL eq_ini_rota4(1,mrot(1),ae,be)	
	CALL linf(mrot(1),mrott,knotr,l)	!les splines
	CALL bvald(mrot(1),mrott,ord_rot,l,1,d)	!dérivée

c contribution au système linéaire,
c la matrice compressée est le jacobien 'diagonal' ie. sans les
c éléments 'non diagonaux' identiquement nuls
	ind=nrot*(l-ord_rot)+1	!ind. prem. col. for nrot next lignes
	DO ie=1,nrot/2		!pour chaque équation
	 ligne=ligne+1 ; b(1,ligne)=be(ie) ; indpc(ligne)=ind
	 DO j=1,ord_rot		!pour chaque spline j
	  DO iv=1,nrot	!pour chaque variable
	   col=nrot*(j-1)+iv
	   DO id=0,1
	    a(ligne,col)=a(ligne,col)+ae(ie,iv,id)*d(id,j)
	   ENDDO	!id
	  ENDDO		!iv variable
	 ENDDO		!j
	ENDDO		!ie équation	

c les conditions limites à la surface		
	 CALL eq_ini_rota4(2,mrot(n_rot),ae,be)	
	 CALL linf(mrot(n_rot),mrott,knotr,l)		!les splines
	 CALL bvald(mrot(n_rot),mrott,ord_rot,l,1,d)	!dérivée
	 
c contribution au système linéaire,
c la matrice compressée est le jacobien 'diagonal' ie. sans les
c éléments 'non diagonaux' identiquement nuls
	ind=nrot*(l-ord_rot)+1	!ind. prem. col. for nrot next lignes
	DO ie=1,nrot/2		!pour chaque équation
	 ligne=ligne+1 ; b(1,ligne)=be(ie) ; indpc(ligne)=ind
	 DO j=1,ord_rot		!pour chaque spline j
	  DO iv=1,nrot	!pour chaque variable
	   col=nrot*(j-1)+iv
	   DO id=0,1
	    a(ligne,col)=a(ligne,col)+ae(ie,iv,id)*d(id,j)
	   ENDDO	!id
	  ENDDO		!iv variable
	 ENDDO		!j
	ENDDO		!ie équation	
	
c----------------fin de la construction du système linéaire--------------

c résolution
	CALL gauss_band(a,b,indpc,nl,nl,bloc,1,inversible)
	IF(.NOT.inversible)THEN
	 PRINT*,'ARRRET, matrice singulière dans initialise_rota4' ; STOP
	ENDIF

c on place la solution dans rota
	ALLOCATE(rota(nrot,dim_rot))
	rota=RESHAPE(b(1,:),(/ nrot,dim_rot /))	
	
	DEALLOCATE(a,b,indpc)	
c	PAUSE'initialise_rota4 fin'	
	
	RETURN
	
	END SUBROUTINE initialise_rota4
