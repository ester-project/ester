
c***************************************************************************

	SUBROUTINE initialise_rota

c routine PRIVATE du module mod_evol
c initialisation provisoire de la rotation rota(Omega,U,theta,Lambda,Psi)
c cette initialisation, effectuée sur le nombre initial de couches avec
c de bsp1dn sera reprise dans diffus pour le nombre définitif de couches

c pour un modèle de ZAMS ou de PMS initial ou encore pour un modèle
c repris en ASCII les masses de rota et de bp coïncident

c Auteur: P.Morel, Département Cassiopée, O.C.A.
c CESAM2k

c------------------------------------------------------------------------------

	USE mod_donnees, ONLY : en_masse, g, Krot, msol, m_ch, m_rot, ne,
	1 nrot, ord_qs, ord_rot, pi, rsol, w_rot
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : bp, dim_ch, dim_rot, knot, knotc,
	1 knotr, mc, mct, mrot, mrott, n_ch, n_qs, n_rot, q, qt, rota, sortie

	IMPLICIT NONE

	REAL (kind=dp), DIMENSION(ne) :: dfqs, fqs
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot
	REAL (kind=dp), SAVE :: cte1
	REAL (kind=dp) :: r

	INTEGER :: i, l=1

c-----------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	cte1=16.d0*pi*rsol**6/9.d0/g/msol**2

c initialisation de l'interpolation de la rotation
c pour plus de clarté, certaines définitions sont volontairement naïves

	SELECT CASE(Krot)
	CASE(3)		!Formalisme de Talon & Zahn 97
	 n_rot=n_ch ; knotr=n_rot+ord_rot ; dim_rot=knotr-ord_rot	
	 ALLOCATE(mrot(n_rot),mrott(knotr),rota(nrot,dim_rot))
	 mrot=mc ; rota=0.d0	
	
c initialisation de rota(1,i)=Omega
	rota(1,:)=ABS(w_rot)

c initialisation de rota
	 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.FALSE.,mrot(1),
	1 l,frot,dfrot)
	 IF(no_croiss)PRINT*,'Pb. en 3 dans initialise_rota'
	
c initialisation de Theta et Psi=-Theta (n_ch et n_qs coïncident)	
	 DO i=1,n_rot
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),l,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. en  4 dans initialise_rota'
	  IF(en_masse)THEN
	   r=SQRT(ABS(fqs(3)))
	  ELSE
	   r=ABS(fqs(3))
	  ENDIF
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,mrot(i),
	1 l,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb. en 5 dans initialise_rota'
	  IF(mrot(i) <= 0.d0)THEN
	   rota(3,i)=0.d0
	  ELSE
	   rota(3,i)=cte1*r**6*(dfrot(1)-dfrot(2))/mrot(i)**2*frot(1)*dfrot(1)
	  ENDIF
	  rota(5,i)=-rota(3,i)
	 ENDDO

c pour rota4
	CASE(4)		!Formalisme de Mathis & Zahn 04
	 CALL tab_vth ; CALL initialise_rota4

	CASE(5)		!Conservation locale du moment cinétique
	 n_rot=n_ch ; ord_rot=m_rot+1 ; knotr=n_rot+ord_rot
	 dim_rot=knotr-ord_rot	
	 ALLOCATE(mrot(n_rot),mrott(knotr),rota(nrot,dim_rot))
	 mrot=mc ; rota=0.d0 ; rota=w_rot
	 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.FALSE.,mrot(1),
	1 l,frot,dfrot)		 
	END SELECT

c	DO i=1,n_rot
c	 WRITE(*,2000)mrot(i),rota(:,i)
c	ENDDO
c	WRITE(*,2000)rota(1,:) ; Pause'ini'

	RETURN

	END SUBROUTINE initialise_rota

