
c******************************************************************

	SUBROUTINE resout_rota(dt,ok)

c routine public du module mod_evol

c résolution par collocation du système d'équa. diff.
c non linéaires de la diffusion du moment cinétique
c par itération Newton-Raphson, formalisme de Matis & Zahn 2004, Krot=4

c Auteur: P.Morel, Département Cassiopée, O.C.A.

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, Krot, nom_rot, nrot, ord_rot,
	1 precix, rot_min
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, bvald, gauss_band, intgauss,
	1 linf, no_croiss
	USE mod_variables, ONLY : dim_rot, knotr, mrot, mrott, n_rot, rota

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt
	LOGICAL, INTENT(out) :: ok

	REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a, b, temp
	REAL(kind=dp), DIMENSION(nrot,0:1) :: y
	REAL(kind=dp), DIMENSION(0:1,ord_rot) :: d
	REAL(kind=dp), PARAMETER :: explose=1.d2
	REAL(kind=dp), SAVE :: dtp=-HUGE(dt), err_lim, err_max
	REAL(kind=dp) :: er, err, errg, nu

	INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: iter_rot
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc
	INTEGER, PARAMETER :: compt_max=10 , itest2=7
	INTEGER, SAVE :: bloc, neq_c, neq_e 
	INTEGER :: compt, fait, ipe, iv, ive, k, ligne, nl, rang,
	1 spl, sple

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: inversible

	CHARACTER (len=2) :: czr

c----------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(10es8.1)
2003	FORMAT(i3,9es8.1)

c initialisations
	IF(init)THEN
	 init=.FALSE.

c bloc: longueur du bloc des coeff. non idt. nuls
	 bloc=nrot*ord_rot

c les niveaux de précision
	 err_max=10.d0*precix
	 err_lim=10.d0*err_max

c vecteurs de contrôle
	 ALLOCATE(iter_rot(nrot)) ; iter_rot=0
	 
c nombre de conditions limites au centre et surface
c controle des itérations sur Omega et Phi pour Krot=4
	 SELECT CASE (Krot)
	 CASE(3)
	  neq_c=4 ; neq_e=3 ; iter_rot(1)=1
	 CASE(4)
	  neq_c=5 ; neq_e=3 ; iter_rot(1)=1 ; iter_rot(7)=7  	 
	 END SELECT	 

	ENDIF

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c n_rot : nombre VARIABLE de points, élément de mod_variables
c nrot : nombre FIXE de fonctions, élément de mod_donnees
c m_rot : ordre FIXE des splines, élément de mod_donnees
c ord_rot : ordre d'interpolation, ord_rot=m_rot+1
c mrot(n_rot) : abscisses VARIABLES, élément de mod_variables

c dim_rot : dimension de l'espace des splines, élément de mod_variables
c défini dans base_rota, dim_rot = knotr-ord_rot

c nl = rang : nombre de lignes du système linéaire
	rang=nrot*dim_rot ; nl=rang+nrot
	
c allocations
	ALLOCATE(a(nl,bloc),b(1,nl),indpc(nl),temp(nrot,dim_rot))

c itérations Newton Raphson: boucle infinie B1
	compt=0
	B1: DO

c initialisations
	 a=0.d0 ; b=0.d0 ; ligne=0

c---------------début de la construction du système linéaire------------

c pour chaque point de collocation
	 DO k=1,ncoll_rot
	  nu=xcoll_rot(k)
  
c formation du jacobien, fait=0
	  fait=0
	  CALL jacob_rota(fait,nrot,nu)
	 ENDDO 		!k

c~~ au centre
	 fait=1 ; nu=mrot(1)
	 CALL jacob_rota(fait,neq_c,nu)	 

c limite externe toujours convectif fait=2
	 fait=2 ; nu=mrot(n_rot)
	 CALL jacob_rota(fait,neq_e,nu)

c----------------fin de la construction du système linéaire--------------

c k ligne qui pose Pb
c	 k=7086 ;  PRINT*,k,indpc(k) ;  WRITE(*,2000)a(k,:),b(1,k)

c ipe ligne qui pose Pb
c	 ipe=739	 
c	 PRINT*,ligne,nl,bloc ; PAUSE'ligne,nl,bloc avant a'
c	 DO k=1,nl
c	 DO k=ipe-nrot,ipe+nrot
c	 DO k=500,510
c	 DO k=5,nl,5
c	 DO k=nl-8,nl
c	  IF(indpc(k) == 721)THEN	   
c	   PRINT*,k,indpc(k)
c	   WRITE(*,2000)a(k,:),b(1,k)
c	  ENDIF
c	  IF(mod(i,200) == 0)PAUSE'a200'
c	  IF(MAXVAL(ABS(a(k,:))) == 0.d0)PAUSE'a=0'
c	 ENDDO
c	 PAUSE'a'
c	 Write(*,2000)b(1,5:nl:5) ; PAUSE'2-nd membres'
c	 PRINT*,nl
c	 DO i=1,nl
c	  IF(MAXVAL(ABS(a(i,:))) == 0.d0)PRINT*,i
c	 ENDDO
c	 PAUSE'lignes nulles'

c résolution du système linéaire
	 CALL gauss_band(a,b,indpc,ligne,rang,bloc,1,inversible)
	 IF(.NOT.inversible)THEN
	  PRINT*,'ARRET, matrice singulière dans resout_rota, nl=',ligne,
	1 ', rang=',rang
	  STOP
	 ENDIF

c	 WRITE(*,2000)b(1,:) ; PAUSE'solution'

c construction du vecteur temporaire des corrections
	 temp=RESHAPE(b,SHAPE(temp))

c	 B5: DO k=1,dim_rot
c	  PRINT*,k
c	  WRITE(*,2000)temp(1,k),rota(1,k),temp(3,k),rota(3,k),mrott(k)**(2./3.)
c	  WRITE(*,2000)temp(7,k),rota(7,k),temp(8,k),rota(8,k),mrott(k)**(2./3.)
c	  WRITE(*,2000)temp(:,k)
c	  WRITE(*,2000)rota(:,k)
c	  IF(MINVAL(ABS(rota(:,k))) == 0.d0)CYCLE B5
c	  WRITE(*,2000)(temp(j,k)/rota(j,k),j=1,nrot)
c	  WRITE(*,2000)rota(:,k)
c	  IF(mod(k,200) == 0)PAUSE'corrections'
c	 ENDDO B5
c	 PAUSE'corrections'

c estimation de la précision
	 err=0.d0 ; errg=0.d0
	 DO spl=1,dim_rot
	  B4: DO iv=1,nrot
	   IF(ABS(rota(iv,spl)) < rot_min(iv))CYCLE B4
	   er=ABS(temp(iv,spl)/rota(iv,spl))
	   IF(er > err)THEN
	    sple=spl ; ive=iv ; err=er
c	    PRINT*,iv,spl,nom_rot(iv)
c	    WRITE(*,2000)er,err,temp(iv,spl),rota(iv,spl),rot_min(iv)
	   ENDIF
	   IF(COUNT(ive == iter_rot) == 1)errg=MAX(errg,er)	   	   
	  ENDDO B4	!var
	 ENDDO		!spl
c	 PRINT*,sple,ive ; WRITE(*,2000)er,err ; PAUSE

c écritures
	 ipe=sple/ord_rot+1		!indice de couche
	 compt=compt+1
	 IF(lmix(mrott(sple)))THEN
	  czr='ZC'
	 ELSE
	  czr='ZR'
	 ENDIF
	 SELECT CASE (Krot)
	 CASE(3)	 
	  WRITE(usl_evol,100)compt,errg,nom_rot(1),czr,err,
	1 nom_rot(ive),sple,ipe,mrott(sple)**(2.d0/3.d0)
100	  FORMAT('resout_rota, itération:',i3,', erreur max:',es8.1,
	1 ', sur ',a,', dans ',a2/,'err. max. globale:',es8.1,
	2 ', sur ',a,', B-spline:',i4,', couche:',i4,', m=',es10.3)		
	  WRITE(usl_evol,101)mrot(ipe),nom_rot,rota(:,sple),temp(:,sple)
101	  FORMAT(/,'variables et corrections au maximum d''erreur, nu=',es10.3,/,
	1 7(3x,a4,3x),/,7es10.3,/,7es10.3,/)
	 CASE(4)	 
	  WRITE(usl_evol,103)compt,errg,nom_rot(1),nom_rot(itest2),czr,err,
	1  nom_rot(ive),sple,ipe,mrott(sple)**(2.d0/3.d0)
103	  FORMAT('resout_rota, itération:',i3,', erreur max:',es8.1,
	1 ', sur ',a,' ou ',a,', dans ',a2/,'err. max. globale:',es8.1,
	2 ', sur ',a,', B-spline:',i4,', couche:',i4,', m=',es10.3)	 	
	  WRITE(usl_evol,102)mrott(sple),nom_rot,rota(:,sple),temp(:,sple)
102	  FORMAT(/,'variables et corrections au maximum d''erreur, nu=',es10.3,/,
	1 8(3x,a4,3x),/,8es10.3,/,8es10.3,/)
	 END SELECT
c	 DO i=1,n_rot
c	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,mrot(i),
c	1 l,bs,bd)
c	  WRITE(*,2000)SQRT(mrot(i))**3,bs
c	 ENDDO
c	 PAUSE'S'

c	 PAUSE'avant écritures'
c	 DO i=1,dim_rot
c	  PRINT*,i
c	  WRITE(*,2000)temp(1,i),rota(1,i),temp(3,i),rota(3,i)
c	  WRITE(*,2000)temp(:,i)
c	  WRITE(*,2000)rota(:,i)
c	  IF(mod(i,200) == 0)PAUSE'corrections'
c	 ENDDO
c	 PAUSE'écritures'

c limitation des corrections
c	 DO spl=1,dim_rot
c	  B5: DO iv=1,nrot
c	   IF(ABS(rota(iv,spl)-temp(iv,spl)) < rot_min(iv))CYCLE B5
c	   temp(iv,spl)=SIGN(MIN(ABS(temp(iv,spl)),0.6d0*ABS(rota(iv,spl))),
c	1  temp(iv,spl))
c	  ENDDO B5	!iv
c	 ENDDO		!spl

c corrections
	 rota=rota-temp
	 
c écriture directe de U (NB. m_rot=ord_rot-1)
c	 WRITE(*,2000)(rota(2,(ord_rot-1)*(iv-1)+1),iv=1,n_rot) ; PAUSE'U direct' 

c mises à 0 en dessous de 1.d-20
c	 WHERE(ABS(rota) < 1.d-20)rota=0.d0

c au centre toutes les variables de la rotation sont identiquement
c nulles sauf Omega
c	 rota(2:nrot,1)=0.d0

c gestion de la précision
c on poursuit si l'erreur est < err_lim après compt_max
c il y a convergence forcée si l'erreur est < err_max après 2*compt_max
c pas de test de précision pour le premier pas temporel
c sortie avec une erreur > explose
c	 PRINT*,compt,compt_max
c	 WRITE(*,2000)errg,precix,err_max,err_lim

	 IF(compt <= 1)CYCLE B1
	 ok= errg <= precix
	 IF(ok .OR. errg > explose)THEN
	  EXIT B1
	 ELSEIF(compt < compt_max)THEN
	  CYCLE B1
	 ELSEIF(errg < err_max)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   PRINT*,'bad conv. in resout_rota, number of iterations :',compt
	  CASE DEFAULT
	   PRINT*,'conv. forcée dans resout_rota, nb. itérations :',compt
	  END SELECT
	  ok=.TRUE. ; EXIT B1
	 ELSEIF(errg < err_lim .AND. compt < 2*compt_max)THEN
	  CYCLE B1
	 ELSEIF(dt /= dtp)THEN
	  ok=.TRUE. ; EXIT B1
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   PRINT*,'no conv. in resout_rota, number of iterations :',compt
	  CASE DEFAULT
	   PRINT*,'pas de conv. dans resout_rota, nb. itérations :',compt
	  END SELECT
	  EXIT B1
	 ENDIF
	ENDDO B1
	dtp=dt

c nettoyage
	DEALLOCATE(a,b,indpc,temp)

c	PAUSE'sortie resout_rota'

	RETURN

	CONTAINS
	 INCLUDE'jacob_rota.f'

	END SUBROUTINE resout_rota
