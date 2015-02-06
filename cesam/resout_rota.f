
c******************************************************************

	SUBROUTINE resout_rota(dt,ok)

c	routine public du module mod_evol
	
c	résolution par éléments finis Galerkin du système d'équa. diff.
c	ord. non linéaires de la rotation par itération Newton-Raphson

c	fait=0 : on forme les conditions limites
c	fait=1:
c	on forme les produits scalaire pour toutes les splines de la base
c	continue - sauf la première te la dernière - par intégration de Gauss
c	(algorithme inspiré de celui de Schumaker, 5.22 p. 203),
c	les quantités intégrées sont nulles

c	 as : dérivées/ Xi, Xi' des coeff. de la spline du prod. scalaire
c	 ad : dérivées/ Xi, Xi' des coeff. de la der. du produit scalaire

c	 as(i,j,k)=coefficient de la k-ième dérivée (k=0,1)
c	 de la j-ième variable de la i-ième équation < . Ni >

c	 ad(i,j,k)=coefficient de la k-ième dérivée (k=0,1)
c	 de la j-ième variable de la i-ième équation < . dNi/dx > 

c	 a(<spline 1 . équation 1: variable 1,2,..,ne>
c	   <spline 1 . équation 2: variable 1,2,..,ne>
c		 :     :     :     :     :     
c	   <spline 1 . équation ne: variable 1,2,..,ne>
c	   <spline 2 . équation ne: variable 1,2,..,ne>
c	   <spline 2 . équation ne: variable 1,2,..,ne>
c           :     :     :     :     :
c	   <spline 2 . équation ne: variable 1,2,..,ne>
c	   <spline 3 . équation ne: variable 1,2,..,ne>
c           :     :     :     :     :
c	   <spline dim . équation ne: variable 1,2,..,ne>)

c	 bs(i)=second membre de la i-ième équation < . Ni > 
c	 bd(i)=second membre de la i-ième équation < . d/dx Ni >

c	 b(<spline 1 . équation 1, 2,.., ne>
c	  <spline 2 . équation 1, 2,.., ne>
c           :     :     :     :
c	  <spline dims . équation 1, 2,.., nvar>)

c Auteur: P.Morel, Département Cassiopée, O.C.A.

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, m_rot, nom_rot, nrot, precix,
	1 rot_min
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, bvald, gauss_band, intgauss,
	1 linf, no_croiss
	USE mod_variables, ONLY : dim_rot, knotr, mrot, mrott, n_rot, rota
	
	IMPLICIT NONE
	REAL (kind=dp), INTENT(in) :: dt
	LOGICAL, INTENT(out) :: ok
		
	REAL(kind=dp), DIMENSION(nrot,nrot,0:1) :: ad, as		
	REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a, b, temp
	REAL(kind=dp), DIMENSION(nrot,0:1) :: y
	REAL(kind=dp), DIMENSION(0:1,m_rot) :: d	
	REAL(kind=dp), DIMENSION(nrot) :: bd, bs
	REAL(kind=dp), DIMENSION(m_rot) :: qg, wg
	REAL(kind=dp), SAVE :: dtp=-HUGE(dt)
	REAL(kind=dp) :: corr, er, err, errg	
		
	INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc, indpc0
	INTEGER, PARAMETER :: compt_max=15, err_lim=10.d0, err_max=0.1d0,
	1 itest=3
	INTEGER, SAVE :: bloc
	INTEGER :: col, colonne, compt, i, id,
	1 idep, ie, ifin, ig, ipe, iv, ive, j, k, l=1, li, ligne,
	2 nl, rang, spi, spl, sple

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: inversible
		
c----------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(10es8.1)

c	 initialisations
c	bloc: longueur du bloc des coeff. non idt. nuls
	IF(init)THEN
	 bloc=nrot*(2*m_rot-1)
	 init=.FALSE.
	ENDIF		 

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees 
c	mc(n_ch) : abscisses VARIABLES élément de mod_variables

c	nl = rang : nombre de lignes du système linéaire
	
	rang=nrot*dim_rot ; nl=rang
		
	ALLOCATE(a(nl,bloc),b(1,nl),indpc(nl),indpc0(nl))

c	indices de première colonne pour les produits scalaires
c	de la base continue
c	ils seront changés dans gauss_band

	indpc0=0 ; ligne=0
	DO li=1,nlim_rot
	 CALL linf(xlim_rot(li),mrott,knotr,l)
	 DO ie=1,nrot
	  ligne=ligne+1 ; indpc0(ligne)=nrot*(l-m_rot)+1
	 ENDDO
	ENDDO
	DO i=2,dim_rot-1      !pour chaque spline de la base des mrot
	 k=MAX(1,i-m_rot+1)      !indice s'il n'y avait qu'une variable
	 DO ie=1,nrot         !avec toutes les variables
	  ligne=ligne+1 ; indpc0(ligne)=nrot*(k-1)+1
	 ENDDO
	ENDDO
	 
c	itérations Newton Raphson: boucle infinie B1

	compt=0
	B1: DO

c	 initialisations

	 a=0.d0 ; b=0.d0 ; ligne=0 ; indpc=indpc0
	 
c---------------début de la construction du sytème linéaire------------	

c	   pour chaque limite

	 DO li=1,nlim_rot
	
c	    localisation
	 
	  CALL linf(xlim_rot(li),mrott,knotr,l) ; spi=l-m_rot

c	    variables et dérivées en xlim_rot(li)
	
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1   xlim_rot(li),l,bs,bd)
	  IF(no_croiss)PRINT*,'resout_rota, Pb. en 1'	
	  y(:,0)=bs ; y(:,1)=bd
	 
c	    formation des coefficients pour la limite
	 
	  CALL eq_diff_rota(0,xlim_rot(li),y,ad,as,bd,bs)
	   		   
c	    les B-splines en , dérivées 0 à 1

	  CALL bvald(xlim_rot(li),mrott,m_rot,l,1,d)
	  	 
c	    contribution au système linéaire

	  DO ie=1,nrot	!pour chaque équation	 
	   ligne=ligne+1
	   DO i=1,m_rot		!pour chaque spline d'indice i
	    b(1,ligne)=b(1,ligne)+bs(ie)*d(0,i)+bd(ie)*d(1,i)

c	    la matrice compressée est le jacobien 'diagonal' ie. sans les
c	    éléments 'non diagonaux' identiquement nuls

	    DO j=1,m_rot		!pour chaque spline
	     DO iv=1,nrot		!pour chaque variable
	      colonne=nrot*(spi+j-1)+iv    
	      col=colonne-indpc(ligne)+1   !colonne matrice compressée
	      DO id=0,1  !pour 0: NiNj, 1: NiN'j, 3: NiN"j etc...
	       a(ligne,col)=a(ligne,col)+as(ie,iv,id)*d(id,j)*d(0,i)
	1       +ad(ie,iv,id)*d(id,j)*d(1,i)
	      ENDDO	!id
	     ENDDO	!iv variable   
	    ENDDO 	!i	  
c	      WRITE(*,2000)xl_chim(li),a(ligne,:),b(1,ligne) ; PAUSE'xl_chim'
	   ENDDO  
	  ENDDO
	 ENDDO		!li	  
	    
c	   les produits scalaires, pour chaque intervalle
c	   pour rota de m_rot+1 à dim_rot-1 (plus=1)
c	   pour chim de m_ch à dim_ch (plus=0)

	 B3: DO k=m_rot,dim_rot

c on évite les multiplicités	 
	  IF(mrott(k) >= mrott(k+1))CYCLE B3
	 
	  IF(k == m_rot)THEN
	   idep=2 ; ifin=m_rot
	  ELSEIF(k == dim_rot)THEN
	   idep=1 ; ifin=m_rot-1
	  ELSE
	   idep=1 ; ifin=m_rot	  
	  ENDIF    
	  
c	    spi : indice-1 de la première spline non id.nulle  
c	    comme mct(l) <= qg(ig) < mct(l+1), l est le même pour toutes
c	    les abscisses de Gauss dans l'intervalle [mct(k),mct(k+1)]

	  CALL linf(mrott(k),mrott,knotr,l) ; spi=l-m_rot	  

c	     poids, wg et abscisses, qg pour intégration Gauss d'ordre m_rot

	  CALL intgauss(mrott(k),mrott(k+1),qg,wg,m_rot)

c	     pour chaque abscisse de GAUSS

	  DO ig=1,m_rot

c	      variables et dérivées

	   CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1    qg(ig),l,bs,bd)

	   IF(no_croiss)PRINT*,'resout_rota, Pb. en 1'
	   y(:,0)=bs ; y(:,1)=bd
	   
c	      formation des coefficients des équations

	   CALL eq_diff_rota(1,qg(ig),y,ad,as,bd,bs)	    
	   		   
c	      les B-splines en qg(ig), dérivées 0 à 1

	   CALL bvald(qg(ig),mrott,m_rot,l,1,d)
	   
c	      contribution au système linéaire
    
	   DO ie=1,nrot		!pour chaque équation
	    DO i=idep,ifin		!pour chaque spline d'indice i
c	     ligne=nrot*(spi+i-1)+ie+nrot
	     ligne=nrot*(spi+i)+ie     
	     b(1,ligne)=b(1,ligne)+wg(ig)*(d(0,i)*bs(ie)+d(1,i)*bd(ie))
c	        WRITE(*,2003)ligne,d(0,i),d(1,i),bs(ie),bd(ie),y(1,0),
c	1       y(1,1),qg(ig),b(1,ligne) ; PAUSE'b'
2003	     FORMAT(i3,9es8.1)
		
c	        la matrice compressée est le jacobien 'diagonal' ie. sans les
c	        éléments 'non diagonaux' identiquement nuls

	     DO j=1,m_rot		!pour chaque spline j
	      DO iv=1,nrot	!pour chaque variable	       
	       colonne=nrot*(spi+j-1)+iv    
	       col=colonne-indpc(ligne)+1   !colonne matrice compressée
c	          PRINT*,'ligne,colonne,col,iv',ligne,colonne,col,iv;PAUSE'1'
	       DO id=0,1  !pour 0: NiNj, 1: NiN'j, 3: NiN"j etc...
	        a(ligne,col)=a(ligne,col)+wg(ig)*
	1	(as(ie,iv,id)*d(0,i)*d(id,j)+ad(ie,iv,id)*d(1,i)*d(id,j))
	       ENDDO    !id
	      ENDDO	!iv variable
	     ENDDO	!j
	    ENDDO	!i
	   ENDDO	!ie équation	   
	  ENDDO 	!ig
	 ENDDO B3	!k
	 
c----------------fin de la construction du système linéaire-------------- 

c	  PRINT*,ligne,nl,bloc ; PAUSE'ligne,nl,bloc avant a'
c	 DO i=1,nl
c	   PRINT*,i,indpc(i)
c	  WRITE(*,2000)a(i,:),b(1,i) ; IF(mod(i,200) == 0)PAUSE'a200'
c	  IF(MAXVAL(ABS(a(i,:))) == 0.d0)PAUSE'a=0'
c	 ENDDO
c	 PAUSE'a'
	
c	 PRINT*,nl
c	 DO i=1,nl
c	  IF(MAXVAL(ABS(a(i,:))) == 0.d0)PRINT*,i
c	 ENDDO
	 	
c	 résolution du système linéaire

	 CALL gauss_band(a,b,indpc,nl,rang,bloc,1,inversible)
	 IF(.NOT.inversible)THEN
	  PRINT*,'ARRRET, matrice singulière dans resout_rota' ; STOP
	 ENDIF
	 
c	 WRITE(*,2000)b(1,:) ; PAUSE'solution'
	
c	 construction du vecteur temporaire des corrections

	 ALLOCATE(temp(nrot,dim_rot)) ; temp=RESHAPE(b,SHAPE(temp))
	 
c	 B5: DO i=1,dim_rot
c	  PRINT*,i
c	  WRITE(*,2000)temp(1,i),rota(1,i),temp(3,i),rota(3,i)
c	  WRITE(*,2000)temp(:,i)
c	  WRITE(*,2000)rota(:,i)
c	  IF(MINVAL(ABS(rota(:,i))) == 0.d0)CYCLE B5
c	  WRITE(*,2000)(temp(j,i)/rota(j,i),j=1,nrot)	  
c	  WRITE(*,2000)rota(:,i)
c	  IF(mod(i,200) == 0)PAUSE'corrections'
c	 ENDDO B5
c	 PAUSE'corrections'
	 	
c	 limitation des corrections

	 corr=1.d0
	 DO spl=1,dim_rot
	  B2: DO iv=1,nrot	  
	   IF(iv == 1 .OR. iv == itest)THEN	  
	    IF(ABS(rota(iv,spl)) < rot_min(iv))CYCLE B2
	    IF(corr*ABS(temp(iv,spl)) > 0.6d0*ABS(rota(iv,spl)))
	1     corr=corr/2.d0
	   ENDIF
	  ENDDO B2	!iv
	 ENDDO		!spl

c	 estimation de la précision

	 err=0.d0 ; errg=0.d0
	 DO spl=1,dim_rot
	  B4: DO iv=1,nrot
	   IF(ABS(rota(iv,spl)) < rot_min(iv))CYCLE B4
	   er=ABS(temp(iv,spl)/rota(iv,spl))
	   IF(er > err)THEN
	    sple=spl ; ive=iv ; err=er
c	       PRINT*,iv,spl,nom_rot(iv)
c	       WRITE(*,2000)er,err,temp(iv,spl),rota(iv,spl),rot_min(iv)
	   ENDIF
	   IF((iv == 1 .OR. iv == itest) .AND. er > errg)errg=er 
	  ENDDO B4	!var
	 ENDDO		!spl
	 
c	 PRINT*,sple,ive ; WRITE(*,2000)er,err ; PAUSE	
	 
c	 écritures

	 ipe=sple/m_rot+1
	 
	 compt=compt+1
	 WRITE(*,100)compt,errg,nom_rot(1),nom_rot(itest),corr,err,
	1  nom_rot(ive),sple,ipe
100	 FORMAT('resout_rota, itération:',i3,', erreur max:',es8.1,
	1 ', sur ',a,' ou ',a,', corr:',es8.1,/,'err. max. globale:',es8.1,
	2 ', sur ',a,', B-spline:',i4,', couche:',i4)
	 
c	 DO i=1,n_rot
c	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,mrot(i),
c	1  l,bs,bd)
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

c corrections et suppression de temp
	 rota=rota-temp*corr ; DEALLOCATE(temp)

c on poursuit si l'erreur est < err_lim après compt_max 
c il y a convergence forcée si l'erreur est < err_max après2*compt_max
c on accepte une mauvaise convergence pour un nouveau pas temporel
	 ok= errg <= precix	 
	 IF(errg > 1.d2)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   PRINT*,'too large error in resout_rota, number of iterations :',
	1   compt	  
	  CASE DEFAULT
	   PRINT*,'correction >> 1 dans resout_rota, nb. itérations :',
	1   compt
	  END SELECT
	  PRINT*,'dt --> dt/2' ; EXIT B1
	 ELSEIF(compt <= 1)THEN
	  CYCLE B1
	 ELSEIF(ok)THEN
	  EXIT B1
	 ELSEIF(compt < compt_max)THEN
	  CYCLE B1	 
	 ELSE
	  IF(errg < err_max)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    PRINT*,'bad conv. in resout_rota, number of iterations :',
	1    compt	  
	   CASE DEFAULT
	    PRINT*,'conv. forcée dans resout_rota, nb. itérations :',
	1    compt
	   END SELECT
	   ok=.TRUE. ; EXIT B1
	  ELSEIF(errg < err_lim .AND. compt < 2*compt_max)THEN
	   CYCLE B1
	  ELSEIF(dt /= dtp)THEN
	   ok=.TRUE. ; EXIT B1 	  
	  ELSE  
	   SELECT CASE(langue)
	   CASE('english')
	    PRINT*,'no conv. in resout_rota, number of iterations :',
	1    compt	  
	   CASE DEFAULT
	    PRINT*,'pas de conv. dans resout_rota, nb. itérations :',
	1    compt
	   END SELECT
	   EXIT B1
	  ENDIF
	 ENDIF
	ENDDO B1
	dtp=dt

c	DO i=1,n_rot
c	 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,mrot(i),
c	1 l,bs,bd)
c	 WRITE(*,2000)mrot(i),bs
c	ENDDO
c	PAUSE'solution'	

c la base continue est désormais inutile		
	DEALLOCATE(a,b,indpc,indpc0)

c	PAUSE'sortie'
	
	RETURN
	
	END SUBROUTINE resout_rota
