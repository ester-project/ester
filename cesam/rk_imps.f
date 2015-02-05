
c**************************************************************************

      SUBROUTINE rk_imps(t_t,ro_t,compx,t,ro,compy,dt,esti,ok,n_pt,
     1 z_vent,dm)

c     routine private du module mod_evol

c     intégration de la composition chimique par divers schémas RKi
c     implémentation selon Hairer et Wanner p.128
c     contrairement a rk_imp, il n'y a pas de test de précision
c     par comparaison avec un second schéma
c     on interpole linéairement T et ro pour les points intermédiaires

c     Auteur: P.Morel, Département J.D. Cassini, O.C.A.,
c     CESAM2k

c     on a noté z les zi de Hairer et Wanner (zi, charges des élem. chim.
c     est variable publique du module mod_donnees)

c entrees
c     t_t(k), t(k) : températures aux temps t et t+dt, couche k
c     ro_t(k), ro(k) : densités aux temps t et t+dt, couche k
c     compx(nchim,k) : composition chimique au temps t couche k
c     dt: pas temporel
c     dm : épaisseur en masse des couches des ZC
c     n_pt : nombre de couches à mélanger

c entrees / sorties
c     compy : composition chimique au temps t+dt
c     ok = .TRUE. : nouveau modèle / la précision est atteinte

c sorties
c     esti : précision estimée par élément

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ab_min, dtmin, ihe4, langue,
	1 nchim, nom_elem, ordre, precix
	USE mod_kind
	USE mod_nuc, ONLY : nuc, vent
	USE mod_numerique, ONLY : gauss_band, neville
      
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: compx 
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: dm, ro, ro_t, t, t_t
	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in) :: n_pt
	LOGICAL, INTENT(in) :: z_vent     
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: compy
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: esti
	LOGICAL, INTENT(out) :: ok

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: a, corr, corr1,
	1 jacob, z
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: roi, ti
	REAL (kind=dp), PARAMETER, DIMENSION(0:1) :: x=(/ 0.d0, 1.d0 /)
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: c   
     
	REAL (kind=dp), DIMENSION(nchim,nchim) :: jac, jacs 
	REAL (kind=dp), DIMENSION(nchim) :: compx0, dcomp, dcomps, ex, scale 	
	REAL (kind=dp), DIMENSION(0:1) :: f    
	REAL (kind=dp), DIMENSION(5) :: e 
     
	REAL (kind=dp) :: be7, bid, b8, corr_max, ero, f17, et, hh, n13, o15

	INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: indpc
	INTEGER, SAVE :: ns, ord_intp, s, tourmax
	INTEGER :: i, ii, imax, j, jj, l, m, k, tour

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: force

c---------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(9es8.1)
2002	FORMAT(12es8.1)
2003	FORMAT(12es8.1)

      IF(init)THEN
       init=.FALSE. ; PRINT*

       SELECT CASE(ordre)
       CASE(1)
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'integration using IRK1 of order 1'	
	CASE DEFAULT
         PRINT*,'intégration par IRK1  d''ordre 1'
	END SELECT	
        s=1       !p pour précision
        ALLOCATE(a(s,s),c(s))
        a(1,1)=1.d0 ; c(1)=1.d0

       CASE(2)
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'integration using Lobatto IIIc of order 2'	
	CASE DEFAULT       
         PRINT*,'intégration par Lobatto IIIc d''ordre 2'
        END SELECT	
        s=2     
        ALLOCATE(a(s,s),c(s))
        a(1,1)=0.5d0 ; a(1,2)=-0.5d0
        a(2,1)=0.5d0 ; a(2,2)=0.5d0
        c(1)=0.d0 ; c(2)=1.d0

       CASE(3)
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'integration using Radau IIa of order 3'	
	CASE DEFAULT       
        PRINT*,'intégration par Radau IIa d''ordre 3'
       END SELECT	
        s=2
        ALLOCATE(a(s,s),c(s))
        a(1,1)=5.d0/12.d0 ; a(1,2)=-1.d0/12.d0
        a(2,1)=0.75d0     ; a(2,2)=0.25d0
        c(1)=1.d0/3.d0 ; c(2)=1.d0

       CASE(4)
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'integration using Lobatto IIIc of order 4'	
	CASE DEFAULT      
         PRINT*,'intégration par Lobatto IIIc d''ordre 4'
        END SELECT	
        s=3
        ALLOCATE(a(s,s),c(s))
        a(1,1)=1.d0/6.d0 ; a(1,2)=-1.d0/3.d0 ; a(1,3)=1.d0/6.d0
        a(2,1)=1.d0/6.d0 ; a(2,2)=5.d0/12.d0 ; a(2,3)=-1.d0/12.d0
        a(3,1)=1.d0/6.d0 ; a(3,2)=2.d0/3.d0  ; a(3,3)=1.d0/6.d0     
        c(1)=0.d0 ; c(2)=0.5d0 ; c(3)=1.d0

       CASE(5)
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'integration using Radau IIa of order 5'		
	CASE DEFAULT   
         PRINT*,'intégration par Radau IIa d''ordre 5'
        END SELECT	
        s=3
        ALLOCATE(a(s,s),c(s))
        a(1,1)=(88.d0-7.d0*SQRT(6.d0))/360.d0
        a(2,1)=(296.d0+169.d0*SQRT(6.d0))/1800.d0
        a(3,1)=(16.d0-SQRT(6.d0))/36.d0
        a(1,2)=(296.d0-169.d0*SQRT(6.d0))/1800.d0
        a(2,2)=(88.d0+7.d0*SQRT(6.d0))/360.d0
        a(3,2)=(16.d0+SQRT(6.d0))/36.d0
        a(1,3)=(-2.d0+3.d0*SQRT(6.d0))/225.d0
        a(2,3)=(-2.d0-3.d0*SQRT(6.d0))/225.d0
        a(3,3)=1.d0/9.d0
        c(1)=(4.d0-SQRT(6.d0))/10.d0
        c(2)=(4.d0+SQRT(6.d0))/10.d0
        c(3)=1.d0

       CASE(6)
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'integration using Lobatto IIIc of order 6'		
	CASE DEFAULT       
         PRINT*,'intégration par Lobatto IIIc d''ordre 6'
        END SELECT	
        s=4
        ALLOCATE(a(s,s),c(s))
        a(1,1)=1.d0/12.d0
        a(2,1)=1.d0/12.d0
        a(3,1)=1.d0/12.d0
        a(4,1)=1.d0/12.d0
        a(1,2)=-SQRT(5.d0)/12.d0
        a(2,2)=1.d0/4.d0
        a(3,2)=(10.d0+7.d0*SQRT(5.d0))/60.d0
        a(4,2)=5.d0/12.d0
        a(1,3)=SQRT(5.d0)/12.d0
        a(2,3)=(10.d0-7.d0*SQRT(5.d0))/60.d0
        a(3,3)=1.d0/4.d0
        a(4,3)=5.d0/12.d0
        a(1,4)=-1.d0/12.d0
        a(2,4)=SQRT(5.d0)/60.d0
        a(3,4)=-SQRT(5.d0)/60.d0
        a(4,4)=1.d0/12.d0
        c(1)=0.d0
        c(2)=(5.d0-SQRT(5.d0))/10.d0
        c(3)=(5.d0+SQRT(5.d0))/10.d0
        c(4)=1.d0
       CASE DEFAULT
        SELECT CASE(langue)
	CASE('english')
         PRINT*,'STOP, there is no formulae for order=',ordre	
	CASE DEFAULT      
         PRINT*,'ARRET, il n''y a pas de formule pour ordre=',ordre
        END SELECT
        STOP
       END SELECT

c      allocations, ns: dimension du vecteur Z, du système N-R

       ns=nchim*s

c      définitions, initialisations

       ord_intp=1     !ordre d'interpolation : linéaire
       tourmax=30 !nombre de tours max pour N-R
c      scale=ab_min(1:nchim)*100.d0  !abondances de référence
 
       ALLOCATE(corr(nchim,s),corr1(1,ns),indpc(ns),jacob(ns,ns),z(nchim,s)) 
       indpc=1

c      WRITE(*,2000)ab_min ; WRITE(*,2000)scale
c      CALL pause('ab_min, scale')

      ENDIF
      
c     les valeurs moyennes dépendent du nombre de couches
      
	ALLOCATE(ti(n_pt,s),roi(n_pt,s))

c     calcul des valeurs moyennes des abondances au temps t   
c     l'ordre des arguments de compx(nchim,n_pt) permet le produit de
c     matrice avec dm sur 1,n_pt

      compx0(:)=MATMUL(compx,dm(1:n_pt)) ; scale=abs(compx0)
c     PRINT*,'n_pt,dm',n_pt ; WRITE(*,2000)dm(1:n_pt)
c     WRITE(*,2000)compx0(1:nchim) ; CALL pause('matmul')    

c     interpolation des températures et densités  
c     les valeurs de T et de ro aux points intermédiaires par interp.
c     ti(n_pt,s)=ti(point de ZM,etape RK), idem pour roi

	DO i=1,s
	 DO k=1,n_pt
	  f(0)=t(k) ; f(1)=t_t(k) ; ti(k,i)=neville(c(i),x,f,ord_intp)
c	  PRINT*,i,c(i),k,ti(k,i)
	 ENDDO
c	 WRITE(*,2000)ti(1:n_pt,i)
	ENDDO
	DO i=1,s
	 DO k=1,n_pt
	  f(0)=ro(k) ; f(1)=ro_t(k) ; roi(k,i)=neville(c(i),x,f,ord_intp)
c	  PRINT*,i,c(i),k,roi(k,i)
	 ENDDO
c	 PRINT*,'roi',i ; WRITE(*,2000)roi(:,i)
	ENDDO
	
c initialisations
c       vecteur Z: z(ns)=z(nchim,s), z(variable,etape RK)
c       Z=T(Z1,..,Znchim), Zi=T(z_i1, z_i2, z_is), i=1,nchim
c       pour la variable=k, etape=l, i=nchim*(l-1)+k ==> z(i)

	z=0.d0

c iterations NR
	tour=0 ; corr_max=1.d5
	B2: DO

c        jacob(ns,ns)=jacob(nchim,s,nchim,s)
c             =jacob(variable,etape,variable,etape)
c        pour la variable=k, etape=l, i=ns*(l-1)+k
c        derivée/variable=m, etape=q, j=ns*(q-1)+m ==> jacob(i,j)
 
	 jacob=0.d0  
	 DO i=1,ns          !initialisations
	  jacob(i,i)=1.d0       !1 dans la diagonale du jacobien
	 ENDDO

c        on met Z dans corr= Z - A*fi*h
	 corr=z

c        formation du sytème
c        a la 1-ière itération Z=0, ensuite Z=Z_s de l'itération précédente
c        ce qui permet d'updater compy, utilise dans l'appel à nuc

	 DO j=1,s           	!pour chaque étape R-K      
c	  PRINT*,'compx0' ; WRITE(*,2000)compx0
c	  PRINT*,'z(:,j)',j ; WRITE(*,2000)z(:,j)
       
	  compy=compx0+z(:,j)     !Vt: compy(l) = X0(l)+z_lj
	
c	  PRINT*,'compy' ; WRITE(*,2000)compy
c         PRINT*,'j,z',j,(z(1:nblem,j)

	  dcomps=0.d0 ; jacs=0.d0   !initialisation

	  DO k=1,n_pt   !pour chaque couche de la ZC

c	   F ( x+Ci h,y0+z )	 
	   CALL nuc(ti(k,j),roi(k,j),compy,dcomp,jac,.TRUE.,2,
	1  e,et,ero,ex,hh,be7,b8,n13,o15,f17) 
       
c s'il y a perte/gain de masse
c	   IF(z_vent)CALL vent(compy,dcomp,jac) 
	   
	   dcomps=dcomps+dcomp*dm(k)*dt !somme du second membre sur la ZC
	   jacs=jacs+jac*dm(k)*dt   !somme sur la ZC du jacobien   
	  ENDDO     !k

c         contribution au jacobien et au second membre : Z - A*fi*h
	  DO l=1,nchim     !chaque element
	   DO i=1,s     !chaque etape R_K

c           second membre corr=Z-A*fi*h, Z est deja dans corr
	    corr(l,i)=corr(l,i)-a(i,j)*dcomps(l) 

c           PRINT*,'l,i,j,corr',l,i,j,corr(l,i)

c           contribution au jacobien il y a déjà 1 dans la diag de jacob
	    DO m=1,nchim
c	     indice=nchim*(s*(nchim*(j-1)+m-1)+i-1)+l
	     ii=nchim*(i-1)+l ; jj=nchim*(j-1)+m       
	     jacob(ii,jj)=jacob(ii,jj)-jacs(l,m)*a(i,j)
	    ENDDO       !m chaque derivee/element 
	   ENDDO        !i chaque etape R_K
	  ENDDO         !l chaque element
	 ENDDO          !j chaque etape R_K

c	 DO i=1,ns
c	  WRITE(*,2000)jacob(i,:)
c	 ENDDO
c	 CALL pause('avant inversion')

c	 pour la solution on met tout à la queue leu-leu
c	 il y a peut être moyen de faire plus élégant 
       
	 DO i=1,s
	  DO l=1,nchim
	   corr1(1,nchim*(i-1)+l)=corr(l,i)
	  ENDDO
	 ENDDO

	 CALL gauss_band(jacob,corr1,indpc,ns,ns,ns,1,ok)
	 DO i=1,s
	  DO l=1,nchim
	   corr(l,i)=corr1(1,nchim*(i-1)+l)
	  ENDDO
	 ENDDO

	 IF(.NOT.ok)THEN
	  PRINT*,'dans rk_imps, jacobien non inversible'
	  PRINT*,'t_t' ; WRITE(*,2000)t_t(1:n_pt) ; PRINT*,'ro_t' 
	  WRITE(*,2000)ro_t(1:n_pt) ; PRINT*,'dt' ; WRITE(*,2000)dt
	  PRINT*,'n_pt,tour,ns',n_pt,tour,ns ; PRINT*,'dm'
	  WRITE(*,2000)dm(1:n_pt) ; PRINT*,'t' ; WRITE(*,2000)t(1:n_pt)
	  PRINT*,'ro' ; WRITE(*,2000)ro(1:n_pt) ; EXIT B2
	 ENDIF

c	 PRINT*,'correction',tour ; WRITE(*,2000)corr(1,1:s)
c	 CALL pause('iteration')

	 corr_max=0.d0
	 B1: DO l=1,nchim
	  IF(ABS(compx0(l)) < ab_min(l))CYCLE B1 
	  DO i=1,s
	   bid=ABS(corr(l,i)/scale(l))
	   corr_max=MAX(corr_max,bid) ; IF(corr_max == bid)imax=l
	  ENDDO
	 ENDDO B1

c	 correction sur Z

c	 DO i=1,s
c	  WRITE(*,2002)z(:,i)
c	 ENDDO

	 z=z-corr

c	 DO i=1,s
c	  WRITE(*,2002)corr(:,i)
c	 ENDDO   
c        PRINT*,'corr_max',corr_max! ; CALL pause('corr_max')

	 tour=tour+1
	 IF(tour > tourmax)THEN
	  IF(corr_max <= 1.d-2 .AND. imax /= ihe4)THEN
	   ok=.TRUE. ; force=.TRUE. ; EXIT B2
	  ELSE      !arret
	   compy=compx0+z(:,s) ; esti=ABS(z(:,s))/scale
	   WRITE(*,118)ordre,nom_elem(imax),t(1),ro(1),dt,corr_max
	   WRITE(2,118)ordre,nom_elem(imax),t(1),ro(1),dt,corr_max
118	   FORMAT('pas de conv. de NR dans rk_imps, ordre:',i3,
	1  ', element: ',a,/,'t=',es10.3,', ro=',es10.3,', dt=',
	2  es10.3,', corr_max=',es10.3)
	   WRITE(*,116) ; WRITE(2,116)
	   WRITE(*,117)compx(1:MIN(12,nchim),1)
	   WRITE(*,117)compx(13:nchim,1) 
	   WRITE(*,117)compy(1:MIN(12,nchim)); WRITE(*,117)compy(13:nchim)
	   WRITE(*,117)esti(1:MIN(12,nchim)); WRITE(*,117)esti(13:nchim)
	   WRITE(2,117)compx(1:MIN(12,nchim),1)
	   WRITE(2,117)compx(13:nchim,1) 
	   WRITE(2,117)compy(1:MIN(12,nchim)); WRITE(2,117)compy(13:nchim)
	   WRITE(2,117)esti(1:MIN(12,nchim)); WRITE(2,117)esti(13:nchim)
	   ok=.FALSE.
c	   CALL pause('non convergence dans rk_imps')
	   IF(dt/2.d0 > dtmin)EXIT B2
	   PRINT*,'le dt < dtmin, abandon' ; STOP
	  ENDIF
	  EXIT B2
	 ELSEIF(corr_max < precix)THEN
	  ok=.TRUE. ; force= .FALSE.; EXIT B2
	 ENDIF
	ENDDO B2    !boucle infinie NR

c     PRINT*,'tour',tour,ordre,ordre ; WRITE(*,2000),t(1),ro(1)
c     CALL pause('tour')

c     solution et estimation des variations relatives
c     la solution est en Z_s

      IF(ok)THEN

       compy=compx0+z(:,s) ; esti=ABS(z(:,s))/scale

c	IF(z_vent)THEN
c	 WRITE(*,2000)compx(1:nchim,1) ; WRITE(*,2000)compy(1:nchim)
c	 PAUSE'solution'
c	ENDIF
c      PRINT*,'tour,ordre,s',tour,ordre,s
c      PRINT*,'sortie rk_imps, n_pt, ok/compx/compy',n_pt,ok
c      WRITE(*,2001)compx0 ; WRITE(*,2001)compy
c      WRITE(*,2001)corr(1,1:ns) ; WRITE(*,2001)scale 
c      CALL pause('rk_imps')
c      PRINT*,'sortie rk_imps, n_pt, ok/compx/compy',n_pt,ok
c      WRITE(*,2003)compx0 ; WRITE(*,2003)compy ; WRITE(*,2003)esti
c      CALL pause('solution')

       IF(force)THEN
        WRITE(*,115)ordre,nom_elem(imax),t(1),ro(1),dt,corr_max
115     FORMAT(/,'on force la  conv. de NR dans rk_imps, ordre:',i3,
     1  ', element: ',a3,/,'t=',es10.3,', ro=',es10.3,', dt=',
     2  es10.3,', corr_max=',es10.3)
        WRITE(*,116)
116     FORMAT('x / y / esti')
        WRITE(*,117)compx(1:MIN(12,nchim),1)
	WRITE(*,117)compx(13:nchim,1)
        WRITE(*,117)compy(1:MIN(12,nchim)) ; WRITE(*,117)compy(13:nchim)
        WRITE(*,117)esti(1:MIN(12,nchim)); WRITE(*,117)esti(13:nchim)
117     FORMAT(12es10.3)
       ENDIF
      ENDIF

      DEALLOCATE(roi,ti)
      
      RETURN

      END SUBROUTINE rk_imps
 
