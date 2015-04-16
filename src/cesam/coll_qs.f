
c***********************************************************************

	SUBROUTINE coll_qs(dt,compt,reprend,err,vare,qmax,corr)

c routine private du module mod_static
	
c résolution du système des équations de l'équilibre quasi-statique

c	Pour le probleme aux limites:
c	résolution du système d'équations différentielles non linéaires
c	de la structure interne par développement sur B-splines
c	par itération Newton avec dérivées analytiques

c	les équations de la SI sont écrites pour tout point dans
c	static_m/_r, le jacobien de la méthode NR est formé et résolu
c	dans resout dans l'espace des splines

c	pour le modèle quasi-statique, au temps t
c	 paramètres indépendants du temps
c	  m_qs : ordre des splines pour l'intégration par collocation 
c	  ord_qs=m_qs+r_qs=m_qs+1 : ordre de la représentation spline
c	  ne : nombre d'inconnues (6 ou 7)
c	 paramètres dépendants du temps
c	  nombre de couches au temps t : n_qs
c	  knot : nombre de points de table (nodaux),
c	  dim_qs=knot-ord_qs : dimension de l'espace des splines
c	  bp(ne,dim_qs) : coefficients des splines
c	  q(n_qs), qt(knot) : ABScisses, et vecteur nodal

c entrées:	
c	compt: compteur du nombre d'itérations
c	dt : pas temporel

c sorties:
c	err : erreur max
c	corr : facteur d'amortissement Newton-Raphson
c	qmax : indice de la couche au max d'erreur
c	vare : variable du max. erreur
c	reprend=.TRUE. : il y a un Pb. il faut diminuer dt

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c-------------------------------------------------------------------------

	USE mod_donnees, ONLY : en_masse, Ipg, langue, m_qs, ne, ord_qs,
	1 pturb, r_qs
	USE mod_kind
	USE mod_numerique, ONLY: bsp1dn, bvald, coll, gauss_band, linf,
	1 no_croiss
	USE mod_variables, ONLY: bp, dim_qs, knot, m23, n_qs, q, qt, r2, xl
	
	IMPLICIT NONE
		
	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in) :: compt		
	REAL (kind=dp), INTENT(out) :: err, corr
	INTEGER, INTENT(out) ::	qmax, vare
	LOGICAL, INTENT(out) :: reprend
		
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ae, dl, ds
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: a, b, y
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: d, temp	
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: be, f, dfdq
	
	REAL (kind=dp), PARAMETER :: nrm=0.15d0, err1=1.d-2
	REAL (kind=dp) :: er
	
	INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: indpc
	INTEGER, SAVE :: lq=1, nbj, ncoll, nl, nv=5, n_qsp
	INTEGER :: eq, cq, i, id, ind, indice, ip, j, ligne, spl, var
	 
	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: inversible
	
	CHARACTER (len=7), SAVE, ALLOCATABLE, DIMENSION(:) :: nom_qs
				
c---------------------------------------------------------------------

2000	FORMAT(8es10.3)

c----------------initialisations---------------------------

	IF(init)THEN
	 init=.FALSE.
	 
c noms des variables	 
	 IF(pturb)nv=6
	 ALLOCATE(nom_qs(nv))
	 SELECT CASE(nv)
	 CASE(6)
	  IF(en_masse)THEN
	   nom_qs=(/ 'ln Ptot',' ln T  ','  r^2  ',' l^2/3 ',' m^2/3 ',
	1  'ln Pgaz' /)
	  ELSE	  
	   nom_qs=(/ 'ln Ptot',' ln T  ','   r   ','   l   ','   m   ',
	1  'ln Pgaz' /)
	  ENDIF
	 CASE(5)
	  IF(en_masse)THEN
           nom_qs=(/ 'ln Pgaz',' ln T  ','  r^2  ',' l^2/3 ',' m^2/3 ' /)
	  ELSE	
	   nom_qs=(/ 'ln Pgaz',' ln T  ','   r   ','   l   ','   m   ' /)
	  ENDIF
	 END SELECT
	 
c définitions et allocations  diverses
	 ncoll=(n_qs-1)*m_qs	!nombre de points de collocation
	 nl=ne*(ncoll+r_qs)	!nombre de lignes du jacobien
	 nbj=ne*ord_qs		!longueur du bloc du jacobien
	 dim_qs=(n_qs-1)*m_qs+r_qs !dimension de l'espace des splines
	 	 
	 ALLOCATE(xl(ne),ae(ne,ne,0:r_qs),be(ne),a(nl,nbj),xcoll(ncoll),
	1 b(1,nl),indpc(nl),y(ne,0:r_qs),d(0:r_qs,ord_qs),
	2 ds(ord_qs,0:r_qs,ord_qs),dl(ne,0:r_qs,ord_qs),f(ne),dfdq(ne))
	 xl(1)=q(1)		!r au centre
	 xl(2)=q(1)		!l au centre
	 xl(3)=q(1)		!m au centre
	 xl(4)=q(n_qs)		!Ptot à l'extérieur
	 xl(5)=q(n_qs)		!T à l'extérieur
	 xl(6)=q(n_qs)		!m à l'extérieur
	 IF(pturb)xl(Ipg)=q(n_qs)	!avec Pturb, Pgaz à l'extérieur

	 CALL coll(q,n_qs,m_qs,xcoll)
c	 PRINT*,'xcoll' ; WRITE(*,2000)xcoll ; PAUSE'xcoll'

c les points de raccord q sont les entiers 1, 2,..n_qs,
c entre deux points de raccord les ABScisses des ord_qs
c points de collocation ne diffèrent que par un nombre entier 	
c les valeurs des ord_qs B-splines non id. nulles aux m_qs
c points de collocation sur [i,i+1] peuvent être tabulés une fois pour toutes	
	 DO i=1,m_qs
	  CALL linf(xcoll(i),qt,knot,lq)
	  CALL bvald(xcoll(i),qt,ord_qs,lq,r_qs,d) ; ds(i,:,:)=d
c	  PRINT*,i ; WRITE(*,2000)ds(i,:,:)
	 ENDDO
c	 PAUSE'ds'
	 
c idem pour les limites
	 DO i=1,ne
	  CALL linf(xl(i),qt,knot,lq)
	  CALL bvald(xl(i),qt,ord_qs,lq,r_qs-1,d) ; dl(i,:,:)=d
c	  PRINT*,i ; WRITE(*,2000)dl(i,:,:)
	 ENDDO	 
c	 PAUSE'dl'
	 
c d est désormais inutile 
	 DEALLOCATE(d)
	 
c quelques paramètres
c err1 : valeur minimale pour correction NR modifiée 
c err4 : correction max apres la 4-ieme itération NR
c nrm : coefficient pour NR modifiée
c iter_max0 : nb. max. d'itérations pour la conv. du modèle initial
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)nrm ; WRITE(2,1001)nrm
1001	  FORMAT('Use of the collocation method',/,
	1 'damping parameter for modified Newton-Raphson method, nrm=',
	2 es10.3)	 	
	 CASE DEFAULT
	  WRITE(*,1)nrm ; WRITE(2,1)nrm
1	  FORMAT('Utilisation de la méthode de collocation',/,
	1 'coefficient pour Newton-Raphson modifiée, nrm=',es10.3)
	 END SELECT

c nombre de couches du modèle précédent 
	 n_qsp=n_qs	 
	ENDIF
	
c------------------initialisations fin ------------------------------	
	 
c si le nombre de couches du modèle quasi-static a changé	
	IF(n_qs /= n_qsp)THEN
	 n_qsp=n_qs
	 	 
c redéfinitions et allocations diverses
	 ncoll=(n_qs-1)*m_qs	!nombre de points de collocation
	 nl=ne*(ncoll+r_qs)	!nombre de lignes du jacobien
	 dim_qs=(n_qs-1)*m_qs+r_qs	!DIMENSION de l'espace des splines
	 	 
	 DEALLOCATE(a,b,xcoll,indpc)
	 ALLOCATE(a(nl,nbj),xcoll(ncoll),b(1,nl),indpc(nl))
	
	 xl(4)=q(n_qs)		!Ptot à l'extérieur
	 xl(5)=q(n_qs)		!T à l'extérieur
	 xl(6)=q(n_qs)		!m à l'extérieur
	 IF(pturb)xl(Ipg)=q(n_qs)	!avec Pturb, Pgas à l'extérieur

	 CALL coll(q,n_qs,m_qs,xcoll)
c	 PRINT*,'xcoll' ; WRITE(*,2000)xcoll ; PAUSE'xcoll'	 
	ENDIF	
			
c initialisations
c ip est l'indice de la couche ou se trouvent les m_qs
c points de collocation, utilisé dans static_m/r pour le
c facteur de répartition fac(ip)
	ligne=0 ; lq=ord_qs ; a=0.d0 ; b=0.d0 ; ip=1
	 
c pour les points de collocation	 
	DO i=1,ncoll,m_qs
	 CALL linf(xcoll(i),q,n_qs,ip)
c	 WRITE(*,2000)REAL(i),REAL(ip),xcoll(i),q(ip),fac(ip)
	 DO j=1,m_qs
	  cq=i+j-1
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,xcoll(cq),lq,f,dfdq)
	  IF(no_croiss)PRINT*,'Pb. en 1 dans coll_qs'
	  y(:,0)=f ; y(:,1)=dfdq ; ind=ne*(lq-ord_qs)+1	
c	  PRINT*,i,j,cq,lq,ind ; WRITE(*,2000)xcoll(cq)
c	  WRITE(*,2000)y(:,0)
c	  WRITE(*,2000)y(:,1) ; PAUSE'y'
	  CALL static(1,cq,0,y,be,ae,compt,dt,reprend,ip) ; IF(reprend)RETURN	  
	  DO eq=1,ne
	   ligne=ligne+1 ; b(1,ligne)=be(eq); indpc(ligne)=ind
	   DO spl=1,ord_qs
	    DO var=1,ne
	     indice=ne*(spl-1)+var
	     DO id=0,r_qs
	      a(ligne,indice)=a(ligne,indice)+ae(eq,var,id)*ds(j,id,spl)
	     ENDDO	!id
	    ENDDO	!var
	   ENDDO	!spl
c	   WRITE(*,2000)a(ligne,:),b(1,ligne) ; PAUSE	   
	  ENDDO		!eq	   	   
	 ENDDO		!j	  
	ENDDO		!i
c	PAUSE'avant lim'
	  
c pour les limites, comme r_qs=1, il y a ne limites	   
	DO i=1,ne
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,xl(i),lq,f,dfdq)
         IF(no_croiss)PRINT*,'Pb de limite dans coll_qs'
	 y(:,0)=f
c	 PRINT*,i ; WRITE(*,2000)xl(i) ; WRITE(*,2000)y(:,0)
c	 WRITE(*,2000)y(:,1) ; PAUSE'lim'	
	 CALL static(2,cq,i,y,be,ae,compt,dt,reprend,ip) ; IF(reprend)RETURN 
	 ligne=ligne+1 ; b(1,ligne)=be(1) ; indpc(ligne)=ne*(lq-ord_qs)+1
	 DO spl=1,ord_qs
	  DO var=1,ne
	   indice=ne*(spl-1)+var
	   a(ligne,indice)=a(ligne,indice)+ae(1,var,0)*dl(i,0,spl)
	  ENDDO	!var
	 ENDDO		!spl
c	 PRINT*,ligne ; WRITE(*,2000)a(ligne,:),b(1,ligne) ; PAUSE'a'
	ENDDO		!i
c	PAUSE'apres lim'

c résolution du système linéaire
	
c	DO i=nl-15,nl
c	 PRINT*,'ligne',i ; WRITE(*,2000)a(i,:),b(1,i) ; PAUSE
c	ENDDO
c	PAUSE	
	 
c	PRINT*,ligne,nl,nbj ; PRINT*,indpc ; PAUSE'avant gauss'
	CALL gauss_band(a,b,indpc,nl,nl,nbj,1,inversible)
	IF(.not.inversible)THEN
	 PRINT*,'jacobien non inversible dans coll_qs' ; STOP
	ENDIF
c	DO i=1,ligne,ne
c	 WRITE(*,2000)(b(1,j),j=i,i+ne-1) ; PAUSE
c	ENDDO
c	PAUSE'apres gauss_band'

c construction du vecteur temporaire des corrections
	ALLOCATE(temp(ne,dim_qs)) ; temp=RESHAPE(b,SHAPE(temp))
	 	
c limitation des corrections
	corr=1.d0
	DO spl=1,dim_qs
	 B1: DO var=1,ne
	  IF(var == 7)CYCLE B1  
	  IF(spl == 1 .AND. (var > 2 .AND. var < 6))CYCLE B1
	  IF(ABS(bp(var,spl)) < err1)CYCLE B1
	  IF(corr*ABS(temp(var,spl)) > nrm*ABS(bp(var,spl)))corr=corr/2.d0
	 ENDDO B1	!iv
	ENDDO		!spl

c estimation de la précision
	err=0.d0
	DO spl=1,dim_qs
	 B2: DO var=1,ne
c	  IF(COUNT(var == iter_qs) == 1 )CYCLE B2
	  IF(COUNT(var == iter_qs) > 0 )CYCLE B2
	  
c au centre, on n'est pas exigeant pour R, L, M	  	  	 
c	  IF(spl <= ord_qs .AND. (var > 2 .AND. var < 6))CYCLE B2
	  IF(spl == 1 .AND. (var > 2 .AND. var < 6))CYCLE B2
	  
	  IF(ABS(bp(var,spl)) < 1.d0)CYCLE B2
	  er=ABS(temp(var,spl)/bp(var,spl))
	  IF(er > err)THEN
	   qmax=spl ; vare=var ; err=er
	  ENDIF	   
	 ENDDO B2	!var
	ENDDO		!spl

c corrections au maximum de l'erreur
	WRITE(usl_static,100)qmax
100	FORMAT(/,'Coeff. var. q-stat., corr. NR. au max. erreur, B-spline#',i5)
	IF(pturb)THEN
	 WRITE(usl_static,101)nom_qs(1:nv),bp(1:5,qmax),bp(7,qmax),
	1 temp(1:5,qmax),temp(7,qmax)
101	 FORMAT(6(1x,a7,2x),/,6es10.3,/,6es10.3)
	ELSE
	 WRITE(usl_static,103)nom_qs(1:nv),bp(1:5,qmax),temp(1:5,qmax)
103	 FORMAT(5(1x,a7,2x),/,5es10.3,/,5es10.3)
	ENDIF

cessai~~~~~~~~ pour éviter R^2, L^2/3,  M^2/3 négatifs
c	If(.FALSE.)THEN
c	If(.TRUE.)THEN
c limitation des corrections
c	 corr=1.d0 ; err=MIN(err,1.d0)
c	 DO spl=1,dim_qs
c	  B5: DO var=1,ne
c	   IF(var == 7)CYCLE B5
c	   IF(spl == 1 .AND. (var > 2 .AND. var < 6))CYCLE B5	  	  
c	   IF(ABS(bp(var,spl)-temp(var,spl)) < 1.d-4)CYCLE B5
c	   temp(var,spl)=SIGN(MIN(ABS(temp(var,spl)),0.6d0*ABS(bp(var,spl))),
c	1  temp(var,spl))
c	  ENDDO B5	!iv
c	 ENDDO		!spl
c	Endif
cessai~~~~~~~~~

c on considère qu'il y a stagnation si compt > iter_stag
c pour éviter des corrections consécutives égales et de signes opposés on impose
c corr < 0.6
c	IF(compt > iter_stag)corr=MIN(corr,0.7d0)
	IF(compt > iter_stag)corr=MIN(corr,0.51d0)	

c nouvelle solution, temp est désormais inutile
	bp=bp-temp*corr ; DEALLOCATE(temp)
c 	DO spl=1,10 !dim_qs
c	 WRITE(*,2000)bp(1:ne,spl) !; PAUSE
c 	ENDDO

c pour éviter des erreurs d'arrondi au centre, pour R, L, M
c	bp(3:5,1)=0.d0
	
c indice de la couche ou se trouve le maximum
	qmax=(qmax-ord_qs+1)/m_qs+1

c pour l'interpolation inverse m23 ou r2 ----> lnp, lnt, r2, m23, l23
c avec inter, extraction de r2 et m23 de la solution
c en_masse = .true. variables lagrangiennes m23=m**23, r2=r**2
c en_masse = .false. variables eulèriennes m23=m, r2=r
c pour lim_zc il faut les valeurs de r2 et m23 de la dernière itération		
	DO i=1,n_qs 
c	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lx,f,dfdx)
	 r2(i)=bp(3,m_qs*(i-1)+1) ; m23(i)=bp(5,m_qs*(i-1)+1)
	ENDDO
	
	RETURN

	END SUBROUTINE coll_qs
