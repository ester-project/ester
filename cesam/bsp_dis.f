
c******************************************************************

	SUBROUTINE bsp_dis(n,x,f,nd,id,fd,eps,nx,m,xt,knot,contour)

c	subroutine public du module mod_numerique

c	interpolation avec discontinuités
c	formation de la base nodale, calcul des coefficients des splines
c	s'exploite avec bsp1dn

c	avec contour=.TRUE. lissage par spline, la base est utilisée comme
c	base duale, on obtient un lissage

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	version: 23 01 04
c	SPLINES/SOURCE95

c entrées
c	n : nombre de tables
c	x(nx) : abscisses
c	nd : nombre de discontinuités
c	id(0:nd+1) : commence en 0 finit en nd+1!!!
c	indices des discontinuités, la première en 1, la dernière en nd
c	fd(n,nd): valeurs des fonctions à droite des discontinuités
c	eps : écart à droite, parameter définit avant l'appel à bsp_dis
c	m : ordre des B-splines
c	nx : nombre de points
c	contour=.TRUE. : lissage par spline, il y a formation du vecteur
c	nodal et introduction des discontinuités dans f
c	contour est OPTIONAL

c entrée/sortie

c	f(n,nx+nd) : fonction à interpoler f(n,nx)/coeff. des splines

c sorties
c	xt(knot) : points de table
c	knot : nombre de points de table, knot=nx+nd+m

c-----------------------------------------------------------------------

	USE mod_kind
      
	IMPLICIT NONE

	REAL(kind=dp), INTENT(in), DIMENSION(:,:) :: fd
	REAL(kind=dp), INTENT(in), DIMENSION(:) :: x
	REAL(kind=dp), INTENT(in) :: eps
	INTEGER, INTENT(in) :: n, m, nd, nx
	LOGICAL, INTENT(in), OPTIONAL :: contour
	REAL(kind=dp), INTENT(inout), DIMENSION(:,:) :: f
	REAL(kind=dp), INTENT(inout), DIMENSION(:) :: xt
	INTEGER, INTENT(inout), DIMENSION(0:) :: id
	INTEGER, INTENT(inout) :: knot
	
	REAL(kind=dp), DIMENSION(n,nx+nd) :: s	
	REAL(kind=dp), DIMENSION(nx+nd,m) :: ax
	REAL(kind=dp), DIMENSION(nx+nd) :: xc
	REAL(kind=dp), DIMENSION(m+1) :: qx
	REAL(kind=dp), DIMENSION(n) :: dfxdx, fx

	INTEGER, DIMENSION(nx+nd) :: indpc
	INTEGER :: lx, i, ij, nc

	LOGICAL :: cont, inversible

c---------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	WRITE(*,2000)x(1:nx) ; WRITE(*,2000)fd(1:nd)

c	si contour est présent et est .TRUE. lissage par contour

	IF(PRESENT(contour))THEN
	 cont=contour
	ELSE
	 cont=.FALSE.
	ENDIF

c	si nd <= 0 pas de discontinuité, appel à bsp1dn ou seulement
c	formation du vecteur nodal si contour=.TRUE.

	IF(nd <= 0)THEN
	 IF(cont)THEN
          CALL noein(x,xt,n,m,knot)
	  IF(no_croiss)THEN
	   PRINT*,'Arrêt dans bsp_dis après appel à noein' ; RETURN
	  ENDIF	  
	 ELSE
	  CALL bsp1dn(n,f,x,xt,nx,m,knot,.FALSE.,x(1),lx,fx,dfxdx)
	  IF(no_croiss)THEN
	   PRINT*,'Arrêt dans bsp_dis après appel à bsp1dn' ; RETURN
	  ENDIF
	 ENDIF

	ELSE 	!avec discontinuités
      
c	 PRINT*,'bsp_dis ',nd,id(0:nd+1)    
       
c	 le vecteur nodal avec discontinuités

	 CALL noeu_dis(eps,id,knot,m,nd,nx,x,xt)
	 IF(no_croiss)THEN
	  PRINT*,'Dans bsp_dis après appel à noeu_dis' ; RETURN
	 ENDIF
	 
c	 PRINT*,nx,nd,knot,id(0:nd+1)
c	 WRITE(*,2000)xt(1:knot) ; PRINT*,'eps',eps ; PAUSE

	 IF(SIZE(f,1) /= n)THEN
	  PRINT*,'dans bsp_dis, la 1-iere dim. de f=',SIZE(f,1),' /=',n
	  PRINT*,'ARRET' ; STOP
	 ELSEIF(SIZE(f,2) /= nx+nd)THEN
	  PRINT*,'dans bsp_dis, la 2-de dim. de f=',SIZE(f,2),' /=',nx+nd
	  PRINT*,'ARRET' ; STOP  
	 ENDIF

	 nc=0   !nc: nombre de points de donnée nx+nd
	 ij=1   !indice de la discontinuité
	 i=1    !indice de la couche
	 DO WHILE(i <= nx)
	  nc=nc+1 ; xc(nc)=x(i)
	  s(:,nc)=f(:,i)    !s : VT, hors discontinuité

c	  sur une discontinuité: petit écart à droite

	  IF(i == id(ij))THEN	   	   
	   nc=nc+1 ; xc(nc)=x(i)+eps*(x(i+1)-x(i))
	   s(:,nc)=fd(:,ij) !sur discontinuités
c	   PRINT*,ij,id(ij),fd(1,ij)
	   ij=MIN(ij+1,nd)
	  ENDIF
	  i=i+1
	 ENDDO      !i
       
	 IF(cont)THEN
	  f=s ; RETURN
	 ENDIF  

c	 PRINT*,'les f, nc=',nc,knot

c	 DO i=1,nx
c	  WRITE(*,2000)(f(1:n,i)
c	 ENDDO  !i
c	 WRITE(*,2000)xc ; WRITE(*,2000)x ; PAUSE'nc'

	 lx=m
	 DO i=1,nc
	  CALL linf(xc(i),xt,knot,lx) ; CALL bval0(xc(i),xt,m,lx,qx)
	  ax(i,:)=qx(1:m) ; indpc(i)=lx-m+1
	 ENDDO

	 CALL gauss_band(ax,s,indpc,nc,nc,m,n,inversible)
	 IF(.NOT.inversible)THEN
	  PRINT*,'matrice non inversible dans sbsp_dis, ARRET'
	  no_croiss=.TRUE. ; RETURN	  
	 ENDIF
	 f=s
	ENDIF

	RETURN  

	END SUBROUTINE bsp_dis
