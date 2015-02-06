
c***********************************************************************

      SUBROUTINE opa_yveline(xchim,t,ro,kappa,dkapdt,dkapdr,dkapdx)

c     routine public du module mod_opa

c     calcul de l'opacité repris du code de Genève
c     polynômes de Lagrange pour LOG T6 et LOG R
c     quadratique pour X et Z
c     Yveline (aout 1991-code de Genève--> 10/95 CESAM)   

c     Daniel (octobre 1996-version pour CESAM 2.x --> CESAM 3.1)

c     P. Morel (2 avril 97) suppression des j1,j2,j3,j4 (risque de
c     difficultés avec j2=0), modif de la logique, réécriture de zr

c     Yveline (décembre 99- calcul de l'indice dans pos_table_op)
c     P.Morel 25 04 00 : évite des sorties de table

c     Adaptation: P.Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k

c-------------------------------------------------------------------

c entrées :
c      xchim(1)=X : comp. chim. par gramme
c      t : température K
c      ro : densité cgs
  
c sorties :
c      kappa : opacité gr / cm2)
c      dkapdt : kappa / d t
c      dkapdr : kappa / d densité              
c      dkapdx : kappa / d xchim(1)
  
c      Z est obtenu par 1-X-Y

c-------------------------------------------------------------------

	USE mod_donnees, ONLY : f_opa, ihe4, langue, ln10, nchim, z0
	USE mod_kind
	USE mod_variables, ONLY : sortie
      
	IMPLICIT NONE
   
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out) :: dkapdt, dkapdr, dkapdx, kappa
	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: vlk
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: vlr, vlt6, vx, vz	

	REAL (kind=dp) :: lr, t6, x, z	
	
	INTEGER, SAVE :: nz, nx, nt, nr

	LOGICAL, SAVE :: init=.TRUE., pass=.TRUE.
	CHARACTER (len=80) :: nom_chemin = "/data1/sdeheuve/local/src/cesam2k_v1.1.8_ESTA/SUN_STAR_DATA/"

c------------------------------------------------------------------------
  
2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 IF(z0 <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1001)z0 ; WRITE(2,1001)z0
1001	   FORMAT('STOP, in opa_yveline 0 > Z0=',es10.3)	  	  
	  CASE DEFAULT	
           WRITE(*,1)z0 ; WRITE(2,1)z0
1	   FORMAT('ERREUR, dans opa_yveline 0 > Z0=',es10.3)
	  END SELECT
	  CALL sortie
	 ENDIF
	ENDIF

c	si Z > 0.1 appel à opa_opal2

	IF(nchim > 1)THEN
	 z=1.d0-SUM(xchim(1:ihe4))
	ELSE
	 z=z0
	ENDIF
	IF(z > 0.098d0)THEN
	 IF(xchim(1) > 9.d-2)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1003)z,xchim(1) ; WRITE(2,1003)z,xchim(1)
1003	   FORMAT('STOP, in opa_yveline overtaking of the limit',
	1   '0.1 > Z = ',es10.3,' with 0.09 < X =',es10.3)	   
	  CASE DEFAULT
	   WRITE(*,3)z,xchim(1) ; WRITE(2,3)z,xchim(1)
3	   FORMAT('ARRET, dans opa_yveline dépassement de la limite',
	1   ' 0.1 > Z =',es10.3,' avec 0.09 < X =',es10.3)
	  END SELECT
	  CALL sortie 	 
	 ENDIF	
	 IF(pass)THEN
	  pass=.FALSE.
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1002)z,t,ro ; WRITE(2,1002)z,t,ro
1002	   FORMAT('call to opa_opal2, in opa_yveline because Z=',es10.3,/,
	1  'T='es10.3,', ro=',es10.3)	  	  
	  CASE DEFAULT	
           WRITE(*,2)z,t,ro ; WRITE(2,2)z,t,ro
2	   FORMAT('appel à opa_opal2, dans opa_yveline car Z=',es10.3,/,
	1  'T='es10.3,', ro=',es10.3)
	  END SELECT
	 ENDIF
	 CALL opa_opal2(xchim,t,ro,kappa,dkapdt,dkapdr,dkapdx,.FALSE.)
	 RETURN
	ELSE
	 t6=t*1.d-6 ; lr=LOG10(ro/t6**3) ; x=xchim(1)   
	 CALL kappa_opal(lr,t6,ro,x,z,kappa,dkapdr,dkapdt,dkapdx)	
	ENDIF

	RETURN

	CONTAINS

c************************************************************************

      SUBROUTINE kappa_opal(lr,t6,ro,x,z,cap,capr,capt6,capx)

c     sous-programme de calcul de l'opacité radiative par interpolation
c     dans les tables de Livermore. Ces tables donnent le log10 de
c     l'opacite pour differentes valeurs de z, x, T6(T/1e6) et
c     log R (R=ro/T6^3). Les valeurs tabulees de l'opacite correspondent
c     a l'opacite radiative on entre en parametre nz, nx, nt, nr qui
c     sont respectivement le nombre de valeurs de z, de x, de T et de R.

c     parametres d'entree de kappa_opal:
c     t6 et log R     
c     x et z: H et metaux dans la couche consideree
c     sortie: kappa, dlogkappa/dt6,dlogkappa/dlogR, dlogkappa/dx    
    
c              Yveline (aout 1991-code de Genève--> 10/95 CESAM)
c              Daniel (octobre 1996-version pour CESAM 2.x --> CESAM 3.1)
c            P.Morel 24 04 00 : modifs pour éviter les sorties de table

c----------------------------------------------------------------------

      REAL (kind=dp), INTENT(in) :: lr, t6, ro, x, z
      REAL (kind=dp), INTENT(out) :: cap, capr, capt6, capx

      REAL (kind=dp), DIMENSION(3,3) :: opa, opar, opat6
      REAL (kind=dp), DIMENSION(3) :: bidon, opax, oparx, opaxx, qx,
     1 qk, qkk, qkl, opat6x
      REAL (kind=dp), SAVE :: z_tab_lim     
      REAL (kind=dp) :: kapq, dkapq, lr_use, lt6, lt6_use, x_use, z_use

      INTEGER, DIMENSION(3,3) :: ilag, ilin
      INTEGER :: iz, iza, izb, ix, ixa, ixb, it, ita, ir, ira, i, j,
     1 irl, itl, irla, itla, jdum, idum

      LOGICAL, SAVE :: i_lect=.TRUE., i_exit=.FALSE., r_sort=.FALSE.,
     1 t_sort=.FALSE., x_sort=.FALSE., z_sort=.FALSE.

c-------------------------------------------------------------------- 

2000  FORMAT(8es10.3)

c     ----- lecture des tables d'opacité ------------------------------

	IF(i_lect)THEN
	 i_lect=.FALSE. ; CALL lect_opal

c on admet un petit dépassement en Z	 
	 z_tab_lim=MIN(1.d0,vz(nz)*1.1d0)
	 
c appel fictif à kappa_cond pour écritures (pas de calcul)
	 bidon(1)=x ; bidon(2)=1.d0-x-z ; bidon(3)=z
	 CALL kappa_cond(bidon,1.d0,1.d-6,kappa,dkapdt,dkapdr,dkapdx)
	ENDIF

c	localisation dans les tables du point pour lequel on cherche
c	l'opacité. Sortie éventuelle des tables.

c	en composition chimique

	x_use=MAX(vx(1),MIN(x,vx(nx)))	
	IF(.NOT.x_sort)THEN
	 x_sort=x_use /= x
	 IF(x_sort)THEN
	  WRITE(*,2)x,vx(1),vx(nx),xchim ; WRITE(2,2)x,vx(1),vx(nx),xchim
2	  FORMAT('opa_yveline, au moins une sortie de table en X',/,
	1 'X=',es10.3,' en dehors de: [',es10.3,',',es10.3,'], Xchim=',/,
	2 (8es10.3))
	 ENDIF
	ENDIF
	CALL pos_table_op(x_use,vx,nx,jdum,ix,idum)
      	
	IF(z > z_tab_lim)THEN
	 WRITE(*,3)z,vz(1),vz(nz),xchim ; WRITE(2,3)z,vz(1),vz(nz),xchim
	 WRITE(*,6) ; WRITE(2,6) ; STOP
6	 FORMAT('ARRET')
	ENDIF	
	z_use=MAX(vz(1),MIN(z,vz(nz)))		
	IF(.NOT.z_sort)THEN
	 z_sort=z_use /= z
	 IF(z_sort)THEN
	  WRITE(*,3)z,vz(1),vz(nz),xchim ; WRITE(2,3)z,vz(1),vz(nz),xchim
3	  FORMAT('opa_yveline, au moins une sortie de table en Z',/,
	1 'Z=',es10.3,' en dehors de: [',es10.3,',',es10.3,'], Xchim=',/,
	2 (8es10.3))
	 ENDIF
	ENDIF
	CALL pos_table_op(z_use,vz,nz,jdum,iz,idum)

c     ---------- en densité et température ----------------------------

c	Pour éviter les sorties de table (PM mars 2003)

	lt6=LOG10(t6)
	lt6_use=MAX(vlt6(1),MIN(lt6,vlt6(nt)))
	IF(.NOT.t_sort)THEN
	 t_sort=lt6_use /= lt6
	 IF(t_sort)THEN
	  WRITE(*,4)lt6,vlt6(1),vlt6(nt),ro,t,x,z
	  WRITE(2,4)lt6,vlt6(1),vlt6(nt),ro,t,x,z
4	  FORMAT('opa_yveline, au moins une sortie de table en LOG10 T6',
	1 /,'LOG10 T6=',es10.3,' en dehors de: [',es10.3,',',es10.3,']',/,
	2 'ro=',es10.3,', t=',es10.3,', x=',es10.3,', z=',es10.3)
	 ENDIF
	ENDIF
	CALL pos_table_op(lt6_use,vlt6,nt,itl,idum,itla) 

	lr_use=MAX(vlr(1),MIN(lr,vlr(nr)))
	IF(.NOT.r_sort)THEN
	 r_sort=lr_use /= lr
	 IF(r_sort)THEN
	  WRITE(*,5)lr,vlr(1),vlr(nr),ro,t,x,z
	  WRITE(2,5)lr,vlr(1),vlr(nr),ro,t,x,z !; PAUSE'sortie'
5	  FORMAT('opa_yveline, au moins une sortie de table en LOG10 R',
	1 /,'LOG10 R=',es10.3,' en dehors de: [',es10.3,',',es10.3,']',/,
	2 'ro=',es10.3,', t=',es10.3,', x=',es10.3,', z=',es10.3)
	 ENDIF
	ENDIF
	CALL pos_table_op(lr_use,vlr,nr,irl,idum,irla)

c     In this version of the opacity tables, a 4-point double argument
c     interpolation is incorporated. This requires 16 given opacity
c     entries, and the input pair (rh,T) should be preferently in the
c     center square. The ideal configuration required is the following
c     (horizontal axis:log R, vertical axis: T6, *: entry given in the
c     table, X: position of the input (log R,T6)

c        *   *   *   *
c
c        *   *   *   *
c              x
c        *   *   *   *
c
c        *   *   *   *

c      the 16-point configuration is ONLY used if all corresponding
c      entries are provided. If a 16-point configuration is impossible
c      a linear interpolation is made if the 4 surrounding points have
c      their opacity entries, otherwise error EXIT occurs. In the case
c      of very low densities, however, extrapolation is done by
c      constant continuation. In this way spurious increases of the
c      opacity due to 3-rd degree extrapolation are suppressed.

c     ----- Peut-on faire une interpolation Lagrangienne ou seulement
c     une interpolation linéaire? -------------------------------------

c     ---------- test pour l'interpolation lagrangienne ---------------

      ir=irla ; it=itla
      IF(ir >= 3 .AND. ir <= nr-1 .AND. it >= 3 .AND. it <= nt-1)THEN
       izb=0
       DO iza=iz-1,iz+1
        izb=izb+1 ; ixb=0
        DO ixa=ix-1,ix+1
         ixb=ixb+1 ; ilag(izb,ixb)=1
         b1: DO ita=it-2,it+1
          DO ira=ir-2,ir+1
           IF(vlk(iza,ixa,ita,ira) > 98.d0)THEN
            ilag(izb,ixb)=0 ; EXIT b1
           ENDIF
          ENDDO
         ENDDO b1
        ENDDO
       ENDDO
      ELSE
       ilag=0
      ENDIF

c     ---------- sinon test pour l'interpolation linéaire -------------

      ir=irl ; it=itl ; izb=0
      DO iza=iz-1,iz+1
       izb=izb+1 ; ixb=0
       DO ixa=ix-1,ix+1
        ixb=ixb+1
        IF(ilag(izb,ixb)==0)THEN
         b3: DO ita=it-1,it
          b2: DO ira=ir-1,ir
           IF(vlk(iza,ixa,ita,ira) > 98.d0)THEN
            i_exit=.TRUE.
            IF(ira==ir .AND. abs(lr-vlr(ir-1)) < 1e-15
     1       .AND. vlk(iza,ixa,ita,ir-1) <= 98.d0)THEN
             ir=ir-1 ; i_exit=.FALSE. ;  CYCLE b2
            ENDIF
           ENDIF
          ENDDO b2
          IF(i_exit)THEN
           IF(it /=2 .AND. ita == it .AND.
     1       ABS(lt6_use-vlt6(it-1)) < 1e-15 .AND.
     2       vlk(iza,ixa,it-1,ira) <= 98.d0)THEN
            it=it-1 ; i_exit=.FALSE. ; CYCLE b3
           ELSE
            WRITE(*,20)ro,t6*1.d-6 ; WRITE(2,20)ro,t6*1.d-6 ; STOP
20          FORMAT(/,'ARRET: sortie de table en ro T, ro=',es10.3,
	1   ' t=',es10.3)	    
           ENDIF
          ENDIF
         ENDDO b3
         ilin(izb,ixb)=1
        ENDIF
       ENDDO
      ENDDO

c     ----- calcul des opacités pour les points voisins ---------------
c
      izb=0
      DO iza=iz-1,iz+1
       izb=izb+1 ; ixb=0
       DO ixa=ix-1,ix+1
        ixb=ixb+1
        IF(ilag(izb,ixb)==1)THEN
         ir=irla ; it=itla
         CALL intl_opal(iza,ixa,ir,it,lr_use,lt6_use,opa(izb,ixb),
     1   opar(izb,ixb),opat6(izb,ixb))
        ELSE
         IF(ilin(izb,ixb)==1)THEN
          ir=irl ; it=itl
          CALL intlin_opal(vlr(ir-1),vlr(ir),
     1    vlt6(it-1),vlt6(it),vlk(iza,ixa,it-1,ir-1),
     3    vlk(iza,ixa,it-1,ir),vlk(iza,ixa,it,ir-1),vlk(iza,ixa,it,ir),
     5    lr_use,lt6_use,opa(izb,ixb),opar(izb,ixb),opat6(izb,ixb))
         ELSE
          PRINT*,'pb opacités: ilin=ilag=0' ; STOP
         ENDIF
        ENDIF
        opa(izb,ixb)=10.d0**opa(izb,ixb)
        opat6(izb,ixb)=opa(izb,ixb)/t6/1d6*(-3.d0*opar(izb,ixb)+
     1  opat6(izb,ixb))
        opar(izb,ixb)=opa(izb,ixb)/ro*opar(izb,ixb)
       ENDDO
      ENDDO

c     interpolation quadratique (K, X)

      DO i=1,3
       qx(i)=vx(ix+i-2)
      ENDDO

      DO i=1,3
       DO j=1,3
        qk(j)=opa(i,j) ; qkk(j)=opar(i,j) ; qkl(j)=opat6(i,j)
       ENDDO
       CALL sub_quad(x,qx,qk,kapq,dkapq) ; opax(i)=kapq ; opaxx(i)=dkapq
       CALL sub_quad(x,qx,qkk,kapq,dkapq) ; oparx(i)=kapq
       CALL sub_quad(x,qx,qkl,kapq,dkapq) ; opat6x(i)=kapq
      ENDDO

c     interpolation quadratique (K, Z)

      DO i=1,3
       qx(i)=vz(iz+i-2)
      ENDDO

      CALL sub_quad(z,qx,opax,cap,dkapq)
      CALL sub_quad(z,qx,opaxx,capx,dkapq)
      CALL sub_quad(z,qx,oparx,capr,dkapq)
      CALL sub_quad(z,qx,opat6x,capt6,dkapq)
      
      RETURN

      END SUBROUTINE kappa_opal

c------------------------------------------------------------------------

      SUBROUTINE intlin_opal(R1,R2,T1,T2,vk11,vk12,vk21,vk22,lr,
     1 T6,opac,opacr,opact6)

      REAL (kind=dp) ,INTENT(in) :: R1, R2, T1, T2, vk11, vk12, vk21, vk22,
     1 lr, T6
      REAL (kind=dp) ,INTENT(out) :: opac, opacr, opact6  

      REAL (kind=dp) :: dr, dt, x, y, dr1, dr2, dt1, ct1, ct2

      dr=R2-R1 ; dt=T2-T1 ; x=lr-R1 ; y=T6-T1 ; dr1=(vk12-vk11)/dr 
      dr2=(vk22-vk21)/dr ; ct1=vk11+x*dr1 ; ct2=vk21+x*dr2
      dt1=(ct2-ct1)/dt ; opac=ct1+y*dt1
      opacr=dr1+(dr2-dr1)*(y/dt) ; opact6=dt1
      
      RETURN

      END SUBROUTINE intlin_opal
C
c------------------------------------------------------------------------
C
      SUBROUTINE intl_opal(iza,ixa,ir,it,lr,T6,opac,opacr,opact6)

c     lagrangian interpolation of opac in (log R,T6) within a table

      REAL (kind=dp), INTENT(in) :: lr,T6
      INTEGER, INTENT(in) :: iza, ixa, ir, it
      REAL (kind=dp), INTENT(out) :: opac, opacr, opact6

c      INTEGER, PARAMETER :: pnz=16, pnx=9, pnt=85, pnr=23
c      REAL (kind=dp), DIMENSION(pnz,pnx,pnt,pnr) :: vlk
c      REAL (kind=dp), DIMENSION(pnx) :: vx
c      REAL (kind=dp), DIMENSION(pnt) :: vlt6, vt6
c      REAL (kind=dp), DIMENSION(pnr) :: vlr
c      INTEGER :: nz, nx, nt, nr
c      COMMON/val_opal/vz, vx, vlt6, vlr, vt6, vlk, nt, nr, nz, nx

      REAL (kind=dp), DIMENSION(4):: xr, xt, s1, s2, pr, ppr, pt, ppt

      INTEGER i, k

      DO i=1,4
       xr(i)=vlr(ir-3+i) ; xt(i)=vlt6(it-3+i)
      ENDDO

      CALL lpol_op(xr,4,lr,pr,ppr) ; CALL lpol_op(xt,4,T6,pt,ppt)

      s1=0.d0 ; s2=0.d0
      DO k=1,4
       DO i=1,4
        s1(k)= s1(k)+pr(i)*vlk(iza,ixa,it-3+k,ir-3+i)
        s2(k)= s2(k)+ppr(i)*vlk(iza,ixa,it-3+k,ir-3+i)
       ENDDO
      ENDDO

      opac=0.d0 ; opacr=0.d0 ; opact6=0.d0
      DO k=1,4
       opac=opac+pt(k)*s1(k) ; opacr=opacr+pt(k)*s2(k)
       opact6=opact6+ppt(k)*s1(k)
      ENDDO
      
      RETURN

      END SUBROUTINE intl_opal

c------------------------------------------------------------------------

      SUBROUTINE sub_quad(x,qx,qk,k,dk)

c     interpolation quadratique

      REAL (kind=dp), DIMENSION(:), INTENT(in) :: qx, qk
      REAL (kind=dp), INTENT(in) :: x
      REAL (kind=dp), INTENT(out) :: k, dk

      REAL (kind=dp) :: x1, x2, x3, y1, y2, y3, denom, aa, bb, cc

      x1=qx(1) ; x2=qx(2) ; x3=qx(3) ; y1=qk(1) ; y2=qk(2) ;  y3=qk(3)
      denom=(x2-x1)*(x3-x1)*(x2-x3) ; aa=(y2-y1)*(x3-x1)-(y3-y1)*(x2-x1)
      aa=aa/denom ; bb=(y3-y1)/(x3-x1)-aa*(x3+x1) ; cc=y1-aa*x1*x1-bb*x1
      k=aa*x*x+bb*x+cc ; dk=2.*aa*x+bb
      
      RETURN

      END SUBROUTINE sub_quad

c----------------------------------------------------------------------

	SUBROUTINE lect_opal

c	sous-programme de lecture des tables d'opacité de Livermore
c	les sorties sont vz=z, vx=x, vlt6=log T6, vlrro=log R et
c	vlk=log kappa (log decimaux). on affecte la valeur 99.
c	aux points pour lesquels l'opacité n'est pas donnée.

c	INTEGER, PARAMETER :: pnz=16, pnx=9, pnt=85, pnr=23
c	REAL (kind=dp), DIMENSION(pnz,pnx,pnt,pnr) :: vlk
c	REAL (kind=dp), DIMENSION(pnz) :: vz
c	REAL (kind=dp), DIMENSION(pnx) :: vx
c	REAL (kind=dp), DIMENSION(pnt) :: vlt6, vt6
c	REAL (kind=dp), DIMENSION(pnr) :: vlr
c	INTEGER :: nz, nx, nt, nr
c	COMMON/val_opal/vz,vx,vlt6,vlr,vt6,vlk,nt,nr,nz,nx

	INTEGER :: jz, jx, jt, jr

	LOGICAL :: ok

	CHARACTER (len=120) :: nom_table

c	ouverture des fichiers d'opacité

	nom_table=TRIM(nom_chemin)//f_opa(1)
	INQUIRE(file=nom_table,exist=ok)
	IF(ok)THEN
	 OPEN(unit=11,form='unformatted',status='old',file=nom_table)
	ELSE
	 WRITE(*,"('fichier d''opacités inconnu: ',a50)")nom_table
	 WRITE(*,*)'ARRET' ; STOP
	ENDIF

c	lecture des tables	

	READ(11) nz, nx, nt, nr	
	WRITE(*,1)TRIM(f_opa(1)),nz,nx,nt,nr
	WRITE(2,1)TRIM(f_opa(1)),nz,nx,nt,nr	
1	FORMAT(/,'--------------opacités Yveline---------------',//,
	1 'tables d''opacités Yveline: ',a,/,i3,' valeurs de Z,',
	2 i3,' valeurs de X,',i3,' valeurs de T6,',i3,' valeurs de LOG R')


	ALLOCATE(vlr(nr),vlt6(nt),vlk(nz,nx,nt,nr),vx(nx),vz(nz))

	READ(11)vlt6 ; READ(11)vlr
	DO jz=1,nz
	 READ(11) vz(jz)
	 DO jx=1,nx
	  READ(11)vx(jx) ;  READ(11)((vlk(jz,jx,jt,jr),jr=1,nr),jt=1,nt)
	 ENDDO
	ENDDO

	vlt6=LOG10(vlt6)

c	fermeture des fichiers d'opacité

	CLOSE(unit=11)
	
	RETURN

	END SUBROUTINE lect_opal

c------------------------------------------------------------------------

      SUBROUTINE lpol_op(x,n,x0,p,pp)

c     polynomes de Lagrange (4 points) et leurs derivees (n=4)
c     x=coordonnees tabulees
c     x0=coordonnee pour laquelle p et pp sont demandes
c     p=valeur des polynomes en x0
c     pp=valeur des derivees de ces polynomes en x0

      REAL (kind=dp), DIMENSION(:), INTENT(in) :: x
      REAL (kind=dp), INTENT(in) :: x0
      INTEGER, INTENT(in) :: n

      REAL (kind=dp), DIMENSION(:), INTENT(inout) :: p, pp

          REAL (kind=dp), DIMENSION(4) :: s

      INTEGER :: i, j, k, l

      DO j=1,n
       p(j)=1.d0 ; pp(j)=0.d0 ; s(j)=1.d0
       DO i=1,n
        IF(i/=j)THEN
         p(j)=p(j)*(x0-x(i))/(x(j)-x(i)) ; s(j)=s(j)*(x(j)-x(i))
        ENDIF
       ENDDO
       DO  k=1,n
        IF(k /= j)THEN
         DO l=1,n
          IF(l /= j)THEN
           IF(l > k)THEN
            pp(j)=pp(j)+(x0-x(k))*(x0-x(l))
           ENDIF
          ENDIF
         ENDDO
        ENDIF
       ENDDO
       pp(j)=pp(j)/s(j)
      ENDDO
      
      RETURN

      END SUBROUTINE lpol_op

c------------------------------------------------------------------------

      SUBROUTINE pos_table_op(val,tab_val,nval,ilin,iquad,ilag)

c     modification du positionnement suivant que l'on veut interpoler
c     linéairement, quadratiquement ou Lagrangien Yveline : 99 12 20
c     repérage dans une table
      
      REAL (kind=dp), INTENT(in), DIMENSION(:) :: tab_val
      REAL (kind=dp), INTENT(in) :: val
      INTEGER, INTENT(in) :: nval

      INTEGER, INTENT(out) :: ilin, iquad, ilag

      REAL (kind=dp) :: t12

      INTEGER :: i

      ilin=0 ; ilag=0 ; iquad=0
      IF((tab_val(2)) > val)THEN
       ilin=2 ; iquad=2 ; ilag=0 ; RETURN
      ENDIF

      IF((tab_val(nval-1)) < val)THEN
       ilin=nval ; iquad=nval-1 ; ilag=0 ; RETURN
      ENDIF

      DO i=3,nval
       t12=(tab_val(i-1)+tab_val(i))/2.d0
       IF(tab_val(i) >= val)THEN
        ilin=i ; ilag=i
        IF(t12 > val)THEN
         iquad=i-1
        ELSE
         iquad=i
        ENDIF
        RETURN
       ENDIF
      ENDDO
      
      RETURN

      END SUBROUTINE pos_table_op

	END SUBROUTINE opa_yveline
