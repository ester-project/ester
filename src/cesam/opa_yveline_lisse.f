
c*****************************************************************

	SUBROUTINE opa_yveline_lisse(xchim,t,ro,kap,dkapt,dkapro,dkapx)

c routine public du module mod_opa

c lissage des opacités Yveline par contour

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c entrées :
c	xchim : composition chimique / gramme
c	t : température
c	ro : densité

c sorties :
c	kap : opacité
c	dkapt : d kap / dt
c	dkapro : d kap / dro
c	dkapx : d kap /dx

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : langue, f_opa, ihe4, ln10, nchim,
	1 nom_chemin, z0
	USE mod_kind
	USE mod_numerique, ONLY : bval1, linf, noein, no_croiss
      
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out) ::  kap, dkapt, dkapro, dkapx

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: kap_tab
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: x_tab, z_tab ,
	1 t6_tab, lr_tab, lrtt, t6tt, xtt, ztt, qt, qr, qx, qz, dt, dr, dx, dz
	REAL (kind=dp), SAVE :: z_tab_lim
	REAL (kind=dp) :: lr, lr_use, ro_use, t6_use, t6, xi, x_use,
	1 zr, z_use

	INTEGER, SAVE :: nx, nz, nt, nr, mx, mr, mz, mt,
	1 l1, l2, l3, l4, knotx, knotz, knott, knotr
	INTEGER :: ir, jt, kx, lz, i, j, k, l 

	LOGICAL, SAVE :: init=.TRUE., lr_sort=.FALSE.,
	1 t6_sort=.FALSE., x_sort=.FALSE., z_sort=.FALSE.
	LOGICAL :: ok

	CHARACTER (len=120) :: nom_table
      
c--------------------------------------------------------------- 

2000	FORMAT(8es10.3)

c	PAUSE'entrée opa_yveline_lisse'

	IF(init)THEN
	 init=.FALSE.

c on lit lr_tab=Log10 R, t6_tab=T6, T6= T/1.d6, kap_tab=Log10 kap
	 nom_table=TRIM(nom_chemin)//f_opa(1)
	 INQUIRE(file=nom_table,exist=ok)
	 IF(.NOT.ok)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1008)nom_table ; WRITE(2,8)nom_table   
1008	   FORMAT('STOP, unknown opacity file: ',a)
	  CASE DEFAULT
	   WRITE(*,8)nom_table ; WRITE(2,8)nom_table   
8	   FORMAT('ARRET, fichier d''opacités inconnu: ',a)
	  END SELECT	  
	  STOP
	 ENDIF
	 
	 OPEN(unit=11,form='unformatted',status='unknown',
	1 file=TRIM(nom_chemin)//f_opa(1))
	 READ(11) nz,nx,nt,nr !; PRINT*, nz,nx,nt,nr ; PAUSE'nz etc..'
	 ALLOCATE(kap_tab(nz,nx,nt,nr),z_tab(nz),x_tab(nx),t6_tab(nt),
	1 lr_tab(nr))
	 READ(11)t6_tab ; READ(11)lr_tab !; CALL pause('t6, lr')
	 DO lz=1,nz
	  READ(11)z_tab(lz)
	  DO kx=1,nx
	   READ(11)x_tab(kx)
	   READ(11)((kap_tab(lz,kx,jt,ir),ir=1,nr),jt=1,nt)
	  ENDDO
	 ENDDO
	 CLOSE(unit=11) !; CALL pause('close')

c Log10 R ==> ln R, Log kap==> ln kap
	 lr_tab=lr_tab*ln10 ; kap_tab=kap_tab*ln10
 
c préparation des tables, m=2 pour interpolation linéaire,
c m=3 ou 4 pour lissage par contour
c	 mr=2 ; mt=2 ; mx=2 ; mz=2
	 mr=4 ; mt=4 ; mx=4 ; mz=4
 
c        PAUSE'avant allocate1'
       
	 ALLOCATE(lrtt(nr+mr),t6tt(nt+mt),xtt(nx+mx),ztt(nz+mz))

	 CALL noein(lr_tab,lrtt,nr,mr,knotr)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 1 dans opa_yveline_lisse' ; STOP
         ENDIF	 
	 CALL noein(t6_tab,t6tt,nt,mt,knott)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 2 dans opa_yveline_lisse' ; STOP
         ENDIF	  
	 CALL noein(x_tab, xtt ,nx,mx,knotx)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 3 dans opa_yveline_lisse' ; STOP
         ENDIF	 	 
	 CALL noein(z_tab, ztt ,nz,mz,knotz)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 4 dans opa_yveline_lisse' ; STOP
         ENDIF	 

	 l1=mt ; l2=mr ; l3=mx ; l4=mz

c on admet un petit dépassement en Z 
	 z_tab_lim=MIN(1.d0,z_tab(nz)*1.1d0)
       
c	 PAUSE'avant allocate2'

	 ALLOCATE(qz(mz+1),qx(mx+1),qt(mt+1),qr(mr+1),dz(mz),dx(mx),dt(mt),
	1 dr(mr))
	 SELECT CASE(langue)
	 CASE('english')	 
	  WRITE(*,1001)TRIM(f_opa(1)),nz,mz,nx,mx,nt,mt,nr,mr
	  WRITE(2,1001)TRIM(f_opa(1)),nz,mz,nx,mx,nt,mt,nr,mr	 
1001	  FORMAT('------Opacities Yveline------------',/,
	1 'm=2 linear interpolation, m > 2 smoothing by Béziers',/,
	2 'tables of the file : ',a,/,i3,
	3 ' number and order of B-splines in Z :',i3,/,i3,
	4 ' number and order of B-splines in X :',i3,/,i3,
	5 ' number and order of B-splines in T6 :',i3,/,i3,
	6 ' number and order of B-splines in Log R ::',i3,/,
	7 '--------------------------------------------------')	 
	 CASE DEFAULT	
	  WRITE(*,1)TRIM(f_opa(1)),nz,mz,nx,mx,nt,mt,nr,mr
	  WRITE(2,1)TRIM(f_opa(1)),nz,mz,nx,mx,nt,mt,nr,mr	 
1	  FORMAT('------Opacités Yveline------------',/,
	1 'interp. linéaire si m=2, lissage par Béziers si m > 2',/,
	2 'tables lues dans le fichier: ',a,/,i3,
	3 ' valeurs de Z, ordre des B-splines utilisé:',i3,/,i3,
	4 ' valeurs de X, ordre des B-splines utilisé:',i3,/,i3,
	5 ' valeurs de T6, ordre des B-splines utilisé:',i3,/,i3,
	6 ' valeurs de Log R, ordre des B-splines utilisé:',i3,/,
	7 '--------------------------------------------------')
	 END SELECT
	ENDIF
	
c calcul de Z
	IF(nchim > 1)THEN
	 zr=1.d0-SUM(xchim(1:ihe4))
	ELSE
	 zr=Z0
	ENDIF
	
c si Z > 0.1 appel à opa_opal2	      
	IF(zr > 0.098d0)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1007)zr,xchim(1) ; WRITE(2,1007)zr,xchim(1)
1007	  FORMAT('STOP, in opa_yveline_lisse overtaking of the limit',
	1 '0.1 > Z = ',es10.3,' with 0.09 < X =',es10.3,/,
	2 'use opa_opal2_co or opa_opal2_cno')	   
	 CASE DEFAULT
	  WRITE(*,7)zr,xchim(1) ; WRITE(2,7)zr,xchim(1)
7	  FORMAT('ARRET, dans opa_yveline_lisse dépassement de la limite',
	1 ' 0.1 > Z =',es10.3,' avec 0.09 < X =',es10.3,/,
	2 'utiliser opa_opal2_co ou opa_opal2_cno')
	 END SELECT
	 STOP 	 
	ENDIF	 

	t6=t*1.d-6 ; lr=LOG(ro/t6**3)
c	WRITE(*,2000)t,ro,t6,lr

	xi=abs(xchim(1))        !détermination de X

c on évite les sorties de table

c en composition chimique
	x_use=MAX(x_tab(1),MIN(xi,x_tab(nx)))	
	IF(.NOT.x_sort)THEN
	 x_sort=x_use /= xi
	 IF(x_sort)THEN
	  WRITE(*,2)xi,x_tab(1),x_tab(nx),xchim
	  WRITE(2,2)xi,x_tab(1),x_tab(nx),xchim
2	  FORMAT(/,'opa_yveline_lisse, au moins une sortie de table en X',
	1 /,'X=',es10.3,' en dehors de: [',es10.3,',',es10.3,'], Xchim=',/,
	2 (8es10.3))
	 ENDIF
	ENDIF

	IF(zr > z_tab_lim)THEN
	 WRITE(*,3)zr,z_tab(1),z_tab(nz),xchim
	 WRITE(2,3)zr,z_tab(1),z_tab(nz),xchim
	 WRITE(*,6) ; WRITE(2,6) ; STOP
6	 FORMAT('ARRET')
	ENDIF
	z_use=MAX(z_tab(1),MIN(zr,z_tab(nz)))	 	 		
	IF(.NOT.z_sort)THEN
	 z_sort=z_use /= zr
	 IF(z_sort)THEN
	  WRITE(*,3)zr,z_tab(1),z_tab(nz),xchim
	  WRITE(2,3)zr,z_tab(1),z_tab(nz),xchim
3	  FORMAT(/,'opa_yveline_lisse, au moins une sortie de table en Z',
	1 /,'Z=',es10.3,' en dehors de: [',es10.3,',',es10.3,'], Xchim=',/,
	2 (8es10.3))
	 ENDIF
	ENDIF
	
c en T et ro
	t6_use=MAX(t6_tab(1),MIN(t6,t6_tab(nt)))	
	IF(.NOT.t6_sort)THEN
	 t6_sort=t6_use /= t6
	 IF(t6_sort)THEN
	  WRITE(*,4)t6,t6_tab(1),t6_tab(nt),ro,t
	  WRITE(2,4)t6,t6_tab(1),t6_tab(nt),ro,t
4	  FORMAT(/,'opa_yveline_lisse, au moins une sortie de table en T6',
	1 /,'T6=',es10.3,' en dehors de: [',es10.3,',',es10.3,']' ,/,
	2 'ro=',es10.3,', t=',es10.3)
	 ENDIF
	ENDIF

	lr_use=MAX(lr_tab(1),MIN(lr,lr_tab(nr)))
	ro_use=EXP(lr_use)*t6_use**3	
	IF(.NOT.lr_sort)THEN
	 lr_sort=lr_use /= lr
	 IF(lr_sort)THEN
	  WRITE(*,5)lr,lr_tab(1),lr_tab(nr),ro,ro_use,t
	  WRITE(2,5)lr,lr_tab(1),lr_tab(nr),ro,ro_use,t
5	  FORMAT(/,'opa_yveline_lisse, au moins une sortie de table en ro',
	1 /,'lr=',es10.3,' en dehors de: [',es10.3,',',es10.3,']' ,/,
	2 'ro=',es10.3,', ro_use=',es10.3,', t=',es10.3)
	 ENDIF
	ENDIF

c recherche des indices
	CALL linf(lr_use,lrtt,knotr,l1) ; CALL linf(t6_use,t6tt,knott,l2)
	CALL linf(x_use,xtt, knotx,l3) ; CALL linf(z_use,ztt, knotz,l4)
c	PRINT*,l1,l2,l3,l4

c interpolation
	CALL bval1(lr_use,lrtt,mr,l1,qr,dr)
	CALL bval1(t6_use,t6tt,mt,l2,qt,dt)
	CALL bval1(x_use, xtt,mx,l3,qx,dx)
	CALL bval1(z_use, ztt,mz,l4,qz,dz)

	kap=0.d0 ; dkapt=0.d0 ; dkapro=0.d0 ; dkapx=0.d0
	DO i=1,mr
	 ir=l1-mr+i
	 DO j=1,mt
	  jt=l2-mt+j
	  DO k=1,mx
	   kx=l3-mx+k
	   DO l=1,mz
	    lz=l4-mz+l
c	    PRINT*,ir,jt,kx,lz
	    kap=kap      +  kap_tab(lz,kx,jt,ir)*qr(i)*qt(j)*qx(k)*qz(l)
	    dkapro=dkapro+  kap_tab(lz,kx,jt,ir)*dr(i)*qt(j)*qx(k)*qz(l)
	    dkapt=dkapt  +  kap_tab(lz,kx,jt,ir)*qr(i)*dt(j)*qx(k)*qz(l)
	    dkapx=dkapx  +  kap_tab(lz,kx,jt,ir)*qr(i)*qt(j)*dx(k)*qz(l)
	   ENDDO
	  ENDDO
	 ENDDO
	ENDDO

c tranformation inverse (log R, T6)-->(T , ro)
	kap=EXP(kap) ; dkapt=kap*1.d-6*(dkapt-3.d0/t6_use*dkapro)
	dkapro=kap/ro_use*dkapro ; dkapx=kap*dkapx

c	PAUSE'sortie opa_yveline_lisse'

	RETURN

	END SUBROUTINE opa_yveline_lisse
