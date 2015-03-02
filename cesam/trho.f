
c*********************************************************************    

	SUBROUTINE trho(nom,tau,teff,logg,t,dtdtau,dtdteff,dtdg,rhom4,
	1 drhom4dteff,drhom4dg,ftau,dftau,d2ftau)

c	routine private du module mod_atm
c	interpolation en (T, Rho) des lois T(tau) de Roger

c	lit un fichier apxx.dat ou amxx.dat et en tire la valeur
c	t(tau,teff,logg), de rho a tau=0.0001,
c	nommee rhom4, des derivees de t et rhom4 par rapport a teff et logg.

c	date de creation: 98 02 03
c	indice sur logtau: i
c	indice sur logg: j
c	indice sur teff: k

c	AUTEUR: R. CAYREL, adaptation a CESAM2k: P.Morel
c	A cause de la recherche de tau* ou T(tau*)=Teff, on ne peut pas
c	déborder des limites en teff et logg

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : langue
	USE mod_kind

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in) :: tau
      CHARACTER (len=*), INTENT(in) :: nom
      REAL (kind=dp), INTENT(inout) :: teff, logg
      REAL (kind=dp), INTENT(out) ::   t, dtdtau, dtdteff, dtdg, rhom4,
     1 drhom4dteff, drhom4dg, ftau, dftau, d2ftau

      REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ty
      REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: tlrhom4
      REAL (kind=dp), SAVE :: teffinit, logginit, stepteff, steplogg,
     1 teffmax, loggmax

      REAL (kind=dp) :: aux, dtdlogg, lrhom4, steplogtau, logtauinit,
     1 u, v, w, y, yu, yv, yw, logtau, dlrhom4dteff, dlrhom4dlogg,
     2 g, yprime, yseconde, yuu

      INTEGER, SAVE :: imax, jmax, kmax
      INTEGER :: i, j, k, kp, jp, ip

      CHARACTER (len=80) :: titre

      LOGICAL, SAVE :: init=.TRUE.
      LOGICAL :: ok     

c------------------------------------------------------------------------

c à la première entrée on lit le tableau des lois t(tau)
c	à toutes les autres on fait les interpolations nécessaires pour 
c	calculer t, rhom4 et leurs dérivées

	IF(init)THEN                       !premier appel
	 init=.FALSE.
	 INQUIRE(file=TRIM(nom),exist=ok)	 
	 IF(.NOT.ok)THEN
	  SELECT CASE(langue)	 
	  CASE('english')
	   WRITE(*,1005)TRIM(nom) ; WRITE(2,1005)TRIM(nom)
1005	   FORMAT('STOP, file of T(tau) law not found :',a)  	       
	  CASE DEFAULT		
	   WRITE(*,5)TRIM(nom) ; WRITE(2,5)TRIM(nom)
5	   FORMAT('ARRET, fichier de loi T(tau) non trouvé :',a)
	  END SELECT
	  STOP
	 ENDIF
	 OPEN(40,form='formatted',file=TRIM(nom),status='old')
	 READ(40,10)titre
10	 FORMAT(a)       
	 READ(40,*) teffinit,logginit,stepteff,steplogg,imax,jmax,kmax
c	 PRINT*, teffinit,logginit,stepteff,steplogg,imax, jmax, kmax
	 ALLOCATE(tlrhom4(jmax,kmax),ty(imax,jmax,kmax))
	 DO k=1,kmax
	  DO j=1,jmax
	   READ(40,11)tlrhom4(j,k)
11	   FORMAT(1x,f9.4)	 
	   READ(40,12) (ty(i,j,k), i=1,imax)
12	   FORMAT(1x,97(f8.5))	 
	  ENDDO
	 ENDDO
	 CLOSE(unit=40)

c on traite les cas aux limites des tables************

	 teffmax=teffinit+(kmax-1)*stepteff
	 loggmax=logginit+(jmax-1)*steplogg
	 tau_min=1.d-4
	ENDIF

	IF(teff > teffmax)THEN
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001)teff,teffmax ; WRITE(2,1001)teff,teffmax
1001	  FORMAT('STOP, call to the T(tau) law with Teff=',es10.3,
	1 ' > Teff max=',es10.3)	
	 CASE DEFAULT		
	  WRITE(*,1)teff,teffmax ; WRITE(2,1)teff,teffmax
1	  FORMAT('ARRET: appel à la loi T(tau) avec Teff=',es10.3,
	1 ' > Teff max=',es10.3)
	 END SELECT
	 STOP
	ELSEIF(teff < teffinit)THEN
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1002)teff,teffinit ; WRITE(2,1002)teff,teffinit
1002	  FORMAT('STOP, call to the T(tau) law with Teff=',es10.3,
	1 ' < Teff min=',es10.3)	
	 CASE DEFAULT		
	  WRITE(*,2)teff,teffinit ; WRITE(2,2)teff,teffinit
2	  FORMAT('ARRET, appel à la loi T(tau) avec Teff=',es10.3,
	1 ' < Teff min=',es10.3)
	 END SELECT
	 STOP
	ELSEIF(logg > loggmax)THEN
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1003)logg,loggmax ; WRITE(2,1003)logg,loggmax
1003	  FORMAT('STOP, call to the T(tau) law with Log(g)=',es10.3,
	1 ' > log(g)max =',es10.3)
	 CASE DEFAULT	
	  WRITE(*,3)logg,loggmax ; WRITE(2,3)logg,loggmax
3	  FORMAT('ARRET, appel à la loi T(tau) avec Log(g)=',es10.3,
	1 ' > log(g)max =',es10.3)
	 END SELECT	
	 STOP
	ELSEIF(logg < logginit)THEN
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1004)logg,logginit ; WRITE(2,1004)logg,logginit
1004	  FORMAT('STOP, call to the T(tau) law with Log(g)=',es10.3,
	1 ' < log(g)min =',es10.3)	
	 CASE DEFAULT		
	  WRITE(*,4)logg,logginit ; WRITE(2,4)logg,logginit
4	  FORMAT('ARRET, appel à la loi T(tau) avec Log(g)=',es10.3,
	1 ' < log(g)min =',es10.3)
	 END SELECT	
	 STOP
	ENDIF

c détermination du point pivot

	IF(teff > (teffmax-stepteff))THEN           ! bord du tableau
	 kp=kmax-2
	ELSEIF(teff < (teffinit+stepteff))THEN    ! bord du tableau
	 kp=2
	ELSE                                    ! cas general
	 aux=(teff-teffinit)/stepteff ; kp=NINT(aux)+1
	ENDIF
	IF(logg > (loggmax-steplogg))THEN
	 jp=jmax-2
	ELSEIF(logg < (logginit+steplogg)) THEN
	 jp=2
	ELSE
	 aux=(logg-logginit)/steplogg ; jp=NINT(aux)+1
	ENDIF
	logtau=LOG10(tau) ; logtauinit=-4. ; steplogtau=0.0625
	IF(logtau > (2.0-steplogtau)) THEN
	 ip=96
	ELSEIF(logtau < (logtauinit+steplogtau))THEN
	 ip=2
	ELSE
	 aux=(logtau-logtauinit)/steplogtau ; ip=NINT(aux)+1
	ENDIF

c fin calcul pivot
c	interpolation sur lrhom4(teff,logg)
c	PRINT*, ' coucou 2', ' jp =',jp,' kp= ',kp

	u=(logtau-(logtauinit+steplogtau*(ip-1)))/steplogtau
	v=(logg-(logginit+steplogg*(jp-1)))/steplogg
	w=(teff-(teffinit+stepteff*(kp-1)))/stepteff
c	PRINT*, ' coucou3 ',' u= ',u, ' v= ',v, ' w= ',w
	CALL quadr2(tlrhom4,v,w,jp,kp,lrhom4,yv,yw)
	rhom4=10.d0**lrhom4 ; g=10.d0**logg
	dlrhom4dteff=(yw/stepteff)
	drhom4dteff=dlrhom4dteff*2.3026d0*rhom4
	dlrhom4dlogg=(yv/steplogg) ; drhom4dg=dlrhom4dlogg*rhom4/g

c interpolation sur t 

	CALL quadr3(ty,u,v,w,ip,jp,kp,y,yu,yv,yw,yuu)

	t=teff*y			!PRINT*,'t=',t
	dtdteff=y+teff*yw/stepteff	!PRINT*,'dtdteff=',dtdteff
	dtdtau=(teff*yu/steplogtau)*0.43429d0/tau !PRINT*,'dtdtau=',dtdtau
	dtdlogg=teff*yv/steplogg ; dtdg=dtdlogg*0.4343d0/g
c	PRINT*,'dtdg=',dtdg

c variables supplémentaires f-tau, df_tau, d2f_tau***

	ftau=4.d0/3.d0*y**4 ; yprime=yu/steplogtau*0.43429d0/tau
	dftau=16.d0/3.d0*y**3*yprime
	yseconde=yuu*(0.43429d0/steplogtau/tau)**2-yprime/tau
	d2ftau=16.d0*y**2*(y*yseconde/3.d0+yprime**2)

	RETURN

	CONTAINS

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	SUBROUTINE quadr2(tableau,u,v,jp,kp,y,yu,yv)

c       interpolation quadratique d'une fonction donnee par tableau
c       au voisinage d'un point pivot predetermine (jp,jk)
c       donne aussi les derivees par rapport aux variables reduites
c       variables reelles divisees par le pas de la variable
c       auteur: R. CAYREL, F95: P. Morel

	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: tableau
	REAL (kind=dp), INTENT(in) :: u,v      
	INTEGER, INTENT(in) :: jp, kp
	REAL (kind=dp), INTENT(out) ::y, yu, yv             
	REAL (kind=dp) :: a1, a2, b11, b22, b12

	a1= (tableau(jp+1,kp)-tableau(jp-1,kp))/2.d0
	a2=(tableau(jp,kp+1)-tableau(jp,kp-1))/2.d0
	b11=(tableau(jp+1,kp)+tableau(jp-1,kp)-2.d0*tableau(jp,kp))/2.d0
	b22=(tableau(jp,kp+1)+tableau(jp,kp-1)-2.d0*tableau(jp,kp))/2.d0
	b12=(tableau(jp+1,kp+1)+tableau(jp-1,kp-1)-tableau(jp-1,kp+1)
	1 -tableau(jp+1,kp-1))/4.d0
	y=tableau(jp,kp)+a1*u+a2*v+b11*u**2+b22*v**2+b12*u*v
c	PRINT*,' a1=',a1,' a2=',a2,' b11= ',b11,' b22=',b22,' b12=',b22
	yu= a1+2.d0*b11*u+b12*v
	yv=a2+2.d0*b22*v+b12*u

	END SUBROUTINE quadr2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	SUBROUTINE quadr3(tableau,u,v,w,ip,jp,kp,y,yu,yv,yw,yuu)

c interpolation quadratique d'une fonction de trois variables donnee par 
c	un tableau. on donne aussi le "pivot" , point du tableau le plus
c	proche du point d'interpolation. le programme calcule aussi les
c	dérivées par rapport aux variables réduites (variables divisées
c	par le pas)

c	auteur: R. Cayrel
c	F95: P.Morel

	REAL (kind=dp), INTENT(in), DIMENSION(:,:,:) :: tableau
	REAL (kind=dp), INTENT(in) :: u,v,w    
	INTEGER, INTENT(in) :: ip, jp, kp
	REAL (kind=dp), INTENT(out) :: y, yu, yv, yw, yuu
	REAL (kind=dp) :: a1, a2, a3, b11, b22, b33, b12, b23, b13

	a1= (tableau(ip+1,jp,kp)-tableau(ip-1,jp,kp))/2.d0
	a2= (tableau(ip,jp+1,kp)-tableau(ip,jp-1,kp))/2.d0
	a3= (tableau(ip,jp,kp+1)-tableau(ip,jp,kp-1))/2.d0
	b11=(tableau(ip+1,jp,kp)+tableau(ip-1,jp,kp)
	1 -2.d0*tableau(ip,jp,kp))/2.d0
	b22=(tableau(ip,jp+1,kp)+tableau(ip,jp-1,kp)
	1 -2.d0*tableau(ip,jp,kp))/2.d0
	b33=(tableau(ip,jp,kp+1)+tableau(ip,jp,kp-1)
	1 -2.d0*tableau(ip,jp,kp))/2.d0
	b12=(tableau(ip+1,jp+1,kp)+tableau(ip-1,jp-1,kp)
	1 -tableau(ip-1,jp+1,kp)-tableau(ip+1,jp-1,kp))/4.d0
	b13=(tableau(ip+1,jp,kp+1)+tableau(ip-1,jp,kp-1)
	1 -tableau(ip-1,jp,kp+1)-tableau(ip+1,jp,kp-1))/4.d0
	b23=(tableau(ip,jp+1,kp+1)+tableau(ip,jp-1,kp-1)
	1 -tableau(ip,jp-1,kp+1)-tableau(ip,jp+1,kp-1))/4.d0
	y=tableau(ip,jp,kp)+a1*u+a2*v +a3*w+b11*u**2
	1 +b22*v**2+b33*w**2+b12*u*v+b13*u*w+b23*v*w
	yu=a1+2.d0*b11*u+b12*v+b13*w ; yv=a2+2.d0*b22*v+b12*u+b23*w
	yw=a3+2.d0*b33*w+b13*u+b23*v ; yuu=2.d0*b11

c	PRINT*,'a1=',a1,'a2=',a2,'a3=',a3
c	PRINT*,'b11=',b11,'b22=',b22,'b33=',b33
c	PRINT*,'b12=',b12,'b13=',b13,'b23=',b23

	END SUBROUTINE quadr3
 
	END SUBROUTINE trho
