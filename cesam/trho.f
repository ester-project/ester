
c********************************************************************

      SUBROUTINE trho(nom,tau,teff,logg,t,dtdtau,dtdteff,dtdg,rhom4,
     1 drhom4dteff,drhom4dg,ftau,dftau,d2ftau,table)
C
C     routine PRIVATE du module mod_atm
c     interpolation en (T, Rho) des lois T(tau) de Roger
C
c     lit un fichier apxx.dat ou amxx.dat et en tire la valeur
c     t(tau,teff,logg), de rho a tau=0.0001,
c     nommee rhom4, des derivees de t et rhom4 par rapport a teff et logg.
C
c     indice sur logtau: i
c     indice sur logg: j
c     indice sur teff: k
C
c     AUTEUR: R. CAYREL, adaptation a CESAM2k: P.Morel, B.Pichon
C     aussi version compatible avec MARCS (B.Pichon)
c     A cause de la recherche de tau* ou T(tau*)=Teff, on ne peut pas
c     déborder des limites en teff et logg on se place sur elles
C
c---------------------------------------------------------------------
C
      USE mod_donnees, ONLY : langue
      USE mod_kind
C
      IMPLICIT NONE
C
      REAL (kind=dp), INTENT(in) :: tau, teff, logg
      CHARACTER (len=*), INTENT(in) :: nom
      Character(Len=5), Intent(In) :: table      
            
      REAL (kind=dp), INTENT(out) ::   t, dtdtau, dtdteff, dtdg, rhom4,
     1 drhom4dteff, drhom4dg, ftau, dftau, d2ftau
C
      REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ty
      REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: tlrhom4
      REAL (kind=dp), SAVE :: teffinit, logginit, stepteff, steplogg,
     1 teffmax, loggmax, logtaumax, steplogtau, logtauinit
C
      REAL (kind=dp) :: dtdlogg, lrhom4, 
     1 u, v, w, y, yu, yv, yw, logtau, dlrhom4dteff, dlrhom4dlogg,
     2 g, yprime, yseconde, yuu, teff_lu, logg_lu, teff_calc, logg_calc
C
      INTEGER, SAVE :: imax, jmax, kmax
      INTEGER :: i, j, k, kp, jp, ip, ios, pos, dimlogtau
C
      CHARACTER (len=80) :: titre 
      !                            1234567890123
      Character(Len=13) :: form = "(1x,xy(f8.5))"
C
      LOGICAL, SAVE :: init=.TRUE., passG=.TRUE., passT=.TRUE.     
      LOGICAL :: ok, outsideG, outsideT  
C
c------------------------------------------------------------------------
C
c à la première entrée on lit le tableau des lois t(tau)
c     à toutes les autres on fait les interpolations nécessaires pour 
c     calculer t, rhom4 et leurs dérivées
C
C tau_min et une variable privee de mod_atm
C
      If ( init ) Then                       !premier appel
         init=.FALSE.
C
         SELECT CASE ( table )
            CASE ( "roger" )
               pos = INDEX(nom,"fesh")
               If ( pos == 0 ) STOP " TRHO PB 0 1 nom de tdetau erronné"
               logtauinit = -4.0
               logtaumax  = +2.0
               steplogtau = 0.0625
            CASE ( "MARCS" ) ; 
               pos = INDEX(nom,"MARCS_Z")
               If ( pos == 0 ) STOP " TRHO PB 0 2 nom de tdetau erronné"
               logtauinit = -5.0
               logtaumax  = +2.0
               steplogtau = 0.1000
            CASE DEFAULT 
               STOP " TRHO : PB 0 0 nom de tdetau erronné"
         END SELECT
C
         INQUIRE(file=TRIM(nom),exist=ok)
         If (.NOT.ok) Then
            SELECT CASE(langue)	 
               CASE('english') ; WRITE(*,1005) TRIM(nom) ; WRITE (2,1005)TRIM(nom)
               CASE DEFAULT    ; WRITE(*,5)TRIM(nom) ; WRITE(2,5) TRIM(nom)
            END SELECT
            STOP
1005        FORMAT('STOP, file of T(tau) law not found :',a)  	       
5           FORMAT('ARRET, fichier de loi T(tau) non trouvé :',a)
         End If
C
         OPEN(40,File=TRIM(nom),Form='Formatted',Status='old',Action="Read",IOStat=ios)
         IF ( ios /= 0 ) STOP " TRHO : PB 0 3 "
C
         READ(40,"(A)") titre
C
C        Write(11,"(1X,F6.1,T10,F3.1,T15,F5.1,2X,F3.1,T30,I2,T35,I2,T40,I2)") Teff_min, Logg_min, Teff_step, Logg_step, dim_Tau+1, dim_Logg+1, dim_Teff+1
         READ(40,*) teffinit,logginit,stepteff,steplogg,imax,jmax,kmax
C        PRINT*, teffinit,logginit,stepteff,steplogg,imax, jmax, kmax
C
         ALLOCATE ( tlrhom4(jmax,kmax) , ty(imax,jmax,kmax) , Stat=ios )
         IF ( ios /= 0 ) STOP " TRHO PB 0 4 "
C
         dimlogtau = NINT( ( logtaumax - logtauinit ) / steplogtau ) + 1
         IF ( dimlogtau /= imax ) STOP " TRHO PB 0 5 "
         Write(form(5:6),"(I2.2)") dimlogtau
C
         DO k = 1 , kmax
            teff_calc = teffinit + REAL(k-1) * stepteff
            DO j = 1 , jmax
               READ(40,"(1x,F9.4,4X,F5.0,6X,F3.1)") tlrhom4(j,k), teff_lu, logg_lu
               logg_calc = logginit + REAL(j-1) * steplogg
C
C il y a des tables a la cayrel mal remplie ....
C
               pos = INDEX(nom,"fesh-10a") + INDEX(nom,"fesh-05")
               If ( ( pos == 0 ) .OR. ( teff_calc < 4749.0 ) ) Then          ! ne pas faire le test
                  If ( ABS(logg_calc-logg_lu) > 0.00001 ) Then
                     Write(*,*) logg_calc, logg_lu, teff_calc, teff_lu
                     STOP " TRHO PB 0 6 table à la cayrel mal remplie"
                  End If
                  If ( ABS(teff_calc-teff_lu) > 0.00001 ) Then
                     Write(*,*) logg_calc, logg_lu, teff_calc, teff_lu
                     STOP " TRHO PB 0 7 table à la cayrel mal remplie"
                  End If
               End If
               READ(40,form) (ty(i,j,k), i=1,imax)
            END DO
         END DO
C
         CLOSE(unit=40)
C
c on traite les cas aux limites des tables************
C
         tau_min    = 1.0 / 10.0d0**(-logtauinit)
         teffmax = teffinit+(kmax-1)*stepteff
         loggmax = logginit+(jmax-1)*steplogg
      End If
C
C fin de l'initialisation 
C

c sortie de table en température
	outsideT=teff > teffmax .OR. teff < teffinit
	IF(outsideT .AND. passT)THEN
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001) teff, teffinit, teffmax
	  WRITE(2,1001) teff, teffinit, teffmax
1001	  FORMAT('At least once, call to the T(tau) law with Teff=',es10.3,/,
	1 'outside of [',es10.3,',',es10.3,'], calculation on the limit')	       
	 CASE DEFAULT
	  WRITE(*,1)teff, teffinit, teffmax
	  WRITE(2,1)teff, teffinit, teffmax
1	  FORMAT('Au moins une fois, appel à la loi T(tau) avec Teff=',es10.3,/,
	1 'en dehors de [',es10.3,',',es10.3,'], on se place sur la limite')
	  passT=.FALSE.  	             
	 END SELECT 
	ENDIF

c sortie de table en gravité
	outsideG=logg > loggmax .OR. logg < logginit
	IF(outsideG .AND. passG)THEN
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1002) logg, logginit, loggmax
	  WRITE(2,1002) logg, logginit, loggmax
1002	  FORMAT('At least once, call to the T(tau) law with Log(g)=',es10.3,/,
	1 'outside of [',es10.3,',',es10.3,'], calculation on the limit')	       
 
	 CASE DEFAULT	
	  WRITE(*,2) logg, logginit, loggmax
	  WRITE(2,2) logg, logginit, loggmax
2	  FORMAT('Au moins une fois, appel à la loi T(tau) avec Log(g)=',es10.3,/,
	1 'en dehors de [',es10.3,',',es10.3,'], on se place sur la limite')
	  passG=.FALSE.
	 END SELECT
	ENDIF
	
c on se place à l'intérieur des limites
	teff_calc=MIN(teffmax,MAX(teff,teffinit))
	logg_calc=MIN(loggmax,MAX(logg,logginit))

c détermination du point pivot: code change par B.Pichon avec des "-1"
c  au lieu des "-2"
      logtau = LOG10(tau)
      ip = MIN( imax-1 , MAX ( NINT((logtau-logtauinit)/steplogtau)+1 , 2 ) )
      jp = MIN( jmax-1 , MAX ( NINT((logg_calc-logginit)/steplogg)+1 , 2 ) )
      kp = MIN( kmax-1 , MAX ( NINT((teff_calc-teffinit)/stepteff)+1 , 2 ) )
C comme avant
C      jp = MIN( jmax-2 , MAX ( NINT((logg_calc-logginit)/steplogg)+1 , 2 ) )
C      kp = MIN( kmax-2 , MAX ( NINT((teff_calc-teffinit)/stepteff)+1 , 2 ) )
c
c fin calcul pivot
c
c interpolation sur lrhom4(teff_calc,logg_calc)
c     PRINT*, ' coucou 2', ' jp =',jp,' kp= ',kp
c
      u=(logtau-(logtauinit+steplogtau*(ip-1)))/steplogtau
      v=(logg_calc-(logginit+steplogg*(jp-1)))/steplogg
      w=(teff_calc-(teffinit+stepteff*(kp-1)))/stepteff
c     PRINT*, ' coucou3 ',' u= ',u, ' v= ',v, ' w= ',w
      CALL quadr2(tlrhom4,v,w,jp,kp,lrhom4,yv,yw)
      rhom4=10.d0**lrhom4
      g=10.d0**logg_calc
      dlrhom4dteff=(yw/stepteff)
      drhom4dteff=dlrhom4dteff*2.3026d0*rhom4
      dlrhom4dlogg=(yv/steplogg)
      drhom4dg=dlrhom4dlogg*rhom4/g
c
c interpolation sur t 
c
      CALL quadr3(ty,u,v,w,ip,jp,kp,y,yu,yv,yw,yuu)
c
      t = teff_calc * y
c     PRINT*,' t = ', t
      dtdteff = y + teff_calc * yw / stepteff
c     PRINT*,'dtdteff=',dtdteff
      dtdtau = (teff_calc*yu/steplogtau) * 0.43429d0 / tau 
c     PRINT*,'dtdtau=',dtdtau
      dtdlogg = teff_calc * yv / steplogg
      dtdg = dtdlogg * 0.4343d0 / g
c     PRINT*,' dtdg = ', dtdg
c
c variables supplémentaires f-tau, df_tau, d2f_tau***
c
      ftau=4.d0/3.d0*y**4 
      yprime=yu/steplogtau*0.43429d0/tau
      dftau=16.d0/3.d0*y**3*yprime
      yseconde=yuu*(0.43429d0/steplogtau/tau)**2-yprime/tau
      d2ftau=16.d0*y**2*(y*yseconde/3.d0+yprime**2)
c


      CONTAINS
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE quadr2(tableau,u,v,jp,kp,y,yu,yv)
c
c       interpolation quadratique d'une fonction donnee par tableau
c       au voisinage d'un point pivot predetermine (jp,jk)
c       donne aussi les derivees par rapport aux variables reduites
c       variables reelles divisees par le pas de la variable
c       auteur: R. CAYREL, F95: P. Morel
c
      REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: tableau
      REAL (kind=dp), INTENT(in) :: u,v      
      INTEGER, INTENT(in) :: jp, kp
      REAL (kind=dp), INTENT(out) ::y, yu, yv             
      REAL (kind=dp) :: a1, a2, b11, b22, b12
c
      a1= (tableau(jp+1,kp)-tableau(jp-1,kp))/2.d0
      a2=(tableau(jp,kp+1)-tableau(jp,kp-1))/2.d0
      b11=(tableau(jp+1,kp)+tableau(jp-1,kp)-2.d0*tableau(jp,kp))/2.d0
      b22=(tableau(jp,kp+1)+tableau(jp,kp-1)-2.d0*tableau(jp,kp))/2.d0
      b12=(tableau(jp+1,kp+1)+tableau(jp-1,kp-1)-tableau(jp-1,kp+1)
     1 -tableau(jp+1,kp-1))/4.d0
      y=tableau(jp,kp)+a1*u+a2*v+b11*u**2+b22*v**2+b12*u*v
c     PRINT*,' a1=',a1,' a2=',a2,' b11= ',b11,' b22=',b22,' b12=',b22
      yu= a1+2.d0*b11*u+b12*v
      yv=a2+2.d0*b22*v+b12*u
c
      END SUBROUTINE quadr2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      SUBROUTINE quadr3(tableau,u,v,w,ip,jp,kp,y,yu,yv,yw,yuu)
c
c interpolation quadratique d'une fonction de trois variables donnee par 
c     un tableau. on donne aussi le "pivot" , point du tableau le plus
c     proche du point d'interpolation. le programme calcule aussi les
c     dérivées par rapport aux variables réduites (variables divisées
c     par le pas)
c
c     auteur: R. Cayrel
c     F95: P.Morel
c
      REAL (kind=dp), INTENT(in), DIMENSION(:,:,:) :: tableau
      REAL (kind=dp), INTENT(in) :: u,v,w    
      INTEGER, INTENT(in) :: ip, jp, kp
      REAL (kind=dp), INTENT(out) :: y, yu, yv, yw, yuu
      REAL (kind=dp) :: a1, a2, a3, b11, b22, b33, b12, b23, b13
c
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
c
c     PRINT*,'a1=',a1,'a2=',a2,'a3=',a3
c     PRINT*,'b11=',b11,'b22=',b22,'b33=',b33
c     PRINT*,'b12=',b12,'b13=',b13,'b23=',b23
c
      END SUBROUTINE quadr3
c 
      END SUBROUTINE trho
