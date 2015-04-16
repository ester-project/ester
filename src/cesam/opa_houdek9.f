
c*************************************************************************

      SUBROUTINE opa_houdek9(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)

c     routine private du module mod_opacite

c     opacités OPAL 95
c     interpolation from Guenter Houdek: bi-rational splines version9
c     This routine calls the file ~/HOUDEK/v9/OPINTPATH_AX
c     which gives the path to the opacity tables

c     The routine also calls three routines:
c     maceps (optimization for the compiler)
c     opinit
c     opints  for the interpolation
c
c     The compilation needs 1 more libraries: opint,
c     linker avec :

c sur DUPARC  /users/user2/morel/v9/lib/libopint.a 
c     tabnam='/users/user2/v9/OPINTPATH_AX'

c ATTENTION
c     en cas de sortie de table on appelle opa_opal2_co
c     le nom et le chemin de la table d'opacité OPINTPATH_AX doit
c         NECESSAIREMENT être
c     indiqué dans f_opa(2) de la NAMELIST nk_opa, du fichier de données
c     mon_modele.don
c     f_opa(2)='/users/user2/morel/v9/OPINTPATH_AX'
c     le nom et le chemin de la table d'opacité opa_opal2_co_etendu.dat
c     doit être indique dans f_opa(1)
c     f_opa(1) ='/users/user2/morel/STAR_DATA/opa_opal2_co_etendu.dat'

c     il y a d'autres possibilités, en particulier supprimer
c     les opacités conductives, voir la notice de Guenter
 
c     CESAM95

c entrées :
c     xchim(1)=X : comp. chim. en H
c     t : température K
c     ro : densité cgs

c sorties :
c     kapa : opacité gr / cm2)
c     dkapdt : kappa / d t
c     dkapdr : kappa / d densite              
c     dkapdx : kappa / d xchim(1)

c     Z est obtenu par 1-X-Y, Y=He4+He3 ou encore Z=Z0

c-------------------------------------------------------------------------

      USE mod_donnees, ONLY : f_opa, ihe4, nchim, nom_chemin, z0
      USE mod_kind

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(:) :: xh
      REAL (kind=dp), INTENT(in) :: t, ro
      REAL (kind=dp), INTENT(out) ::  kappa, dkapdt, dkapdr, dkapdx

      REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: Lxh
      REAL (kind=dp) :: tlg, rlg, opalg, opr, opt, opx, opz, x, z,
     1 eps, Lt, Lro

      INTEGER, PARAMETER :: iord=4, imode=3
      INTEGER :: iexp, ier

      LOGICAL, SAVE :: init=.TRUE., hxp=.TRUE., hrp=.TRUE., hzp=.TRUE.,
     1 htp=.TRUE., hrm=.TRUE., htm=.TRUE.
      LOGICAL :: hou, ok

      CHARACTER (len=80) :: tabnam

c-----------------------------------------------------------------

2000  FORMAT(8es10.3)

      IF(init)THEN         ! reads the tables at the first call
       init=.FALSE.
   
	WRITE(*,1)TRIM(f_opa(2)) ; WRITE(2,1)TRIM(f_opa(2))
1	FORMAT('opacités à Z variable, tables de OPAL 95 :',/,
	1 'interpolation par splines birationnelles de G. Houdek',/,
	2 'version 9 adaptation à CESAM2k de P.Morel, données de ',a)

c lecture et initialisation des tables
       tabnam=TRIM(nom_chemin)//TRIM(f_opa(2))
       INQUIRE(file=TRIM(tabnam),exist=ok)
       IF(.NOT.ok)THEN
        WRITE(*,2)TRIM(f_opa(2)) ; WRITE(2,2)TRIM(f_opa(2))
2	FORMAT('ARRET, le nom de la table d''opacité Houdek9 est inconnu : ',a)
	STOP
       ENDIF

c      PRINT*,tabnam ; CALL pause('tabnam')

       CALL maceps(eps)   !precision machine
       CALL opinit(eps,iord,tabnam,imode) !lecture des tables

       x=0.7d0 ; z=.02d0 ; tlg=4.d0 ; rlg=0.d0
       CALL opints(x,z,tlg,rlg,opalg,opr,opt,opx,opz,iexp,ier)
       PRINT* ; PRINT*,'fin des initialisations des opacités Houdek 95'
       PRINT*
       ALLOCATE(Lxh(nchim))
      ENDIF
 
c transformation of the inputs of OPA  to those of OPINTF:
      x=abs(xh(1))
      IF(nchim > 1)THEN
       z=abs(1.d0-xh(1)-xh(ihe4)-xh(ihe4-1))
      ELSE
       z=z0
      ENDIF
      
c variables de OPAL en LOG10
      tlg=LOG10(t) ; rlg=LOG10(ro/(t*1.d-06)**3) ; hou=.TRUE.      

c encadrement des inputs
      IF(x > 0.8d0)THEN
       IF(x > 0.81d0)THEN
        IF(hxp)THEN
         WRITE(2,11) ; WRITE(*,11)
11       FORMAT(/,'opa_houdek9, on rencontre au moins une fois:',/,
     1   'X > 0.81, appel à opa_opal2_co')
         hxp=.FALSE.
        ENDIF
        hou=.FALSE.
       ENDIF
       x=0.8d0
      ENDIF

      IF(z > 0.1d0)THEN
       IF(hzp)THEN
        WRITE(2,12) ; WRITE(*,12)
12      FORMAT(/,'opa_houdek9, on rencontre au moins une fois:',/,
     1  'Z > 0.10, appel à opa_opal2_co')
        hzp=.FALSE.
       ENDIF
       hou=.FALSE.
      ENDIF

      IF(tlg < 3.d0)THEN
       IF(tlg < 2.9d0)THEN
        IF(htm)THEN
         WRITE(2,13) ; WRITE(*,13)    
13       FORMAT(/,'opa_houdek9, on rencontre au moins une fois:',/,
     1   'T < 1000, appel à opa_opal2_co')
         htm=.FALSE.
        ENDIF
        hou=.FALSE.
       ENDIF
       tlg=3.d0
      ELSEIF(tlg > 8.7d0)THEN
       IF(tlg > 8.8d0)THEN
        IF(htp)THEN
         WRITE(2,14) ; WRITE(*,14)    
14       FORMAT(/,'opa_houdek9, on rencontre au moins une fois:',/,
     1   'T > 6.4d8, appel à opa_opal2_co')
         htp=.FALSE.
        ENDIF
        hou=.FALSE.
       ENDIF
       tlg=8.7d0
      ENDIF

      IF(rlg < -8.d0)THEN
       IF(rlg < -8.05d0)THEN
        IF(hrp)THEN
         WRITE(2,15) ; WRITE(*,15)    
15       FORMAT(/,'opa_houdek9, on rencontre au moins une fois:',/,
     1   'rlg <  -8, appel à opa_opal2_co')
         hrp=.FALSE.
        ENDIF
        hou=.FALSE.
       ENDIF
       rlg=-8.d0
      ELSEIF(rlg > 1.d0)THEN
       IF(rlg > 1.5d0)THEN
        IF(hrm)THEN
         WRITE(2,16) ; WRITE(*,16)
16       FORMAT(/,'opa_houdek9, on rencontre au moins une fois:',/,
     1   'rlg > 1, appel à opa_opal2_co')
         hrm=.FALSE.
        ENDIF
        hou=.FALSE.
       ENDIF
       rlg=1.d0
      ENDIF
 
c CALL OF THE SUBROUTINE OF G. Houdek:
      IF(hou)THEN
       CALL opints(x,z,tlg,rlg,opalg,opr,opt,opx,opz,iexp,ier)
       IF(ier == 1)THEN
        WRITE(*,*)'opa_houdek9, sortie de table et extrapolation'
       ENDIF

       kappa =10.d0**opalg            ! opacity
       dkapdt =(kappa /t )*opt         ! d kappa / d T
       dkapdr =(kappa /ro )*opr        ! d kappa / d ro
       IF(x > 0.d0)THEN
        dkapdx =(kappa /x )*opx         ! d kappa / d X
       ELSE
        dkapdx =0.d0
       ENDIF

c      WRITE(*,*)'kappa,dkapdr,dkapdt,dkapdx'
c      WRITE(*,2000)kappa,dkapdr,dkapdt,dkapdx
c      CALL pause5('les kappa'

      ELSE
       Lxh(1:nchim)=xh(1:nchim) ; Lt=t ; Lro=ro
       CALL opa_opal2(Lxh,Lt,Lro,kappa,dkapdt,dkapdr,dkapdx,.FALSE.)   
      ENDIF

      RETURN

      END SUBROUTINE opa_houdek9
