
c**************************************************************************

	SUBROUTINE vent(xchim,dxchim,jac)
	
c 	routine private du module mod_nuc

c modifie les dérivées temporelles de la composition chimique
c	si le vent a une composition différente de celle de la ZC externe	

c vent n'est appelé que si lvent=.TRUE.,
c 	dans lit_nl lvent est initialisé lvent= mdot /= 0
c	le paramètre de vent p_vent est lu dans la NAMELIST nl_vent

c entrées
c	xchim : composition chimique

c entrées/sortie
c	dxchim : dérivée temporelle modifiée
c	jac : jacobien

c	Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ihe4, i_ex, langue, lvent, mdot,
	1 modif_chim, nchim, nom_elem, nom_fich2, nucleo, p_vent, xvent
	USE mod_kind
	USE mod_variables, ONLY : mstar

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim		
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: jac
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: dxchim
	
	REAL (kind=dp), DIMENSION(28) :: vt
	REAL (kind=dp), SAVE ::	coef	
	REAL (kind=dp) :: vt_H, vt_He, vt_Li, vt_Be, vt_B, vt_C,
	1 vt_N, vt_O, vt_F, vt_Ne, vt_Na, vt_Mg, vt_Al, vt_Si,
	2 vt_P, vt_S, vt_Cl, vt_Ar, vt_K, vt_Ca, vt_Sc, vt_Ti,
	3 vt_V,vt_Cr, vt_Mn, vt_Fe, vt_Co, vt_Ni
	
	INTEGER :: i, ili7=0
		
	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL	:: ok
	
	CHARACTER (len=50) :: chain
	
	NAMELIST/nl_vent/ vt_H, vt_He, vt_Li, vt_Be, vt_B, vt_C,
	1 vt_N, vt_O, vt_F, vt_Ne, vt_Na, vt_Mg, vt_Al, vt_Si,
	2 vt_P, vt_S, vt_Cl, vt_Ar, vt_K, vt_Ca, vt_Sc, vt_Ti,
	3 vt_V,vt_Cr, vt_Mn, vt_Fe, vt_Co, vt_Ni, p_vent
		
c----------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.

c la composition chimique du vent ne peut être différente de celle de la
c couche externe que si modif_chim=.TRUE.

	 IF(.NOT.modif_chim)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1001)modif_chim ; WRITE(2,1001)modif_chim
1001	   FORMAT('modif_chim=',l1,/,
	1  'The wind has the chemical composition of the outer CZ')
	  CASE DEFAULT
	   WRITE(*,1)modif_chim ; WRITE(2,1)modif_chim
1	   FORMAT('modif_chim=',l1,/,
	1  'Le vent a la composition chimique de la ZC externe')
	  END SELECT	 
	  lvent=.FALSE. ; RETURN
	 ENDIF
	 
c recherche de la namelist de la composition chimique du vent
	 chain=TRIM(nom_fich2)//'.vent'	 
	 INQUIRE(file=TRIM(chain),exist=ok)	 
	 IF(.NOT.ok)THEN	!ok 1
	  chain='vent'	 
	  INQUIRE(file=TRIM(chain),exist=ok)
	  IF(.NOT.ok)THEN	!ok 2
	   SELECT CASE(langue)
	   CASE('english')	
	    WRITE(*,1007) ; WRITE(2,1007)
1007	    FORMAT('No file for the abundances by mass in the wind',/,
	1   'CESAM takes the mass ratios of the outer CZ')
	   CASE DEFAULT
	    WRITE(*,7) ; WRITE(2,7)
7	    FORMAT('Pas de rapports de masses pour le vent',/,
	1   'CESAM utilise ceux de la ZC externe')
	   END SELECT
	   lvent=.FALSE. ; RETURN
	  ENDIF			!ok 2
	 ENDIF			!ok 1
	 	  
c lecture de la namelist de la composition chimique du vent	  
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	 WRITE(*,*) ; WRITE(2,*)	
	 READ(3,nl_vent) ; WRITE(*,nl_vent) ; WRITE(2,nl_vent)
	 WRITE(*,*) ; WRITE(2,*)	 	  
 
c normalisation	 
	 vt(1)=vt_H ; vt(2)=vt_He ; vt(3)=vt_Li ; vt(4)=vt_Be
	 vt(5)=vt_B ; vt(6)=vt_C ; vt(7)=vt_N ; vt(8)=vt_O
	 vt(9)=vt_F ; vt(10)=vt_Ne ; vt(11)=vt_Na ; vt(12)=vt_Mg
	 vt(13)=vt_Al ; vt(14)=vt_Si ; vt(15)=vt_P ; vt(16)=vt_S
	 vt(17)=vt_Cl ; vt(18)=vt_Ar ; vt(19)=vt_K ; vt(20)=vt_Ca
	 vt(21)=vt_Sc ; vt(22)=vt_Ti ; vt(23)=vt_V ; vt(24)=vt_Cr
	 vt(25)=vt_Mn ; vt(26)=vt_Fe ; vt(27)=vt_Co ; vt(28)=vt_Ni

	 IF(MAXVAL(ABS((vt))) == 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1003) ; WRITE(2,1003)
1003	   FORMAT('nullity of the mass ratios in the wind'/,
	1   'CESAM takes the chemical composition of the outer CZ')
	  CASE DEFAULT
	   WRITE(*,3) ; WRITE(2,3)
3	   FORMAT('les rapports des abondances du vent sont nuls',/,
	1 'CESAM utilise ceux de la ZC externe')
	  END SELECT
	  lvent=.FALSE. ; RETURN
	 ELSE 

c normalisation par masse  	 
	  vt=vt/SUM(vt)
	  
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1002)TRIM(chain) ; WRITE(2,1002)TRIM(chain)
1002	   FORMAT('normalized abundances by mass from the file : ',a)
	  CASE DEFAULT
	   WRITE(*,2)TRIM(chain) ; WRITE(2,2)TRIM(chain)
2	   FORMAT('abondances normalisées par masse du fichier : ',a)
	  END SELECT
	  
c	  PRINT*,'norm par masse' ; WRITE(*,2000)vt	  
	  
c affectation des abondances par masse des principaux isotopes
c	  ATTENTION l'algorithme n'est totalement indépendant
c	  de l'ordre et de la présence des éléments
	  
	  ALLOCATE(xvent(nchim)) ; xvent=0.d0

c hydrogène et deutérium	  
	  IF(nom_elem(2) == ' H2 ')THEN
	   xvent(1)=vt(1)*(1.d0-h2sh1) ; xvent(2)=vt(1)*h2sh1
	  ELSE
	   xvent(1)=vt(1)
	  ENDIF
	  
c hélium 3 et 4 on suppose H2 dans He4
	  xvent(ihe4-1)=vt(2)*he3she4z ; xvent(ihe4)=vt(2)*(1.d0-he3she4z)

c autres éléments 
	  i=ihe4
	  B1: DO	   
	   i=i+1 ; IF(i > nchim)EXIT B1
	   IF(nom_elem(i) == 'Li7 ')THEN
	    xvent(i)=vt(3) ; ili7=i ; CYCLE B1
	   ELSEIF(nom_elem(i) == 'Li6 ')THEN
	    IF(ili7 /= 0)THEN
	     xvent(i)=vt(3)*li6sli7 ; xvent(ili7)=vt(3)*(1.d0-li6sli7)
	     CYCLE B1
	    ELSE
	     xvent(i)=vt(3) ; CYCLE B1
	    ENDIF  
	   ELSEIF(nom_elem(i) == 'Be9 ')THEN
	    xvent(i)=vt(4) ; CYCLE B1
	   ELSEIF(nom_elem(i) == 'B11 ')THEN
	    xvent(i)=vt(5) ; CYCLE B1
	   ELSEIF(nom_elem(i) == 'C12 ')THEN
	    xvent(i)=vt(6)*(1.d0-c13sc12) ; xvent(i+1)=vt(6)*c13sc12
	    i=i+1 ; CYCLE B1
	   ELSEIF(nom_elem(i) == 'N14 ')THEN
	    xvent(i)=vt(7)*(1.d0-n15sn14) ; xvent(i+1)=vt(7)*n15sn14
	    i=i+1 ; CYCLE B1  
	   ELSEIF(nom_elem(i) == 'O16 ')THEN
	    xvent(i)=vt(8)*(1.d0-o17so16) ; xvent(i+1)=vt(8)*o17so16
	    i=i+1 ; CYCLE B1
	   ELSEIF(nom_elem(i) == 'Ne20')THEN
	    xvent(i)=vt(10) ; CYCLE B1
	   ELSEIF(nom_elem(i) == 'Fe56')THEN
	    xvent(i)=vt(26) ; CYCLE B1	    
	   ENDIF	     
	  ENDDO B1	  
	  
c affectation de l'abondances par masse pour l'extra élément Ex
	  xvent(i_ex)=1.d0-SUM(xvent)
	  
c abondances par mole
	  xvent(1:nchim)=xvent(1:nchim)/nucleo(1:nchim)
	  
c écritures	  
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1004)p_vent ; WRITE(2,1004)p_vent	   
1004	   FORMAT('Mass fraction of the star affected by the wind',es10.3,
	1  /,'abundances by mole in the wind :')
	  CASE DEFAULT	  
	   WRITE(*,4)p_vent ; WRITE(2,4)p_vent
4	   FORMAT('fraction de masse stellaire concernée par le vent',
	1  es10.3,/,'abondances par mole du vent :')
	  END SELECT
	  WRITE(*,6)(nom_elem(i),i=1,MIN(nchim,12))
	  WRITE(2,6)(nom_elem(i),i=1,MIN(nchim,12))	  
	  WRITE(*,5)(xvent(i),i=1,MIN(nchim,12))
	  WRITE(2,5)(xvent(i),i=1,MIN(nchim,12))
5	  FORMAT(1x,12es10.3)
	  IF(nchim > 12)THEN	  	
	   WRITE(*,6)(nom_elem(i),i=13,MIN(nchim,24))
	   WRITE(2,6)(nom_elem(i),i=13,MIN(nchim,24))
6	   FORMAT(12(3x,a4,3x))
	   WRITE(*,5)(xvent(i),i=13,MIN(nchim,24))
	   WRITE(2,5)(xvent(i),i=13,MIN(nchim,24))
	  ENDIF	  
	 ENDIF
	 coef=1.d6/p_vent
	 RETURN	 
	ENDIF		!init

c	PRINT*,'mdot,m_vent' ; PRINT*,'dxchim' ; WRITE(*,2000)dxchim
	 
	IF(lvent)THEN
	 DO i=1,nchim
	  dxchim(i)=dxchim(i)+coef*mdot/mstar*(xvent(i)-xchim(i))
	  jac(i,i)=jac(i,i)-coef*mdot/mstar
	 ENDDO
	ENDIF
	
c	PRINT*,'dxchim' ; WRITE(*,2000)dxchim ; PAUSE'vent'
	
	RETURN
	
	END SUBROUTINE vent
		
