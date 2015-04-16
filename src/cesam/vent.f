
c**************************************************************************

	SUBROUTINE vent(xchim,dxchim,jac)
	
c routine public du module mod_nuc

c modifie les dérivées temporelles de la composition chimique
c si le vent (perte/gain de masse) a une composition différente de celle de
c la ZC externe	

c vent est initialisé dans cesam, il n'est appelé par eq_diff_chim ou rk_imps,
c que si l_vent=.TRUE.

c entrées
c	xchim : composition chimique

c entrées/sortie
c	dxchim : dérivée temporelle modifiée
c	jac : jacobien

c	Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ihe4, i_ex, langue, mdot,
	1 modif_chim, nchim, nom_elem, nom_fich2, nucleo, xvent
	USE mod_kind
	USE mod_variables, ONLY : sortie

	IMPLICIT NONE
	
	REAL (kind=dp), OPTIONAL, INTENT(in), DIMENSION(:) :: xchim		
	REAL (kind=dp), OPTIONAL, INTENT(inout), DIMENSION(:,:) :: jac
	REAL (kind=dp), OPTIONAL, INTENT(inout), DIMENSION(:) :: dxchim
	
	REAL (kind=dp), DIMENSION(28) :: vt
	REAL (kind=dp), SAVE ::	flux	
	REAL (kind=dp) :: vt_H, vt_He, vt_Li, vt_Be, vt_B, vt_C,
	1 vt_N, vt_O, vt_F, vt_Ne, vt_Na, vt_Mg, vt_Al, vt_Si,
	2 vt_P, vt_S, vt_Cl, vt_Ar, vt_K, vt_Ca, vt_Sc, vt_Ti,
	3 vt_V,vt_Cr, vt_Mn, vt_Fe, vt_Co, vt_Ni
	
	INTEGER :: i, ili7=0
		
	LOGICAL, SAVE :: init=.TRUE.
	
	CHARACTER (len=50) :: chain
	
	NAMELIST/nl_vent/ vt_H, vt_He, vt_Li, vt_Be, vt_B, vt_C,
	1 vt_N, vt_O, vt_F, vt_Ne, vt_Na, vt_Mg, vt_Al, vt_Si,
	2 vt_P, vt_S, vt_Cl, vt_Ar, vt_K, vt_Ca, vt_Sc, vt_Ti,
	3 vt_V,vt_Cr, vt_Mn, vt_Fe, vt_Co, vt_Ni
		
c----------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 
c il y possibilité de vent si mdot est différent de 0
	 l_vent=mdot /= 0.d0	 
	 IF(.NOT.l_vent)RETURN
	 
c recherche d'une namelist de composition chimique du vent
	 chain=TRIM(nom_fich2)//'.vent'	 
	 INQUIRE(file=TRIM(chain),exist=l_vent)	 
	 IF(.NOT.l_vent)THEN	!l_vent 1
	  chain='vent'	 
	  INQUIRE(file=TRIM(chain),exist=l_vent)
	  IF(.NOT.l_vent)THEN	!l_vent 2
	   SELECT CASE(langue)
	   CASE('english')	
	    WRITE(*,1007) ; WRITE(2,1007)
1007	    FORMAT('CESAM takes account implicitly of the stellar wind')
	   CASE DEFAULT
	    WRITE(*,7) ; WRITE(2,7)
7	    FORMAT('CESAM tient compte implicitement du vent stellaire')
	   END SELECT
	   RETURN
	  ENDIF			!l_vent 2
	 ENDIF			!l_vent 1
	 	  
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

c normalisation par masse  	 
	 vt=vt/SUM(vt)
	  
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1002)TRIM(chain) ; WRITE(2,1002)TRIM(chain)
1002	  FORMAT('normalized abundances by mass from the file : ',a)
	 CASE DEFAULT
	  WRITE(*,2)TRIM(chain) ; WRITE(2,2)TRIM(chain)
2	  FORMAT('abondances normalisées par masse du fichier : ',a)
	 END SELECT
	  
c	 PRINT*,'norm par masse' ; WRITE(*,2000)vt	  
	  
c affectation des abondances par masse des principaux isotopes
c ATTENTION l'algorithme n'est totalement indépendant
c de l'ordre et de la présence des éléments	  
	 ALLOCATE(xvent(nchim)) ; xvent=0.d0

c hydrogène, deutérium et hélium 2 et 4	  
	 IF(nom_elem(2) == ' H2 ')THEN
	  xvent(1)=vt(1)*(1.d0-h2sh1) ; xvent(2)=vt(1)*h2sh1
	  xvent(ihe4-1)=vt(2)*he3she4 ; xvent(ihe4)=vt(2)*(1.d0-he3she4)
	   
	 ELSE
	  xvent(1)=vt(1)
c on suppose H2 dans He3
	  xvent(ihe4-1)=vt(2)*he3she4z ; xvent(ihe4)=vt(2)*(1.d0-he3she4z)	   
	 ENDIF

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
	   xvent(i)=vt(8)*(1.d0-o17so16-o18so16) ; xvent(i+1)=vt(8)*o17so16
	   xvent(i+2)=vt(8)*o18so16
	   i=i+1 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'O18 ')THEN
	   xvent(i)=vt(8)*o18so16 ; CYCLE B1	    
	  ELSEIF(nom_elem(i) == 'Ne20')THEN
	   xvent(i)=vt(10)*(1.d0-ne22sne20) ; xvent(i+1)=vt(10)*ne22sne20
	   i=i+1 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Mg24')THEN	    
	   xvent(i)=vt(12)*(1.d0-mg25smg24*(1.d0+mg26smg25)) ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Mg25')THEN	    
	   xvent(i)=vt(12)*mg25smg24 ; CYCLE B1
	  ELSEIF(nom_elem(i) == 'Mg26')THEN	    
	   xvent(i)=vt(12)*mg25smg24*(1.d0+mg26smg25) ; CYCLE B1    
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
	  WRITE(*,1004) ; WRITE(2,1004)   
1004	  FORMAT('abundances per mole in the wind :')
	 CASE DEFAULT	  
	  WRITE(*,4) ; WRITE(2,4)
4	  FORMAT('abondances par mole du vent :')
	 END SELECT
	 WRITE(*,6)(nom_elem(i),i=1,MIN(nchim,12))
	 WRITE(2,6)(nom_elem(i),i=1,MIN(nchim,12))	  
	 WRITE(*,5)(xvent(i),i=1,MIN(nchim,12))
	 WRITE(2,5)(xvent(i),i=1,MIN(nchim,12))
5	 FORMAT(1x,12es10.3)
	 IF(nchim > 12)THEN	  	
	  WRITE(*,6)(nom_elem(i),i=13,MIN(nchim,24))
	  WRITE(2,6)(nom_elem(i),i=13,MIN(nchim,24))
6	  FORMAT(12(3x,a4,3x))
	  WRITE(*,5)(xvent(i),i=13,MIN(nchim,24))
	  WRITE(2,5)(xvent(i),i=13,MIN(nchim,24))
	 ENDIF
	 
c mesure de sécurité : le fichier vent ne peut être pris en compte que si
c modif_chim = .TRUE.
	 IF(.NOT.modif_chim)THEN
	  SELECT CASE(langue)
	  CASE('english')	
	   WRITE(*,1008) ; WRITE(2,1008)
1008	   FORMAT('STOP, with a file vent, modif_chim = .TRUE. is needed',/,
	1  'erase the file vent, or modify the input file')
	  CASE DEFAULT
	   WRITE(*,8) ; WRITE(2,8)
8	   FORMAT('ARRET, avec un fichier vent il faut modif_chim = .TRUE.',/,
	1  'supprimer le fichier vent, ou modifier le fichier de données')
	  END SELECT
	  CALL sortie	 
	 ENDIF	 
	 
	 RETURN 
	ENDIF		!init
	 
c flux d'élément chimique par unité de masse et de temps dans la ZC externe
	flux=1.d6*mdot/mzc_ext
	DO i=1,nchim
	 dxchim(i)=dxchim(i)+flux*(xvent(i)-xchim(i))
	 jac(i,i)=jac(i,i)-flux
	ENDDO
	
c	PRINT*,'dxchim' ; WRITE(*,2000)dxchim ; PAUSE'vent'
	
	RETURN
	
	END SUBROUTINE vent
		
