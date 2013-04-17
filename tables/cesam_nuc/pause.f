
c**********************************************************************

	SUBROUTINE pause(mess)
	
c	routine public du module mod_numerique	
   
c	Il suffit de remplacer l'instruction :  PAUSE
c	par : Call PAUSE ()
c	ou : Call PAUSE ( " Pour faire joli : un message qq " ) 
   
c	Auteur: B. Pichon, Departement J.D. Cassini, O.C.A.
c	Modif: P. Morel, Departement J.D. Cassini, O.C.A.

c------------------------------------------------------------------------	
    
	CHARACTER(len=*), OPTIONAL, INTENT(in) :: mess
	
	CHARACTER (len=1) :: suite
	CHARACTER (len=50) :: text
	      
c-----------------------------------------------------------------

	IF (PRESENT(mess))then
	 text='PAUSE: '//ADJUSTL(mess)
	ELSE
	 text='PAUSE'
	ENDIF
	PRINT*,TRIM(text)
		
	PRINT*,'pour arrêter entrer q, poursuivre entrer c'
	PRINT*,'to stop enter q, to continue enter c'
	
	READ*,suite
	IF(suite == 'q')STOP
        
      RETURN
        
	END SUBROUTINE pause
