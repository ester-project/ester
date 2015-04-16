
c************************************************************

	SUBROUTINE lit_binaire(chaine,dt)
	
c subroutine public du module mod_exploit
	
c lecture d'un fichier binaire	

c on vérifie l'existence du fichier *_.don puis on le lit, 
c on recherche le fichier *_B.dat, en cas d'échec le fichier *_B.rep
c en cas d'échec on demande d'entrer un autre fichier

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c-----------------------------------------------------------------

	USE mod_donnees, ONLY : alpha, aradia, clight, diffusion,
	1 en_masse, f_eos, f_opa, g, Krot, lim_ro, lsol,
	2 msol, mtot, m_ch, m_qs, m_rot, m_tds, nchim, ne,
	3 nom_abon, nom_atm, nom_conv, nom_ctes, nom_des, nom_diffm,
	4 nom_difft, nom_diffw, nom_elem, nom_etat, nom_fich2, nom_nuc,
	5 nom_nuc_cpl, nom_opa, nom_output, nom_pertm, nom_pertw, nom_tdetau,
	6 nrot, ord_qs, ord_rot, pi, precision, rot_solid, rsol, r_qs, w_rot, zi
	USE mod_kind
	USE mod_variables, ONLY : age, bp, chim, chim_gram, dim_qs, dim_rot,
	1 jlim, knot, knotc, knotr, knot_tds, lconv, mc, mct, model_num,
	2 mrot, mrott, mstar, mw_tot, m_zc, m23, n_ch, n_qs, n_rot, n_tds,
	3 q, qt, rota, rstar, r_ov, r_zc, r2, tds, wrot, x_tds, xt_tds
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(out) :: dt
	CHARACTER (len=255) , INTENT(out):: chaine
	INTEGER :: cas, dim_ch
	
	LOGICAL ::  ok
	
c-----------------------------------------------------------------

2000	FORMAT(8es10.3)

c lecture du fichier *.don
	WRITE(*,20)	
20	FORMAT(/,'Entrer le nom générique du modèle, Ex: soleil')
	READ*,nom_fich2 ; WRITE(*,1)TRIM(nom_fich2)
1	FORMAT('Lecture en binaire des données du modèle: ',a)

c vérification de l'existence du fichier *.don	
	chaine=TRIM(nom_fich2)//'.don'
	INQUIRE(file=TRIM(chaine),exist=ok)
	IF(ok)THEN
	 OPEN(unit=30,form='formatted',status='old',file=TRIM(chaine))
	ELSE
	 PRINT*,'Arrêt: fichier de données inconnu: ',TRIM(chaine)
	 STOP
	ENDIF
	 	
c lecture du fichier de données
	CALL lit_nl(wrot)

c identification du fichier binaire à lire
	PRINT*
	PRINT*,'Pour lire un fichier binaire  *_B.pms, entrer 1' 
	PRINT*,'Pour lire un fichier binaire  *_B.zams, entrer 2' 
	PRINT*,'Pour lire un fichier binaire  *_B.hom, entrer 3' 
	PRINT*,'Pour lire un fichier binaire  *_B.post, entrer 4' 
	PRINT*,'Pour lire un fichier binaire  *_B.cohe, entrer 5' 
	PRINT*,'Pour lire un fichier binaire  *_B.coca, entrer 6' 
	PRINT*,'Pour lire un fichier binaire  *_B.rep, entrer 7' 
	PRINT*,'Pour lire un fichier binaire  *_B.dat, entrer 8' 
	PRINT*,'Pour lire un fichier binaire  ????_B.???, entrer 9'
	READ*,cas
	
	SELECT CASE(cas)
	CASE(1)
	 chaine=TRIM(nom_fich2)//'_B.pms'
	CASE(2)
	 chaine=TRIM(nom_fich2)//'_B.zams'	
	CASE(3)
	 chaine=TRIM(nom_fich2)//'_B.hom'	
	CASE(4)
	 chaine=TRIM(nom_fich2)//'_B.post'		
	CASE(5)
	 chaine=TRIM(nom_fich2)//'_B.cohe'	
	CASE(6)
	 chaine=TRIM(nom_fich2)//'_B.coca'		
	CASE(7)
	 chaine=TRIM(nom_fich2)//'_B.rep'	
	CASE(8)
	 chaine=TRIM(nom_fich2)//'_B.dat'		
	CASE(9)
	 PRINT*,'entrer le nom COMPLET du fichier binaire Ex: bidLe.xV_Q'
	 READ(*,2)chaine
2	 FORMAT(a)	 	
	CASE DEFAULT
	 PRINT*,'Erreur, vous avez saisi : ',cas ; STOP
	END SELECT
		
	INQUIRE(file=TRIM(chaine),exist=ok)
	IF(ok)THEN
	 PRINT*,'On utilise le fichier binaire: ',TRIM(chaine)
	ELSE
	 PRINT*,'Arrêt: fichier binaire inconnu: ',TRIM(chaine) ; STOP
	ENDIF

c lecture des paramètres dans le fichier binaire original			
	OPEN(unit=4,form='unformatted',status='old',file=chaine)       	
	READ(4)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,
	1 n_rot,m_rot,knotr,n_tds,knot_tds	
		
c reprise du modèle, on ajuste les dim. à celles du modèle repris
	 ord_qs=m_qs+r_qs ; dim_qs=knot-ord_qs ; dim_ch=knotc-m_ch
	 ord_rot=m_rot+r_qs ; dim_rot=knotr-ord_rot ; m_tds=knot_tds-n_tds
	 ALLOCATE(bp(ne,dim_qs),q(n_qs),qt(knot),nom_elem(nchim),
	1 chim(nchim,dim_ch),mc(n_ch),mct(knotc),mrot(n_rot),
	2 mrott(knotr),tds(1,knot_tds-m_tds),x_tds(n_tds),
	3 xt_tds(knot_tds),m23(n_qs),r2(n_qs),rota(nrot,dim_rot))    
	 REWIND(unit=4)
	 READ(4)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1 m_rot,knotr,n_tds,knot_tds,mtot,alpha,w_rot,
	2 lim_ro,diffusion,rot_solid,precision,en_masse,f_eos,
	3 f_opa,nom_ctes,nom_pertm,nom_pertw,nom_tdetau,nom_atm,
	4 nom_conv,nom_nuc,nom_nuc_cpl,nom_diffm,nom_difft,nom_diffw,
	5 nom_etat,nom_opa,nom_elem,
	6 bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,m23,r2,m_zc,
	7 r_zc,r_ov,age,dt,mstar,rstar,mw_tot,wrot,jlim,lconv,model_num
        CLOSE(unit=4)
	
	PRINT*,'Fin de lecture du fichier binaire: ',TRIM(chaine)

	RETURN

	END SUBROUTINE lit_binaire	
