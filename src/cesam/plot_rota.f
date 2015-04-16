
c****************************************************************

	SUBROUTINE plot_rota

c routine subordonnée de ecrit_rota
c dessins des variables de la rotation s'il n'y a pas de dessin en ligne

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c-----------------------------------------------------------------

	USE mod_donnees, ONLY : nchim
	USE mod_numerique, ONLY : min_max
	USE mod_variables, ONLY : age
	
	IMPLICIT NONE		

	REAL (kind=sp), ALLOCATABLE, DIMENSION(:,:) :: ychim	
	REAL (kind=sp), ALLOCATABLE, DIMENSION(:) :: deff, chi, chi_mu,
	1 chi_t, dh, dlldlm, dlnro_nu, dlnro_t, domega, du, dv, eps_mu, eps_nuc,
	2 eps_T,  f_eps, flux, grad_mu, grand_k, lambda, m, mu, n2t, n2mu,
	2 omega, psi, r, ro, t, theta, u, xchim
	REAL (kind=sp) :: ymax, ymin, xmax, xmin		

	INTEGER :: ic=1

	CHARACTER (len=80) :: chaine, device, number, text
	CHARACTER (len=100) :: titre

c allocations
	ALLOCATE(chi(n_rot),chi_mu(n_rot),chi_t(n_rot),deff(n_rot),
	1 dlldlm(n_rot),dh(n_rot),dlnro_nu(n_rot),dlnro_t(n_rot),domega(n_rot),
	2 du(n_rot),dv(n_rot),eps_mu(n_rot),eps_nuc(n_rot),
	3 eps_T(n_rot),flux(n_rot),f_eps(n_rot),grad_mu(n_rot),grand_k(n_rot),
	4 lambda(n_rot),m(n_rot),mu(n_rot),n2mu(n_rot),n2t(n_rot),omega(n_rot),
	5 psi(n_rot),r(n_rot),ro(n_rot),t(n_rot),theta(n_rot),u(n_rot),
	6 xchim(n_rot),ychim(nchim,n_rot))
	
c les variables
	r(:)=ecrit(1,:)
	m(:)=ecrit(2,:)
	omega(:)=ecrit(3,:)/ecrit(3,n_rot)
	u(:)=ecrit(4,:)
	theta(:)=ecrit(5,:)
	lambda(:)=ecrit(6,:)
	psi(:)=ecrit(7,:)
	flux(:)=ecrit(8,:)
	deff(:)=ecrit(9,:)
	dh(:)=ecrit(10,:)
	dv(:)=ecrit(11,:)
	t(:)=ecrit(12,:)
	ro(:)=ecrit(13,:)
	grad_mu(:)=ecrit(14,:)
	f_eps(:)=ABS(ecrit(15,:))
	eps_T(:)=ecrit(16,:)
	eps_nuc(:)=ecrit(17,:)
	dlldlm(:)=ABS(ecrit(18,:))
	eps_mu(:)=ecrit(19,:)
	chi(:)=ecrit(20,:)
	chi_t(:)=ecrit(21,:)
	chi_mu(:)=ecrit(22,:)
	grand_k(:)=ecrit(23,:)
	dlnro_nu(:)=ecrit(24,:)
	dlnro_t(:)=ecrit(25,:)
	domega(:)=ecrit(26,:)
	du(:)=ecrit(27,:)
	n2t(:)=ecrit(28,:)
	n2mu(:)=ecrit(29,:)
	mu(:)=ecrit(30,:)
	ychim(1:nchim,:)=ecrit(ncoeff+1:nb_var,:)

c on passe en LOG10 en évitant les 0	
	WHERE(f_eps > 1.d-9)
	 f_eps=LOG10(f_eps)
	ELSEWHERE
	 f_eps=-9.e0
	END WHERE
	WHERE(eps_nuc > 1.d-9)
	 eps_nuc=LOG10(eps_nuc)
	ELSEWHERE
	 eps_nuc=-9.e0
	END WHERE
	WHERE(dlldlm > 1.d-9)
	 dlldlm=LOG10(dlldlm)
	ELSEWHERE
	 dlldlm=-9.e0
	END WHERE

c LOG10 de Deff, Dh, Dv en évitant les valeurs nulles
	deff=MAX(deff,0.001) ; dh=MAX(dh,0.001) ; dv=MAX(dv,0.001)
	deff=LOG10(deff) ; dh=LOG10(dh) ; dv=LOG10(dv)

c dessins en R / R*
	xmin=-0.05 ; xmax=1.05 ; r=r/MAXVAL(r) ; text='R/R\d\(0856)\u'
c	xmin=-0.05 ; xmax=MAXVAL(r)*1.05 ; text='R/R\d\(2281)\u'  !R/Rsol
	
	device=TRIM(nom_fich2)//'-rot.ps/ps'

c	PRINT*,'device pour rotation ? /xw, /PS, /CPS'
c	READ*,device
	CALL pgbegin(0,device,1,1)

c pour le titre	 
	 WRITE(number,10)model_num
10	 FORMAT(i4.4)
	 chaine=TRIM(nom_fich2)//', modele : '//number
	
c dessin de omega
	CALL pgsci(1)	 
	titre='Profil de \gW\ds\u du modèle '//TRIM(chaine)	
	CALL min_max(omega,n_rot,ymax,ymin) 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'\gW\gW\ds\u',titre)	 
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,omega)
		
c dessin de domega/dnu
	CALL pgsci(1)	 
	titre='Profil de d\gW/d\gn du modèle '//TRIM(chaine)	
	CALL min_max(domega,n_rot,ymax,ymin) 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'d\gW/d\gn',titre)	 
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,domega)

c dessin de U	 
	CALL pgsci(1)
	titre='Profil de U du modèle '//TRIM(chaine)
	ymin=-10. ; ymax=-2.
	u=ABS(U) ; WHERE(u == 0.)u=1.d-8 ; u=LOG10(u)
	CALL min_max(u,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'Log\d10\u(U cm s\u-1)',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,u)
	
c dessin de dU/dnu	 
	CALL pgsci(1)	 
	titre='Profil de dU/d\gn du modèle '//TRIM(chaine)
	CALL min_max(du,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'dU/d\gn ',titre)	
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,du)	

c dessin de Theta/Psi	 
	CALL pgsci(1)
	CALL min_max(theta,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	IF(Krot ==3)THEN		 
	 titre='Profil de \gH du modèle '//TRIM(chaine)
	 CALL pglabel(text,'\gH',titre)	 
	ELSE		 
	 titre='Profil de \gQ du modèle '//TRIM(chaine)
	 CALL pglabel(text,'\gQ',titre)	 
	ENDIF	
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,theta)
	 
c dessin de lambda	 
	CALL pgsci(1)	 
	titre='Profil de \gL du modèle '//TRIM(chaine)
	CALL min_max(lambda,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'\gL',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,lambda)
	 
c dessin de Psi/Phi	 
	CALL pgsci(1)	 
	titre='Profil de \gQ du modèle '//TRIM(chaine)
	CALL min_max(psi,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	IF(Krot ==3)THEN		 
	 titre='Profil de \gQ du modèle '//TRIM(chaine)
	 CALL pglabel(text,'\gQ',titre)	 
	ELSE		 
	 titre='Profil de \gF du modèle '//TRIM(chaine)
	 CALL pglabel(text,'\gF cm\u2',titre)	 
	ENDIF	
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,psi)
	
c dessin de flux	 
	CALL pgsci(1)	 
	titre='Profil du flux de moment cinétique du modèle '//TRIM(chaine)
	CALL min_max(flux,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'flux de moment cinétique',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,flux)
	
c dessin de Deff	 
	CALL pgsci(1)
	titre='Profil de LOG10 Deff du modèle '//TRIM(chaine)
	CALL min_max(deff,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'LOG10 Deff',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,deff)

c dessin de Dh	 
	CALL pgsci(1)
	titre='Profil de LOG10 Dh du modèle '//TRIM(chaine)
	CALL min_max(dh,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'LOG10 Dh',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,dh)

c dessin de Dv	 
	CALL pgsci(1)
	titre='Profil de LOG10 Dv du modèle '//TRIM(chaine)
	CALL min_max(dv,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'LOG10 Dv',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,dv)
	
c dessin de f_epsilon	 
	CALL pgsci(1)
	titre='Profil de LOG10 f\d\ge\u du modèle '//TRIM(chaine)
	CALL min_max(f_eps,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'LOG10 f\d\ge',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,f_eps)

c dessin de mu
	CALL pgsci(1)
	titre='Profil de \gm du modèle '//TRIM(chaine)
	CALL min_max(mu,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'\gm ',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ; CALL pgline(n_rot,r,mu)

c dessin de grad_mu	 
	CALL pgsci(1)
	titre='Profil de \(0583)\d\gm\u du modèle '//TRIM(chaine)
	CALL min_max(grad_mu,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'\(0583)\d\gm',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,grad_mu)
	
c dessin de f_epsilon
	CALL pgsci(1)
	titre='Profil de LOG10 f\d\ge\u du modèle '//TRIM(chaine)
	CALL min_max(f_eps,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'LOG10 f\d\ge',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,f_eps)	
	
c dessin de n2t	 
	CALL pgsci(1)
	titre='Profil de N2t du modèle '//TRIM(chaine)
	CALL min_max(n2t,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'N2t',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,n2t)
	
c dessin de n2mu	 
	CALL pgsci(1)
	titre='Profil de N2\gm du modèle '//TRIM(chaine)
	CALL min_max(n2mu,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'N2\gm',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,n2mu)

c dessin de ro	 
	CALL pgsci(1)
	titre='Profil de \gr du modèle '//TRIM(chaine)
	CALL min_max(ro,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'\gr',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,ro)
	
c dessin de d ln ro / d nu	 
	CALL pgsci(1)
	titre='Profil de dln\gr/d\gn du modèle '//TRIM(chaine)
	CALL min_max(dlnro_nu,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'dln\gr/d\gn',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,dlnro_nu)
	
c dessin de d ln ro / d t	 
	CALL pgsci(1)
	titre='Profil de dln\gr/dt du modèle '//TRIM(chaine)
	CALL min_max(dlnro_t,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'dln\gr/dt',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,dlnro_t)	
	
c dessin de eps_nuc	 
	CALL pgsci(1)
	titre='Profil de LOG10 \ge\dnuc\u du modèle '//TRIM(chaine)
	CALL min_max(eps_nuc,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'LOG10 \ge\dnuc',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,eps_nuc)

c dessin de eps_T	 
	CALL pgsci(1)
	titre='Profil de \ge\dT\u du modèle '//TRIM(chaine)
	CALL min_max(eps_T,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'\ge\dT',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,eps_T)	

c dessin de eps_mu	 
	CALL pgsci(1)
	titre='Profil de \ge\d\gm\u du modèle '//TRIM(chaine)
	CALL min_max(eps_mu,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'\ge\d\gm',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,eps_mu)	
	
c dessin de dlldlm	 
	CALL pgsci(1)
	titre='Profil de LOG10(\ge\dg\u+\ge\dnuc\u)/\ge_bar du modèle '//TRIM(chaine)
	CALL min_max(dlldlm,n_rot,ymax,ymin)	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'LOG10(\ge\dg\u+\ge\dnuc)/\ge_bar',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,dlldlm)

c dessin de chi
	CALL pgsci(1)
	titre='Profil de \gx du modèle '//TRIM(chaine)
	CALL min_max(chi,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'\gx',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,chi)

c dessin de chi_T
	CALL pgsci(1)
	titre='Profil de \gx\dT\u du modèle '//TRIM(chaine)
	CALL min_max(chi_t,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'\gx\dT',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,chi_t)

c dessin de chi_mu
	CALL pgsci(1)
	titre='Profil de \gx\d\gm\u du modèle '//TRIM(chaine)
	CALL min_max(chi_mu,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'\gx\d\gm',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,chi_mu)

c dessin de grand_K
	CALL pgsci(1)
	titre='Profil de K du modèle '//TRIM(chaine)
	CALL min_max(grand_k,n_rot,ymax,ymin)
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)
	CALL pglabel(text,'K',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,grand_k)
	
c dessin de X	 
	CALL pgsci(1)
	titre='Profil de X du modèle '//TRIM(chaine)
	xchim(:)=ychim(1,:)
	ymax=1.05 ; ymin=-0.05	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'X',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,xchim)
	
c dessin de Y	 
	CALL pgsci(1)
	titre='Profil de Y du modèle '//TRIM(chaine)
	xchim(:)=ychim(ihe4,:)
	ymax=1.05 ; ymin=-0.05	 
	CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	CALL pglabel(text,'Y',titre)
	ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,xchim)
	
c dessin de li7
	IF(ili7 > 1)THEN
	 CALL pgsci(1)
	 titre='Profil de Li7 du modèle '//TRIM(chaine)
	 xchim(:)=ychim(ili7,:)
	 CALL min_max(xchim,n_rot,ymax,ymin)	 
	 CALL pgenv(xmin,xmax,ymin,ymax,0,0)	
	 CALL pglabel(text,'Li7',titre)
	 ic=MAX(2,mod(ic+1,9)) ; CALL pgsci(ic) ;  CALL pgline(n_rot,r,xchim)
	ENDIF
	
	DEALLOCATE(deff,dh,dv,flux,grad_mu,lambda,m,omega,psi,r,ro,
	1 t,theta,u,xchim,ychim)		

	CALL pgend

	RETURN
	
	END SUBROUTINE plot_rota
