
C******************************************************************

	SUBROUTINE etat_gong1(p,t,xchim,deriv,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c	routine public du module mod_etat

c	�quation d'�tat pour GONG �tape 1

c	modifs:
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb

c	Auteur: P.Morel, D�partement J.D. Cassini, O.C.A.
c	CESAM2k

c entree :
c	p : pression
c	t : temp�rature
c	xchim : composition chimique
c	deriv=.false. : �vite le calcul de certaines d�riv�es

c sortie :
c	ro : densit� et d�riv�es
c	u : �nergie interne et d�riv�es

c----------------------------------------------------------------

	use mod_donnees, only : z0, aradia, ah, ahe4, granr
	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim	
	REAL (kind=dp), INTENT(in) :: p, t
	LOGICAL, INTENT(in) :: deriv
		
	REAL (kind=dp), INTENT(out) :: ro, drop, drot, drox,
	1 delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
	2 gradad, dgradadp, dgradadt, dgradadx, alfa, beta, gamma1,
	3 u, dup, dut, dux
     
	REAL (kind=dp), save :: cte1, aradias3

	REAL (kind=dp) :: x, mum1, dmum1x, y, zisai, zai, unsai,
	1 cpp, drott, drotp, drotx, dutp, dutt, dutx	

	LOGICAL, SAVE ::  init=.TRUE.

c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	WRITE(*,*)'entr�e etat_gong1 p,t,xchim(1),ah,ahe4,zai'
c	WRITE(*,2000)p,t,xchim(1),ah,ahe4,zai

	IF(init)THEN	!initialisations
	 init=.FALSE.	!Z proportion en elements lourds doit avoir
			!ete initialise par appel fictif a opa
	 cte1=3.d0/2.d0*granr
	 aradias3=aradia/3.d0
	 
	 zisai=0.5d0	![ Zi / Ai ]h		page 9

	 unsai=0.0625d0	![ 1 / Ai ]h
	 zai=(zisai+unsai)*z0	![ Zi + 1 / Ai ]h

	 WRITE(2,1)
	 WRITE(*,1)
1	 FORMAT(/,' EOS GONG1: tot. ioni., no Prad, no d�g�n�.',//)
	ENDIF

	x=xchim(1)
	y=1.d0-x-z0

	mum1=2.d0*x/ah+3.d0*y/ahe4+zai		!2.2

	ro=p/granr/mum1/t	!2.1
	drop= ro/p
	drot=-ro/t
	u=cte1*mum1*t		!�nergie interne par gramme 2.3
	dup=0.d0
	dut=u/t
	delta=-t/ro*drot
	cpp=p/ro/t*delta
	cp=dut+cpp	
	gradad=p/ro*delta/cp/t		!gradient adiabatique
	alfa=p/ro*drop
	beta=1.d0-aradias3*t**4/p
	gamma1=1.d0/(alfa-delta*gradad)

	IF(deriv)THEN	
	 dmum1x=2.d0/ah-3.d0/ahe4
	 drox=-ro*dmum1x/mum1
	 dux=u*dmum1x/mum1
	 drott=-2.d0*drot/t
	 drotp=-drop/t
	 drotx=-drox/t

	 dutp=0.d0
	 dutt=0.d0
	 dutx=dux/t

	 deltap=delta*(-drop/ro+drotp/drot)
	 deltat=delta*(-drot/ro+drott/drot+1.d0/t)
	 deltax=delta*(-drox/ro+drotx/drot)

	 dcpp=dutp+cpp*(-drop/ro+deltap/delta+1.d0/p)
	 dcpt=dutt+cpp*(-drot/ro+deltat/delta-1.d0/t)
	 dcpx=dutx+cpp*(-drox/ro+deltax/delta)

	 dgradadp=gradad*(-drop/ro+deltap/delta-dcpp/cp+1.d0/p)
	 dgradadt=gradad*(-drot/ro+deltat/delta-dcpt/cp-1.d0/t)
	 dgradadx=gradad*(-drox/ro+deltax/delta-dcpx/cp)
	ENDIF

c	WRITE(*,*)'sortie gong1 ro,drop,drot,drox,drott,drotp,drotx/u'
c	WRITE(*,2000)ro,drop,drot,drox,drott,drotp,drotx
c	WRITE(*,2000)u,dup,dut,dux,dutt,dutp,dutx

	END SUBROUTINE etat_gong1
