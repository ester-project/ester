
c*************************************************************************

      SUBROUTINE diffm_mp(p,t,r,m,ro,drox,w,gradrad,dgradradx,xi,d,dd,v,dv)

c     routine private du module mod_evol

c     calcul des coefficients de diffusion microscopique
c     d'après Michaud-Proffit, Vienne, Inside the stars IAU 137, p.250
c     on tient compte de Z: Vy=-xVx/(1-x-z)

c     dérivée de diff_z16, diff_der, diff_mipr

c     les coefficients de la matrice de diffusion ne sont non nuls que sur la
c     diagonale et la première colonne.

c     Dimensions dans le programme appelant
c     d(nchim,nchim), dd(nchim,nchim,nchim), v(nchim), dv(nchim,nchim)

c     convention de notation :
c     équation de diffusion dXi/dt=dFi/dm + nuclear, i=1,nchim
c     Fi=4 pi r**2 ro (4 pi r**2 ro D.dX/dm - Vi Xi)

c     d=D=(di1, di2,... din) avec Dij coefficient de d Xj / dm
c     dans le produit scalaire D.dX/dm=sum d_ij d Xj / dm

c     pour ligne d'indice i
c     v(i) coefficient de x_i,
c     dv(i,k)=dérivée v_i / x_k
c     seule la premiere colonne de dv
c     est non nulle (pas de dérivées / Xi, i /= 1)
c     d(i,j)=coefficient d_ij de d x_j / dm
c     dd(i,j,k)= dérivée de d_ij / x_k

c     Auteur: P.Morel
c     Conseils: J. Matias, DASGAL, G. Alécian, Evry
c     adaptation à CESAM: P. Morel, Département J.D. Cassini, O.C.A.

c     CESAM2k

c    19 08 96 : remplacement ln lambda par Cij (c.a.d. cs)
c    11 12 98 : cohérence avec le cas général, modification de v, dv, d, dd
c    18 12 98 : abondances par mole en entrée
c    05 01 99 : yVy=-xVx au lieu de (1-x-z)Vy=xVx

c entrees
c     p, t, r, l, m, ro: données au point de calcul
c     xi: composition chimique, par mole
c     kap: opacité
c     gradad, gradrad: gradients

c sorties
c     d0, dd : coefficients d_ij de d x_j / d m et dérivées / x_k
c     v0, dv : coefficients v_i de x_i et dérivées / x_k

c     ces quantités doivent etre initialisées à 0 dans le programme appellant

c----------------------------------------------------------------

      USE mod_donnees, ONLY : ah, amu, echarg, g, ihe4, kbol, langue, msol,
     1 nchim, nucleo, pi, rsol, zi
      USE mod_kind

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(:) :: xi
      REAL (kind=dp), INTENT(in) :: gradrad, dgradradx, drox, m, p, r,
     1 ro, t, w
      REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd   
      REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d, dv
      REAL (kind=dp), INTENT(inout), DIMENSION(:) :: v    
      
      REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: as, sa
      REAL (kind=dp), DIMENSION(nchim,2) :: cs, csr, csx, lnlamb,
     1   lnlambr, lnlambx
      REAL (kind=dp), DIMENSION(nchim) :: xchim
      REAL (kind=dp), PARAMETER :: thetae=1.d0
      REAL (kind=dp), SAVE :: b, cte1 
      REAL (kind=dp) :: c, cx, c1, c1x, c2, c2x, den, denx, dxsyx, dxsyy,
     1   grav, gravx, num, numx, q, qx, vth, vthx, xsy

      INTEGER :: i, k
      
      LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------------

2000  FORMAT(8es10.3)
2001  FORMAT(10es8.1)

      IF(init)THEN
       init=.FALSE.
       SELECT CASE(langue)
       CASE('english')
        WRITE(*,1010) ; WRITE(2,1010)
1010     FORMAT(/,'Diff. micro. of Michaud & Proffit, in RZ',/,
     1  'H1 et He4 are diffused, the other species are trace',/,
     2  'for He4: X Vx = - Y Vy ==> Vy = - X / (1 - X - Z) Vx',/)       
       CASE DEFAULT
        WRITE(*,10) ; WRITE(2,10)
10      FORMAT(/,'Diff. micro. de Michaud & Proffit, dans ZR, der. anal. ',/,
     1  'H1 et He4 sont diffusés, les autres sont élément trace',/,
     2  'pour He4: X Vx = - Y Vy ==> Vy = - X / (1 - X - Z) Vx',/)
       END SELECT
       
c      Case mixture of H et He & thetae=1.   :   zi=1  zj=2

       cte1=g*msol/rsol**2
       b=15.d0/16.d0*SQRT(2.d0*amu/5.d0/pi)*SQRT(kbol)**5/echarg**4
           
c      masses reduites

       ALLOCATE(as(nchim,2), sa(nchim,2))
       DO i=1,nchim
        as(i,1)=nucleo(i)*nucleo(1)/(nucleo(i)+nucleo(1))
        sa(i,1)=SQRT(as(i,1))
        as(i,2)=nucleo(i)*nucleo(ihe4)/(nucleo(i)+nucleo(ihe4))
        sa(i,2)=SQRT(as(i,2))
       ENDDO
       
c      convention de signe:
c      équation de diffusion dX/dt=dF/dm + nuclear
c      F=4 pi r**2 ro (4 pi r**2 ro D.dX/dm - v X)
c      le terme en gravité vient de dlnP/dr, il entre dans v (eq.17)
c      et vth (eq.19).
c      - Pour v, il y a déjà 2 signes - (eq 13 et 17)
c      mettre - v revient bien a mettre grav < 0
c      - Pour vth, comme grav > 0, et qu'il y a - dans v (eq.17)
c      il faut alors mettre - vth pour rétablir le signe
c      en - v dans eq_diff_chim
           
      ENDIF       !initialisation
      
c     les équations et les dérivées analytiques en xi

c     les formules de Michaud-Proffit utilisent les abondances par gramme
c     on transforme les xi par mole en xchim par gramme
c     à l'issue du calcul les dérivées sont prises / xi par mole

      xchim=xi*nucleo ; xsy=xchim(1)/xchim(ihe4) ; dxsyx=xsy/xchim(1)
      dxsyy=-xsy/xchim(ihe4)

c     logarithme de Coulomb et intégrale cij pour H=X et He4=Y

      DO i=1,nchim    
       CALL coulomb(zi(i),zi(1),thetae,ro,xchim(1),
     1  t,lnlamb(i,1),lnlambr(i,1),lnlambx(i,1),
     2  cs(i,1),csr(i,1),csx(i,1))
       lnlambx(i,1)=lnlambx(i,1)+lnlambr(i,1)*drox    !dérivée/X
       csx(i,1)=csx(i,1)+csr(i,1)*drox
         
       CALL coulomb(zi(i),zi(ihe4),thetae,ro,xchim(1),
     1 t,lnlamb(i,2),lnlambr(i,2),lnlambx(i,2),
     2 cs(i,2),csr(i,2),csx(i,2))
       lnlambx(i,2)=lnlambx(i,2)+lnlambr(i,2)*drox    !dérivée/X
       csx(i,2)=csx(i,2)+csr(i,2)*drox
      ENDDO
       
c     lnlambdaij est remplacé par cij, remarque sous la formule 17     
       
      c=b*SQRT(t)**5/ro/cs(1,2)/(0.7d0+0.3d0*xchim(1))
      cx=-c*(0.3d0/(0.7d0+0.3d0*xchim(1))+
     1 drox/ro+csx(1,2)/cs(1,2))   !dérivée /X
      q=2.d0/SQRT(5.d0)/ro*b*SQRT(t)**5
      qx=-q*drox/ro               !dérivée /X
      IF(r*m > 0.d0)THEN
       grav=cte1*ro*m/p/r**2+2.d0/3.d0*r*w**2
      ELSE
       v=0.d0 ; dv=0.d0 ; d=0.d0 ; dd=0.d0 ; RETURN
      ENDIF      
      gravx=grav*drox/ro  !dérivée/X

      vth=0.54d0*b/ro*(4.75d0*xchim(1)+2.25d0)*SQRT(t)**5/(5.d0+cs(1,2))*
     1   gradrad*grav        !formule 19
      vthx=0.54d0*b/ro*4.75d0*SQRT(t)**5/(5.d0+cs(1,2))*gradrad*grav
     1 +vth*(-drox/ro+dgradradx/gradrad
     2 +gravx/grav-csx(1,2)/(5.d0+cs(1,2)))

      DO i=1,nchim
       IF(i == 1)THEN !pour l'hydrogène 
        v(1)=c*(1.25d0+1.125d0*gradrad)*grav*(1.d0-xchim(1))
        dv(1,1)=v(1)*(cx/c+1.125d0*dgradradx/(1.25d0+1.125d0*gradrad)
     1  +gravx/grav-1.d0/(1.d0-xchim(1)))  !dérivée /X
        IF(c < 0.d0 .OR. c > 0.d0)THEN     
         d(1,1)=c*(3.d0+xchim(1))/(1.d0+xchim(1))/(3.d0+5.d0*xchim(1))
         dd(1,1,1)=d(1,1)*(cx/c+1.d0/(3.d0+xchim(1))-1.d0/(1.d0+xchim(1))
     1   -5.d0/(3.d0+5.d0*xchim(1)))
        ELSE
         d(1,1)=0.d0 ; dd(1,1,1)=0.d0
        ENDIF
                         
       ELSEIF(i == ihe4)THEN  !pour He4: XVx=-Y Vy ==> Vy=-X/(1-X-Z)Vx
        d(ihe4,ihe4)=d(1,1)*xsy   !coeff de d Y
        dd(ihe4,ihe4,1)=dd(1,1,1)*xsy+d(1,1)*dxsyx    !der/X1
        dd(ihe4,ihe4,ihe4)=d(1,1)*dxsyy   !der/Y
        v(ihe4)=-v(1)*xsy
        dv(ihe4,1)=-dv(1,1)*xsy-v(1)*dxsyx            !der/X1
        dv(ihe4,ihe4)=-v(1)*dxsyy         !der/Y
      
       ELSE   !pour les autres elements: elements test
        num=sa(i,1)*cs(i,1)-sa(i,2)*cs(i,2)
        numx=sa(i,1)*csx(i,1)-sa(i,2)*csx(i,2)
        den=xchim(1)*num+sa(i,2)*cs(i,2)
        denx=num+xchim(1)*numx+sa(i,2)*csx(i,2)
        d(i,i)=q/zi(i)**2/den
        dd(i,i,1)=d(i,i)*(qx/q-denx/den)  !dérivée d(i,i)/X1
           
        c1=1.d0/(1.d0+xchim(1))   -10.d0/(5.d0*xchim(1)+3.d0)
        c1x=-1./(1.d0+xchim(1))**2+50.d0/(5.d0*xchim(1)+3.d0)**2
           
        c2=d(i,i)*c1+xchim(1)*(num/den-0.23d0)*d(1,1)
        c2x=dd(i,i,1)*c1+d(i,i)*c1x+
     1  (num/den-0.23d0)*d(1,1)+xchim(1)*(numx-num*denx/den)/den*d(1,1)
     2  +xchim(1)*(num/den-0.23d0)*dd(1,1,1)
           
        d(i,1)=c2*ah/nucleo(i)*xchim(i) ; dd(i,1,1)=c2x*ah/nucleo(i)*xchim(i)
        dd(i,1,i)=c2*ah/nucleo(i)      
      
        c1=1.d0+zi(i)-nucleo(i)*(5.d0*xchim(1)+3.d0)/4.d0
        c1x=-nucleo(i)*5.d0/4.d0 ; c2=xchim(1)*(num/den-0.23d0)*v(1)*ah
        c2x=c2*(1.d0/xchim(1)+dv(1,1)/v(1)
     1   +(numx-num*denx/den)/den/(num/den-0.23d0))
        v(i)=d(i,i)*c1*grav+c2-vth
        dv(i,1)=dd(i,i,1)*c1*grav+d(i,i)*c1x*grav+d(i,i)*c1*gravx+c2x-vthx   
       ENDIF
      ENDDO
 
c     dérivées par rapport aux abondances par gramme

      DO k=1,nchim
       dv(:,k)=dv(:,k)*nucleo(k) ; dd(:,:,k)=dd(:,:,k)*nucleo(k)
      ENDDO
      
      RETURN 

      END SUBROUTINE diffm_mp
