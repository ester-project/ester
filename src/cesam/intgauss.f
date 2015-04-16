
c*************************************************************************

	SUBROUTINE intgauss(a,b,x,w,n)

c	subroutine public du module mod_numerique

c	initialisation des poids et des abscisses pour l'integration de gauss

c entrees:
c	a, b : intervalle [a,b] d'integration
c	n : nombre de points de Gauss

c sorties:
c	x : abscisses
c	w : poids

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	CESAM95

c----------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: a, b
	INTEGER, INTENT(in) :: n
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: w,x

	INTEGER :: i

c----------------------------------------------------------------------

	SELECT CASE(n)
	CASE(1)
	 x(1)=0.d0 ; w(1)=2.d0

	CASE(2)
	 x(1)= 0.57735 02691 89626d0 ; x(2)=-0.57735 02691 89626d0
	 w(1)= 1.d0 ; w(2)=1.d0

	CASE(3)
	 x(1)= 0.d0 ; x(2)= 0.77459 66692 41483d0
	 x(3)=-0.77459 66692 41483d0
	 w(1)= 0.88888 88888 88889d0 ; w(2)= 0.55555 55555 55556d0
	 w(3)= 0.55555 55555 55556d0

	CASE(4)
	 x(1)= 0.33998 10435 84856d0 ; x(2)= 0.86113 63115 94053d0
	 x(3)=-0.33998 10435 84856d0 ; x(4)=-0.86113 63115 94053d0
	 w(1)= 0.65214 51548 62546d0 ; w(2)= 0.34785 48451 37454d0
	 w(3)= 0.65214 51548 62546d0 ; w(4)= 0.34785 48451 37454d0
	 
	CASE(5)
	 x(1)= 0.d0 ; x(2)= 0.53846 93101 05683d0
	 x(3)= 0.90617 98459 38664d0 ; x(4)=-0.53846 93101 05683d0
	 x(5)=-0.90617 98459 38664d0
	 w(1)= 0.56888 88888 88889d0 ; w(2)= 0.47862 86704 99366d0
	 w(3)= 0.23692 68850 56189d0 ; w(4)= 0.47862 86704 99366d0
	 w(5)= 0.23692 68850 56189d0
	 
	CASE(6)
	 x(1)= 0.23861 91860 82197d0 ; x(2)= 0.66120 93864 66265d0
	 x(3)= 0.93246 95142 03152d0 ; x(4)=-0.23861 91860 82197d0
	 x(5)=-0.66120 93864 66265d0 ; x(6)=-0.93246 95142 03152d0
	 w(1)= 0.46791 39345 72691d0 ; w(2)= 0.36076 15730 48139d0
	 w(3)= 0.17132 44923 79190d0 ; w(4)= 0.46791 39345 72691d0
	 w(5)= 0.36076 15730 48139d0 ; w(6)= 0.17132 44923 79190d0
	
	CASE DEFAULT
	 PRINT*,'Erreur sur le nombre de points dans intgauss, n=',n
	 PRINT*,'n doit etre compris entre 1 et 6, bornes incluses'
	 PRINT*,'ARRET' ; STOP
	END SELECT
		
	DO i=1,n
	 x(i)=(b-a)/2.d0*x(i)+(b+a)/2.d0 ; w(i)=w(i)/2.d0*(b-a)
	ENDDO

	RETURN

	END SUBROUTINE intgauss
