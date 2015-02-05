
c**********************************************************************

      CHARACTER (len=10) FUNCTION arb_rom(i)

c     fonction public du module mod_numerique

c     transforme l'entier 0 <= i <= 30 de notation arabe ==> romaine
c     Exemple: 21 ==> XXI
c     la notation romaine comporte 10 caractères
c     au delà de i=30, la sortie est CCC 

c     P. Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k

c-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(in) :: i

c------------------------------------------------------------------------ 

      SELECT CASE(i)
      CASE(0)
       arb_rom='    0     '
      CASE(1)
       arb_rom='    I     '
      CASE(2)
       arb_rom='    II    '
      CASE(3)
       arb_rom='   III    '
      CASE(4)
       arb_rom='    IV    '
      CASE(5)
       arb_rom='    V     '
      CASE(6)
       arb_rom='    VI    '
      CASE(7)
       arb_rom='   VII    '
      CASE(8)
       arb_rom='   VIII   '
      CASE(9)
       arb_rom='    IX    '
      CASE(10)
       arb_rom='    X     '
      CASE(11)
       arb_rom='    XI    '
      CASE(12)
       arb_rom='   XII    '
      CASE(13)
       arb_rom='   XIII   '
      CASE(14)
       arb_rom='   XIV    '
      CASE(15)
       arb_rom='   XV     '
      CASE(16)
       arb_rom='   XVI    '
      CASE(17)
       arb_rom='   XVII   '
      CASE(18)
       arb_rom='  XVIII   '
      CASE(19)
       arb_rom='  IXX     '
      CASE(20)
       arb_rom='   XX     '
      CASE(21)
       arb_rom='   XXI    '
      CASE(22)
       arb_rom='  XXII    '
      CASE(23)
       arb_rom='  XXIII   '
      CASE(24)
       arb_rom='  XXIV    '
      CASE(25)
       arb_rom='   XXV    '
      CASE(26)
       arb_rom='  XXVI    '
      CASE(27)
       arb_rom='  XXVII   '
      CASE(28)
       arb_rom='  XXVIII  '
      CASE(29)
       arb_rom='  IXXX    '
      CASE(30)
       arb_rom='   XXX    '
      CASE DEFAULT
       arb_rom='   CCC    '
      END SELECT

      RETURN

      END FUNCTION arb_rom
