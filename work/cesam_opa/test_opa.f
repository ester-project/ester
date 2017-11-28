       program test_opa
! Test les opacity de CESAM2k

        USE mod_opa
        USE mod_donnees, ONLY : langue, f_opa, nchim, nom_opa,x0,y0,z0,ihe4,NOM_CHEMIN

        IMPLICIT NONE

       real*8 ro,t,xchim(3)
       real*8 kap,dkapt,dkapro,dkapx
       
       langue='PSE'
       NOM_CHEMIN='/home/rieutord/Ester/work/cesam_opa/'
       f_opa(1)='opa_yveline.bin'
       nom_opa='opa_yveline'
       nchim=3
       !x0=0.70d0
       x0=0.706454d0
       x0=7.0645714230d-1
       z0=0.02d0
       !y0=1d0-x0-z0
       xchim(1)=x0
       xchim(2)=1d0-x0-z0; ihe4=2
       xchim(3)=0.02d0
       !xchim(3)=8.827E-05
       xchim(4)=3.425E-03
       xchim(5)=4.128E-05
       xchim(6)=1.059E-03
       xchim(7)=4.167E-06
       xchim(8)=9.640E-03
       xchim(9)=3.903E-06
       xchim(10)=5.827E-03

!7.000E-01, He3 : 8.827E-05, He4 : 2.799E-01, C12 : 3.425E-03, C13 : 4.128E-05
!N14 : 1.059E-03, N15 : 4.167E-06, O16 : 9.640E-03, O17 : 3.903E-06, Ex : 5.827E-03


       t=2.8d+03
       ro=1.38d-07

       CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)

       print*,'opacity = ',kap

       end
