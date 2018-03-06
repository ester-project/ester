C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: ionization.h 354 2006-04-13 02:37:41Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@beluga.phys.uvic.ca.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C
C       End of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C*******************************************************************************
      integer*4 nions
      parameter (nions = 295)
C     charge of parent ion in ion order.
      integer*4  nion(nions+2)
      data nion/
C       H
     &  1,
C       He
     &  1,2,
C       C
     &  1,2,3,4,5,6,
C       N
     &  1,2,3,4,5,6,7,
C       O
     &  1,2,3,4,5,6,7,8,
C       Ne
     &  1,2,3,4,5,6,7,8,9,10,
C       Na
     &  1,2,3,4,5,6,7,8,9,10,11,
C       Mg
     &  1,2,3,4,5,6,7,8,9,10,11,12,
C       Al
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,
C       Si
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,
C       P
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
C       S
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
C       Cl
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
C       A
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
C       Ca
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,  
C       Ti
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,            
C       Cr
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,            
C       Mn
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,          
C       Fe
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,    
C       Ni
     &  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,
C       H2, H2+
     &  1,2/              
C       The energy of all relevant species in units of cm**{-1} (the unit
C       used for spectroscopic energy [inverse wavelength] measures, and thus
C       the unit of energy least susceptible to systematic errors due to
C       changes in adopted physical constants.)
C       To convert this energy unit to ergs multiply by hc.
C       Thus, E (in ergs)/kt = c2 E (in cm^{-1})/T where c2 = hc/k = 1.4....
C       Note the energies are organized in the following order: energy of
C       first ion relative to neutral, energy of second ion relative to first
C       ion, energy of third ion relative to second, ...., for each element.
C       These relative energies are what are measured spectroscopically, but
C       for EOS purposes we need the energy relative to the neutral monatomic
C       state (or ultimately the lowest ro-vibrational state of H_2 in the
C       case of the hydrogen element) and these conversions are done later.
      real*8 bi(nions+2)
      data bi/
C       H opal value (R inf) (temporary values)
!     &  109737.31534000d0,
C       He opal values
!     &  197977.09060489d0, 438971.20882307d0,
C       H AEL value
     &  1.09678764d5,
C       He AEL values
     &  1.9831076d5, 4.3890885d5,  
C       C
     &  90820.42d0, 196664.7d0, 386241.0d0, 520178.4d0, 3162395.d0,
     &  3952061.3d0,      
C       N
     &  117225.7d0, 238750.50d0, 382703.8d0, 624866.d0, 789537.2d0,
     &  4452758.d0, 5380089.d0,    
C       O
     &  109837.02d0, 283240.d0, 443084.7d0, 624382.0d0, 918657.d0,
     &  1114010.d0, 5962800.d0, 7028394.d0,   
C       missing F
C       Ne
     &  173929.70d0, 330391.0d0, 511800.d0, 783300.d0, 1018000.d0,
     &  1273800.d0, 1671792.d0, 1928462.d0, 9645005.d0, 10986876.d0,   
C       Na
     &  4.144944d4, 381395.0d0, 577800.d0, 797800.d0, 1116200.d0,
     &  1388500.d0, 1681500.d0, 2130800.d0, 2418700.d0, 11817061.d0,
     &  13297676.d0,              
C       Mg
     &  6.167102d4, 121267.61d0, 646410.d0,
     &  881100.d0, 1139400.d0, 1504300.d0, 1814300.d0, 2144700.d0,
     &  2645200.d0, 2964400.d0, 14210261.d0, 15829951.d0,    
C       Al from NIST
     &  4.827837d4, 
     &  1.518627d5, 229445.7d0, 967804.d0,
     &  1240684.d0, 1536300.d0, 1949900.d0, 2295900.d0, 2662650.d0,
     &  3216100.d0, 3564960.d0, 16268000.d0, 16824529.d0,  
C       Si
     &  6.574776d4, 131838.4d0, 270139.3d0, 364093.1d0, 1345100.d0,
     &  1653900.d0, 1988400.d0, 2445300.d0, 2831900.d0, 3237800.d0,
     &  3839800.d0, 4222400.d0, 19661693.d0, 21560630.d0,  
C       P from NIST
     &  8.458083d4,
     &  1.594515d5, 243600.7d0, 414922.8d0, 524462.9d0, 1777820.d0,
     &  2125800.d0, 2497100.d0, 3001400.d0, 3423200.d0, 3867100.d0,
     &  4523000.d0, 4934000.d0, 22719920.d0, 24759942.5d0,  
C       S from NIST
     &  83559.1d0,
     &  188232.7d0,
     &  280600.d0, 380870.d0, 585514.1d0, 710194.7d0, 2266000.d0,
     &  2651500.d0, 3061300.d0, 3609000.d0, 4071300.d0, 4552500.d0,
     &  5260000.d0, 5702400.d0, 26001523.d0, 28182526.d0,  
C       Cl from Moore
     &  104591.0d0,
     &  192070.d0,
     &  319500.d0, 
     &  431226.d0, 547000.d0, 782600.d0, 921051.d0, 2809100.d0, 
     &  3226700.d0, 3674900.d0, 4268900.d0, 4774700.d0, 5296700.d0,
     &  6047200.d0, 6528300.d0, 29507950.d0, 31829012.d0,  
C       A from Moore
     &  127109.9d0,
     &  222848.2d0,
     &  328600.d0,
     &  482400.d0,
     &  605100.d0, 734040.d0, 1002730.d0, 1157080.d0, 3407300.d0, 
     &  3860900.d0, 4347000.d0, 4986600.d0, 5533800.d0, 6095500.d0,
     &  6894200.d0, 7404400.d0, 33237173.d0, 35699936.d0,  
C       missing K
C       Ca from NIST
     &  49305.95d0,
     &  95751.87d0,
     &  410642.d0,
     &  542600.d0,
     &  681600.d0,
     &  877400.d0, 
     &  1026000.d0, 1187600.d0, 1520640.d0, 1704047.d0, 4774000.d0,
     &  5301000.d0, 5861000.d0, 6595000.d0, 7215000.d0, 7860000.d0,
     &  8770000.d0, 9338000.d0, 41366000.d0, 44117410.d0,  
C       missing Sc
C       Ti from NIST
     &  55010.d0,
     &  109494.d0,
     &  221735.6d0,
     &  348973.3d0,
     &  800900.d0,
     &  964100.d0,
     &  1136000.d0,
     &  1374000.d0, 
     &  1549000.d0, 1741500.d0, 2137900.d0, 2351080.d0, 6354300.d0,
     &  6961000.d0, 7597000.d0, 8420000.d0, 9120000.d0, 9850000.d0,
     &  10860000.d0, 11497000.d0, 50401000.d0, 53440800.d0,  
C       missing V
C       Cr from NIST
     &  54575.6d0,
     &  132966.d0,
     &  249700.d0,
     &  396500.d0,
     &  560200.d0,
     &  731020.d0,
     &  1291900.d0,
     &  1490000.d0,
     &  1688000.d0,
     &  1971000.d0,
     &  2184000.d0, 2404000.d0, 2862000.d0, 3098520.d0, 8151000.d0, 
     &  8850000.d0, 9560000.d0, 10480000.d0, 11260000.d0, 12070000.d0,
     &  13180000.d0, 13882000.d0, 60344000.d0, 63675900.d0,  
C       Mn from NIST
     &  59959.4d0,
     &  126145.0d0,
     &  271550.d0,
     &  413000.d0,
     &  584000.d0,
     &  771100.d0,
     &  961440.d0,
     &  1569000.d0,
     &  1789000.d0,
     &  2003000.d0,
     &  2307000.d0,
     &  2536000.d0, 2771000.d0, 3250000.d0, 3509820.d0, 9152000.d0,
     &  9872000.d0, 10620000.d0, 11590000.d0, 12410000.d0,
C       fake missing Mn XXI I.P.
     &  13000000.d0,
     &  14420000.d0, 15162000.d0, 65660000.d0, 69137400.d0,  
C       Fe
     &  63737.d0, 130563.d0, 247220.d0, 442000.d0, 605000.d0,
     &   799000.d0, 1008000.d0, 1218380.d0, 1884000.d0, 2114000.d0,
     &  2341000.d0, 2668000.d0, 2912000.d0, 3163000.d0, 3686000.d0,
     &   3946280.d0, 10180000.d0, 10985000.d0, 11850000.d0, 12708000.d0,
     &  13620000.d0, 14510000.d0, 15797000.d0, 16500000.d0, 71203000.d0,
     &  74829600.d0,          
C       missing Co
C       Ni from NIST
     &  61600.d0,
     &  146541.56d0,
     &   283800.d0,  443000.d0,  613500.d0,  870000.d0, 1070000.d0,
     &  1310000.d0, 1560000.d0, 1812000.d0, 2589000.d0, 2840000.d0,
     &  3100000.d0, 3470000.d0, 3740000.d0, 4020000.d0, 4606000.d0,
     &   4896200.d0, 12430000.d0, 13290000.d0, 14160000.d0, 15280000.d0,
     &  16220000.d0, 17190000.d0, 18510000.d0,
C       fake missing Ni XXVI
     &  20000000.d0, 
     &  82984000.d0,
     &  86909400.d0,
C       H2 I.P. from Huber and Herzberg, footnote b
     &  1.244172d5,
C       H2+ I.P. to be filled in later.
     &  0.d0/    
      real*8 h2diss
C       h2 dissociation energy cm^{-1} taken from Huber and Herzberg
C       footnote a.
      parameter (h2diss = 3.61183d4)
