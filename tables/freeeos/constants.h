C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: constants.h 354 2006-04-13 02:37:41Z airwin $
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
C       physical (and mathematical) constants required by the EOS
      real*8 cr, compton, avogadro, pi, c_e, cd, clight,
     &  electron_mass, h_mass, ct, cpe, c2, alpha20, alpha2,
     &  boltzmann, echarge, ergsperev, stefan_boltzmann, prad_const,
     &  planck, rydberg, bohr, ergspercmm1, proton_mass
C      n.b. all constants are in cgs electrostatic units (sorry about that!)
C      math constant:
      parameter (pi = 3.141592653589793d0)
C      primary constants from 2002 CODATA (the allascii.txt downloaded from
C      http://physics.nist.gov/cuu/Constants/index.html by clicking on the
C      "all values" link. Those are apparently copied from
C      Peter J. Mohr and Barry N. Taylor, CODATA Recommended Values of the
C      Fundamental Physical Constants: 2002, published in Rev. Mod. Phys.
C      vol. 77(1) 1-107 (2005).
      parameter (clight = 2.99792458d10)    !c
      parameter (planck = 6.6260693d-27)    !h
      parameter (cr = 8.314472d7)      !R (gas constant)
      parameter (avogadro = 6.0221415d23)    !Avogadro number N_0
      parameter (rydberg = 109737.31568525d0)  !R(infinity)
      parameter (proton_mass = 1.00727646688d0) ! units are AMU
C      parameter (echarge = 4.8032068d-10)    !Fritz's value in esu's
C      2002 CODATA value transformed to electrostatic units.
      parameter (echarge = 1.60217653d-19*1.d-1*clight)
C      secondary constants (all calculated from primary constants, although
C      I have often specified the rounded (and therefore possibly
C      thermodynamically inconsistent) CODATA value as a commented out
C      value.
C      a Joule is a Coulomb-volt, thus an electron volt is the number of
C      Coulombs in a fundamental charge Joules, or 10^7 times that
C      quantity ergs.
C      parameter (ergsperev = 1.60217653d-19*1.d7)
      parameter (ergsperev = (echarge/(1.d-1*clight))*1.d7)
C      2002 CODATA electron mass in amu.
C      parameter (electron_mass = 5.4857990945d-04)
C      parameter (electron_mass = 9.1093826d-28*avogadro)
      parameter (electron_mass = avogadro*rydberg*
     &    (planck/(echarge*echarge))*
     &    (planck/echarge)*(planck/echarge)*clight/(2.d0*pi*pi))
C      parameter (h_mass = 1.007825035d0)    !amu Wapstra and Audi, 1985
      parameter (h_mass = proton_mass + electron_mass)
C      Boltzmann's constant k calculated consistently with R
      parameter (boltzmann = cr/avogadro)
C      Compton wavelength:
C      parameter (compton = 2.426310238d-10)  !2002 CODATA.
      parameter (compton = planck*avogadro/(electron_mass*clight))
      parameter (c_e = 8.d0*pi/(compton*compton*compton))  !n_e = c_e*rhostar
      parameter (cd = c_e/avogadro)  !8 pi H/lambda_c^3
      parameter (ct = cr/(electron_mass*clight*clight))  !k/m c^2
      parameter (cpe = cr*cd/ct)  !8 pi m c^2/lambda_c^3
C      second radiation constant
C      parameter (c2 = 1.4387752d0)  !2002 CODATA.
      parameter (c2 = planck*clight/boltzmann)
      parameter (alpha20 = c2*c2*cr/(2.d0*pi*clight))
C      alpha2 = (2piHk/h^2)^-3/H^2.
      parameter (alpha2 = alpha20*alpha20*alpha20*
     &  (avogadro/(clight*clight))*(avogadro/clight))
      parameter (stefan_boltzmann = 2.d0*pi*pi*pi*pi*pi/
     &  (15.d0*c2*c2*c2)*clight*boltzmann)
      parameter (prad_const = 4.d0*stefan_boltzmann/(3.d0*clight))
C      parameter (bohr = 0.5291772108d-8)  !Bohr radius 2002 CODATA
      parameter (bohr = 0.5d0*echarge*echarge/(rydberg*clight*planck))
C      number of ergs per cm^-1 = planck*clight = c2*boltzmann
      parameter (ergspercmm1 = c2*boltzmann)
