C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: eff_fit.h 354 2006-04-13 02:37:41Z airwin $
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
C      Common block to communicate with functions that are used
C      to evaluate the basis functions used for pstar, J, and sqrt(I) = K
C      EFF-style fits.
C      nf is the number of f values
C      ng is the number of g values
      integer*4 maxnf, maxng, nf, ng
      parameter (maxnf = 200)
      parameter (maxng = 200)
      real*8 fbasis(maxnf), gbasis(maxng)
C      mf is the number of f basis functions
C      mg is the number of g basis functions
C      mflg and mglg are the equivalent quantities for the f and g basis
C      functions that are multiplied by ln(1+g).  If either of mflg or
C      mglg is zero, then there are no ln(1+g) terms in the approximation.
C      mflg2 and mglg2 are the equivalent quantities for the f and g basis
C      functions that are multiplied by ln^2(1+g).  If either of mflg2 or
C      mglg2 is zero, then there are no ln(1+g) terms in the approximation.
C      mflf and mglf are the equivalent quantities for the f and g basis
C      functions that are multiplied by ln(1+f).  If either of mflf or
C      mglf is zero, then there are no ln(1+f) terms in the approximation.
C      maxmf and maxmg are the maximum values of the various forms of mf and
C      mg.  (For the most economical fits, it turns out that mf and mg are
C      always the largest of the various forms.)
      integer*4 mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf,
     &  maxmf, maxmg
      parameter (maxmf = 10)
      parameter (maxmg = 10)
      common /eff_fit_block/fbasis, gbasis, nf, ng,
     &  mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf
