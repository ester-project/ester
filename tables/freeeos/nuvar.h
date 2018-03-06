C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: nuvar.h 354 2006-04-13 02:37:41Z airwin $
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
C 	common block containing nuvar = n/(rho*navogadro) and related
C       information.
C       nuvar is a function of dv, tl, but dv is a function of fl, and tl.
C       for ifnr = 0:
C       nuvarf is the fl derivative of nuvar(fl, tl)
C       nuvart is the tl derivative of nuvar(fl, tl).
C       for ifnr = 3:
C       nuvarf is the fl derivative of nuvar(dv, tl) = 0
C       nuvart is the tl derivative of nuvar(dv, tl).
C       for ifnr = 1 or ifnr = 3:
C       nuvar_dv is the dv derivative of nuvar(dv, tl).

      integer*4 maxionstage_nuvar, nelements_nuvar
C       must increase this if add element to mix with atomic number greater
C       than 28
      parameter(maxionstage_nuvar = 28)
C       must be the same as nelements in awieos_detailed.f
      parameter(nelements_nuvar = 20)
C       nuvar_nelements is identical to npartial_elements (assigned
C         to that value in ionize.)  Like npartial_elements, it is
C         the number of different elements saved with the compact
C         elements index.

C       nuvar_index_element is identical to partial_elements (assigned
C         to the same values as that array in ionize.)  Like that array
C         it transforms from the compact element index to the element
C          index.

C       nuvar_atomic_number is the atomic number (or the number of neutral
C         and ionized species up to next to bare ion).  It is assigned in
C         ionize to be the same values as iatomic_number although it
C         uses the compact element index rather than the iatomic_number
C         ordinary element index.

C       normally just use npartial_elements, partial_elements, and
C       iatomic_number, but nuvar_nelements, nuvar_index_element, and
C       nuvar_atomic_number assigned in ionize and saved here for
C       the case of external programmes which will not have
C       npartial_elements, partial_elements, and
C       iatomic_number available.

C       N.B. nelements_nuvar index is compact indexed to relevant elements.

      integer*4 nuvar_index_element(nelements_nuvar),
     &  nuvar_atomic_number(nelements_nuvar), nuvar_nelements
C       N.B. first index refers to neutral, first ion, etc. up to next to
C         bare ion (usually).  However, nuvar is also used in eos_sum_calc
C         and eos_free_calc where the bare ion is required (but no derivatives)
C         so leave space for the bare ion value just in the nuvar case.
      real*8
     &  nuvar(maxionstage_nuvar+1,nelements_nuvar),
     &  nuvarf(maxionstage_nuvar,nelements_nuvar),
     &  nuvart(maxionstage_nuvar,nelements_nuvar),
     &  nuvar_dv(maxionstage_nuvar,maxionstage_nuvar,nelements_nuvar)
      common/nuvarblk/nuvar, nuvarf, nuvart, nuvar_dv,
     &  nuvar_index_element, nuvar_atomic_number, nuvar_nelements
