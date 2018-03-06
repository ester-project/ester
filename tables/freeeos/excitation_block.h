C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: excitation_block.h 625 2007-07-18 23:49:08Z airwin $
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
C  common block for communications between excitation_pi and excitation_sum
      logical*4 ifhe1_special
      parameter(ifhe1_special = .true.)
      integer*4 nx, max_nmin_max, max_izhi, nions_excitation
C       number of auxiliary variables that qstar depends on
      parameter (nx = 5)
C       maximum value of nmin_max
      parameter(max_nmin_max = 10)
C      maximum value of izhi corresponds to nickel (28) plus room for one
C      extra corresponding to special H2+ calculation.
      parameter(max_izhi = 29)
C      must be greater than or equal to nions
      parameter(nions_excitation = 295)
      real*8 c2t, x(nx), x_old_excitation(nx),
     &  qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,
     &  qstar(max_nmin_max, max_izhi),
     &  qstart(max_nmin_max, max_izhi),
     &  qstarx(nx, max_nmin_max, max_izhi),
     &  qstart2(max_nmin_max, max_izhi),
     &  qstartx(nx, max_nmin_max, max_izhi),
     &  qstarx2(nx, nx, max_nmin_max, max_izhi),
     &  qmhd_he1, qmhd_he1t, qmhd_he1x(nx),
     &  qmhd_he1t2, qmhd_he1tx(nx), qmhd_he1x2(nx,nx),
     &  psum, psumf, psumt, psum_dv(nions_excitation+2),
     &  ssum, ssumf, ssumt, usum,
     &  free_sum, free_sumf, free_sum_dv(nions_excitation+2)
C      Variables for keeping track of conditions for h2/h2+ partition
C      function calculations and qstar_calc calculations for
C      excitation_pi and excitation_sum.
      real*8 tl_old_excitation
      integer ifpi_fit_old_excitation,
     &  ifh2_old_excitation, ifh2plus_old_excitation
      logical ifpl_logical_old_excitation,
     &  ifmhd_logical_old_excitation,
     &  ifapprox_old_excitation, ifdiff_x_excitation
      common /excitation_block/ c2t, x, x_old_excitation,
     &  qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,
     &  qstar, qstart, qstarx, qstart2, qstartx, qstarx2,
     &  qmhd_he1, qmhd_he1t, qmhd_he1x,
     &  qmhd_he1t2, qmhd_he1tx, qmhd_he1x2,
     &  psum, psumf, psumt, psum_dv,
     &  ssum, ssumf, ssumt, usum,
     &  free_sum, free_sumf, free_sum_dv,
     &  tl_old_excitation, ifpi_fit_old_excitation,
     &  ifh2_old_excitation, ifh2plus_old_excitation,
     &  ifpl_logical_old_excitation,
     &  ifmhd_logical_old_excitation,
     &  ifapprox_old_excitation, ifdiff_x_excitation
