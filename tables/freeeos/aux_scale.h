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
C      scaling factors for auxiliary variables to reduce underflow
C      errors.
C      The following parameters must be equal to
C      maxnextrasum, nxextrasum, and iextraoff.
      integer*4 maxnextrasum_scale, nxextrasum_scale,
     &  iextraoff_scale, naux_scale
      parameter (maxnextrasum_scale = 9)
      parameter (nxextrasum_scale = 4)
      parameter (iextraoff_scale = 7)
      parameter (naux_scale = iextraoff_scale +
     &  maxnextrasum_scale + nxextrasum_scale + 1)
C      consistent floating-point underflow limit used for scaled
C      auxiliary variables.
      real*8 aux_underflow
      parameter(aux_underflow = 1.d-250)
C      currently no scaling of the first 5 auxiliary variables (since
C      those variables rarely used and the 4th one is logarithmic
C      in any case).
      real*8 aux_scale(naux_scale),
     &  sum0_scale, sum2_scale,
     &  extrasum_scale(maxnextrasum_scale),
     &  xextrasum_scale(nxextrasum_scale)
      equivalence
     &  (aux_scale(iextraoff_scale-1), sum0_scale),
     &  (aux_scale(iextraoff_scale  ), sum2_scale),
     &  (aux_scale(iextraoff_scale+1), extrasum_scale),
     &  (aux_scale(iextraoff_scale+maxnextrasum_scale+1),
     &  xextrasum_scale)
      common/aux_scale_block/aux_scale
