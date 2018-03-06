C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 2006 Alan W. Irwin
C
C       $Id: bfgs.h 412 2007-01-25 02:52:03Z airwin $
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

C      Common block variables and definitions used for bfgs algorithm
C      using notation from
C      R. Fletcher, "Practical Methods of Optimization" 2nd ed.
C      John Wiley and Sons, 1987.
C      N.B. in following comments this reference is referred to as PMOO.

C      maximum dimensionality of solution vector, gradient, etc.
      integer*4 nmax_bfgs
      parameter(nmax_bfgs=300)
C      fletcher_estimate is .true. when fletcher_estimate of alpha1 is
C      chosen.
      logical fletcher_estimate
C      alpha_previous is used to store the alpha solution of the previous
C      line search. This quantity is only used when the GSL estimate of
C      alpha1 is chosen.
C      deltaf_previous is used to store the delta f of the previous line
C      search.  This quantity is only used when the Fletcher estimate of
C      alpha1 is chosen.
      real*8 alpha_previous, deltaf_previous
      common/block_bfgs/alpha_previous, deltaf_previous,
     &  fletcher_estimate
