C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: statistical_weights.h 354 2006-04-13 02:37:41Z airwin $
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
C       statistical weights of ions in ion order from Table IV.
      integer*4 nions_stat
      parameter(nions_stat = 295)
      integer*4 iqion(nions_stat)
      data iqion/
     &  1,              !H
     &  2,1,              !He
     &  6,1,2,1,2,1,            !C
     &  9,6,1,2,1,2,1,          !N
     &  4,9,6,1,2,1,2,1,          !O
     &  6,9,4,9,6,1,2,1,2,1,          !Ne
     &  1,6,9,4,9,6,1,2,1,2,1,        !Na
     &  2,1,6,9,4,9,6,1,2,1,2,1,        !Mg
     &  1,2,1,6,9,4,9,6,1,2,1,2,1,        !Al
     &  6,1,2,1,6,9,4,9,6,1,2,1,2,1,        !Si
     &  9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,      !P
     &  4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,      !S
     &  9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,      !Cl
     &  6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,      !A
     &  2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,    !Ca
     &  28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,  !Ti
     &  6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,  !Cr
     &  7,6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,2,1,!Mn
     &  30,25,6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,2,1,
     &  2,1,  !Fe
     &  10,15,30,25,6,25,28,21,10,1,6,9,4,9,6,1,2,1,6,9,4,9,6,
     &  1,2,1,2,1/  !Ni
C       statistical weights of neutrals in element order from Table IV.
      integer*4 nelements_stat
      parameter(nelements_stat = 20)
      integer*4 iqneutral(nelements_stat)
      data iqneutral/2,1,9,4,9,1,2,1,6,9,4,9,6,1,1,21,7,6,25,21/
