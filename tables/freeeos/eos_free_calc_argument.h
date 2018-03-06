C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 2006 Alan W. Irwin
C
C       $Id: eos_free_calc_argument.h 625 2007-07-18 23:49:08Z airwin $
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

C      Common block to hold variables that used to be in the
C      eos_free_calc argument list.  By communicating these variables by
C      this common block, the eos_free_calc argument list is minimalized
C      so that bfgs_iterate (or any other BFGS routine I am testing)
C      can call it easily without carrying along a lot of extra arguments.
      subroutine eos_free_calc(
     &  ifestimate, morder, ifexchange_in,
     &  ifrad,
     &  nux, nuy, nuz, mion_end,
     &  n_partial_aux,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  partial_ions, 
     &  full_sum0, full_sum1, full_sum2, charge, charge2,
     &  dv_pl,
     &  ifdv, ifdvzero, ifcoulomb_mod, if_dc,
     &  ifsame_zero_abundances,
     &  ifexcited,
     &  naux, inv_aux,
     &  inv_ion, max_index,
     &  partial_elements, n_partial_elements, ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass,
     &  ifsame_abundances, ifmtrace, iatomic_number,
     &  ifpi_local,
     &  ifpl, ifmodified, ifh2, ifh2plus,
     &  izlo, izhi, bmin,
     &  nmin, nmin_max, nmin_species, nmax,
     &  eps, tl, bi, h2diss, plop, plopt, plopt2,
     &  r_ion3, r_neutral, nelements, 
     &  ifelement, nion, nions,
     &  nextrasum, maxnextrasum,
     &  match_variable,
     &  fl, rhostar, pstar, sstar, ustar,
     &  f, wf, eta, n_e, rmue,
     &  h, hd, he, rl, h2plus,
     &  sum0, sum2, extrasum, xextrasum,
     &  index_ref)

C      input quantitites:
C         nux, nuy, nuz are number/volume of maximum possible ionization
C         electrons (divided by rho*NA) for hydrogen, helium, and metals.
C         naux is the total number of auxiliary variables.
C         ifexcited > 0 means use excited states (must have Planck-Larkin or
C           ifpi = 3 or 4).
C            0 < ifexcited < 10 means use approximation to explicit summation
C           10 < ifexcited < 20 means use explicit summation
C           mod(ifexcited,10) = 1 means just apply to hydrogen (without molecules)
C             and helium.
C           mod(ifexcited,10) = 2 same as 1 + H2 and H2+.
C           mod(ifexcited,10) = 3 same as 2 + partially ionized metals.
C         ifsame_under = logical variable passed to ionize to control whether
C           to use the same underflow limits or recalculate them.
C           used for removing small discontinuities in eos when attempting
C           to obtain converged solution.
C         inv_ion(nions+2) maps ion index to contiguous ion index used in NR 
C           iteration dot products.
C         max_index is the maximum contiguous index used in NR iteration
C         partial_elements(n_partial_elements+2) index of elements treated as 
C           partially ionized consistent with ifelement.
C         ion_end(nelements) keeps track of largest ion index for each element.
C         ifionized = 0 means all elements treated as partially ionized
C         ifionized = 1 means trace metals fully ionized.
C         ifionized = 2 means all elements treated as fully ionized.
C         if_pteh controls whether pteh approximation used for Coulomb sums (1)
C           or whether use detailed Coulomb sum (0).
C         if_mc controls whether special approximation used for metal part
C           of Coulomb sum for the partially ionized metals.
C         if_mc = 1, use approximation
C         if_mc = 0, don't use approximation
C         ifreducedmass = 1 (use reduced mass in equilibrium constant)
C         ifreducedmass = 0 (use electron mass in equilibrium constant)
C         ifsame_abundances = 1, if call made with same abundances as previous
C           it is the responsibility of the calling programme to set
C           this flag where appropriate.  ifsame_abundances = 0 means
C           slower execution of ionize, but is completely safe against
C           abundance changes.
C         ifmtrace just to keep track of changes in ifelement.
C         iatomic_number(nelements), the atomic number of each element.
C         ifpi = 0, 1, 2, 3, 4 if no, pteh, geff, mhd, or Saumon-style
C           pressure ionization
C         ifpl = 0, 1 if there is no or if there is planck-larkin
C         ifmodified > 0 or not affects ifpi_fit inside excitation_sum.
C         ifh2 non-zero implies h2 is included.
C         ifh2plus non-zero implies h2plus is included.
C         izlo is lowest core charge (1 for H, He, 2 for He+, etc.)
C           for all species.
C         izhi is highest core charge for all species.
C         bmin(izhi) is the minimum bion for all species with a given core charge.
C         nmin(izhi) is the minimum excited principal quantum number for
C           all species with a given core charge.
C         nmin_max(izhi) is the largest minimum excited principal
C           quantum number for all species with a given core charge.
C         nmin_species(nions+2) is the minimum excited principal quantum number
C           organized by species.
C         nmax is maximum principal quantum number included in sum
C           (to be compatible with opal which used nmax = 4 rather than infinity).
C           if(nmax > 300 000) then treated as infinity in qryd_approx.
C           otherwise nmax is meant to be used with qryd_calc only (i.e.,
C           case for mhd approximations not programmed).
C         eps(nelements) is the ratio of mass abundance/atomic mass
C         tl = log(t)
C         tc2 = c2/t
C         bi(nions_in+2) ionization potentials (cm**-1) in ion order:
C           h, he, he+, c, c+, etc.  H2 and H2+ are last two.
C         n.b. the code logic now demands there are iatomic_number ionized
C           states + 1 neutral state for each element.  Thus, the
C           ionization potentials must include all ions, not just
C           first ions (as for trace metals in old code).
C         h2diss, h2 dissociation energy
C         plop(nions+2), plopt, plopt2 planck-larkin occupation probabilities
C           + derivatives.
C         r_ion3(nions+2), *cube of the*
C           effective radii (in ion order but going from neutral
C           to next to bare ion for each species) of MDH interaction between
C           all but bare nucleii species and ionized species.  last 2 are
C           H2 and H2+ (used externally)
C         r_neutral(nelements+2) effective radii for MDH neutral-neutral 
C           interactions.  Last two (used externally) are for for H2 and H2+ (the only ionic
C           species in the MHD model with a non-zero hard-sphere radius).
C         ifelement(nelements) = 1 if element treated as partially ionized
C         (depending on abundance, and trace metal treatment), 0 if
C         element treated as fully ionized or has no abundance.
C         dvzero(nelements) is a zero point shift that is added to dv in ionize,
C           but which is zero for hydrogen.
C         dv(nions) change in equilibrium constant (explained above)
C         dvf(nions) = d dv/d ln f
C         dvt(nions) = d dv/d ln t
C         nion(nions), charge on ion in ion order (must be same order as bi)
C           e.g., for H+, He+, He++, etc.
C         re, ref, ret: fermi-dirac integral and derivatives
C         output quantities:
C         ne, nef, net is nu(e) = n(positive ions)/(Navogadro*rho) and its 
C           derivatives calculated from all elements.
C         sion, sionf, siont is the ideal entropy/R per unit mass and its 
C           derivatives with respect to lnf and ln t.
C         n.b. if ionized, then the contribution from all elements is
C         calculated in ionize.f.  If partially ionized, then ionize separates
C         the hydrogen from the rest of the components.  If molecular, the
C         hydrogen component from ionize is completely ignored, and
C         recalculated in this routine.
C         n.b.  sion returns a component of entropy/R per unit mass,
C         sion = the sum over all non-electron species of
C         -nu_i*[-5/2 - 3/2 ln T + ln(alpha) - ln Na 
C           + ln (n_i/[A_i^{3/2} Q_i]) - dln Q_i/d ln T],
C         where nu_i = n_i/(Na rho), Na is the Avogadro number,
C         A_i is the atomic weight,
C         alpha =  (2 pi k/[Na h^2])^(-3/2) Na
C         (see documentation of alpha^2 in constants.h), and
C         Q_i is the internal ideal partition function
C         of the non-Rydberg states.  (Currently, we calculate this
C         partition function by the statistical weights of the ground states
C         of helium and the combined lower states (roughly
C         approximated) of each of the metals.
C         helium is subsequently corrected for detailed excitation
C         of the non-Rydberg states, and this crude approximation for the metals
C         is not currently corrected.  Thus, in all *current* monatomic
C         cases ln Q_i is a constant, and dln Q_i/d ln T is zero, but this
C         will change for the metals eventually.
C         currently, hydrogen is treated exactly for all cases (full ionization,
C         partial ionization, partial molecular formation).
C         From the equilibrium constant approach and the monatomic species
C         treated in this subroutine (molecules treated outside) we have
C         -nu_i ln (n_i/[A_i^{3/2} Q_i]) =
C         -nu_i * [ln (n_neutral/[A_neutral^{3/2} Q_neutral) +
C         (- chi_i/kT + dv_i)]
C         if we sum this term over all species
C         of an element without molecules we obtain
C         s_element = - eps * ln (eps*n_neutral/sum(n))
C         - eps ln (alpha/(A_neutral^{3/2} Q_neutral)) -
C         - sum over all species of the element of nu_i*(-chi/kT + dv(i))
C         where we have ignored the term
C         eps * (-5/2 - 3/2 ln T + ln rho)
C         (taking into account the first 4 terms above).
C         If hydrogen molecules are included,
C         the result is the same except for the addition of the
C         d ln Q_i/d ln T term (which will also appear for the metals eventually)
C         and the 5/2 factor is multiplied by sum over all species which is
C         corrected to eps by subtracting 5/2 (nu(H2) + nu(H2+)) from sion below.
C         note a final correction of the s zero point occurs
C         in awieos_detailed.f which puts back the ignored terms for both
C         the case of molecules and no molecules.
C         uion is the ideal internal energy (cm^-1 per unit mass divided by
C           avogadro).
C         h, hf, ht is n(H+) and its derivatives with respect to ln f and ln t.
C         hd, hdf, hdt is n(He+) and its derivatives with respect to 
C           ln f and ln t.
C         he, hef, het is n(He++) and its derivatives with respect to 
C           ln f and ln t.
C         the following are returned only if .not.( ifpteh.eq.1.or.if_mc.eq.1)
C         sum0 and f,t,dv derivatives, sum over positive charge number densities
C           with uniform weights.
C         sum2 and f,t,dv derivatives, sum over positive charge number densities
C           weighted by charge^2.
C         the following are returned only if ifpi = 3 or 4.
C         extrasum(maxnextrasum) weighted sums over n(i).
C           for iextrasum = 1,nextrasum-2, sum is only over
C           neutral species + H2+ (the only species with
C           non-zero radii according to the MHD model) and
C           weight is r_neutral^{iextrasum-1}.
C           for iextrasum = nextrasum-1, sum is over all ionized species including
C           bare nucleii, but excluding free electrons, weight is Z^1.5.
C           for iextrasum = nextrasum, sum is over all species excluding bare nucleii
C           and free electrons, the weight is rion^3.
C         extrasumf(maxnextrasum) = partial of extrasum/partial ln f
C         extrasumt(maxnextrasum) = partial of extrasum/partial ln t
C         the following are returned only if ifpl = 1
C         sumpl1 and sumpl2 = weighted sums over non-H nu(i) = n(i)/(rho/H).
C           for sumpl1 sum is over Planck-Larkin occupation
C           probabilities + d ln w/d ln T
C           for sumpl2, sum is over Planck-Larkin d ln w/d ln T
C         sumpl1f and sumpl1t = derivatives of first sum wrt lnf and lnt.
C         rl, rf, rt are the ln mass density and derivatives
C         h2, h2f, h2t are n(H2) and derivatives.
C         h2plus, h2plusf, h2plust are n(H2+) and derivatives.
C         h2plus_dv is the h2plus derivatives wrt dv.  n.b. this vector only
C           returned if molecular hydrogen is calculated.  The calling
C           routine (awieos_detailed) uses it only if 
C           if_mc.eq.1.and.ifh2plus.gt.0.
C         the following are returned only if ifexcited.gt.0 and
C           ifpi.eq.3.or4.
C         xextrasum(4) is the *negative* sum nuvar/(1 + qratio)*
C         partial qratio/partial extrasum(k), for k = 1, 2, 3, and nextrasum-1.
C         xextrasumf(4), xextrasumt(4), xextrasum_dv(nions+2,4) are the
C         fl and tl derivatives (ifnr.eq.0) or dv derivatives (ifnr.eq.1)

C      fixed dimensions for common block variables
      integer*4 nions_efca_dim, naux_efca_dim,
     &  nelements_efca_dim, maxcorecharge_efca_dim,
     &  nstar_efca_dim, nstarsmall_efca_dim,
     &  maxnextrasum_efca_dim, nxextrasum_efca_dim
C      maximum number of ions.
      parameter(nions_efca_dim = 295)
C      maximum number of auxiliary variables (including the possibility of fl)
      parameter(naux_efca_dim = 21)
C      maximum number of elements
      parameter(nelements_efca_dim = 20)
C      maximum core charge for non-bare ion.
      parameter(maxcorecharge_efca_dim = 28)
C      dimensions of rhostar, pstar
      parameter(nstar_efca_dim = 9)
C      dimensions of ustar and sstar
      parameter(nstarsmall_efca_dim = 3)
C      number of MDH auxiliary variables excluding excitation
      parameter(maxnextrasum_efca_dim = 9)
C      number of MDH excitation auxiliary variables
      parameter(nxextrasum_efca_dim = 4)
C      declarations for arguments corresponding to eos_free_calc_in and
C      eos_free_calc_out arguments:
      logical
     &  ifestimate, ifdvzero
      integer*4
     &  morder, ifexchange_in,
     &  ifrad,
     &  mion_end,
     &  n_partial_aux, partial_ions(nions_efca_dim+2),
     &  ifdv(nions_efca_dim+2), ifcoulomb_mod, if_dc,
     &  ifsame_zero_abundances,
     &  ifexcited,
     &  inv_aux(naux_efca_dim),
     &  inv_ion(nions_efca_dim+2), max_index,
     &  partial_elements(nelements_efca_dim+2),
     &  n_partial_elements,
     &  ion_end(nelements_efca_dim+2),
     &  ifionized, if_pteh, if_mc, ifreducedmass,
     &  ifsame_abundances, ifmtrace,
     &  iatomic_number(nelements_efca_dim), ifpi_local,
     &  ifpl, ifmodified, ifh2, ifh2plus,
     &  izlo, izhi,
     &  nmin(maxcorecharge_efca_dim),
     &  nmin_max(maxcorecharge_efca_dim),
     &  nmin_species(nions_efca_dim+2), nmax,
     &  ifelement(nelements_efca_dim), nion(nions_efca_dim+2),
     &  nextrasum, index_ref(nelements_efca_dim)
      real*8
     &  nux, nuy, nuz,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  full_sum0, full_sum1, full_sum2,
     &  charge(nions_efca_dim+2), charge2(nions_efca_dim+2),
     &  plop(nions_efca_dim+2), plopt(nions_efca_dim+2),
     &  plopt2(nions_efca_dim+2),
     &  dv_pl(nions_efca_dim+2),
     &  bmin(maxcorecharge_efca_dim),
     &  eps(nelements_efca_dim), tl, bi(nions_efca_dim+2), h2diss,
     &  r_ion3(nions_efca_dim+2), r_neutral(nelements_efca_dim+2),
     &  match_variable,
     &  fl, rhostar(nstar_efca_dim), pstar(nstar_efca_dim),
     &  sstar(nstarsmall_efca_dim), ustar(nstarsmall_efca_dim),
     &  f, wf, eta, n_e, rmue,
     &  h, hd, he, rl, h2plus,
     &  sum0, sum2,
     &  extrasum(maxnextrasum_efca_dim),
     &  xextrasum(nxextrasum_efca_dim)
C      In order of real*8, integer*4, and logical
      common/block_eos_free_calc_argument/
     &  nux, nuy, nuz,
     &  sum0_mc, sum2_mc, hcon_mc, hecon_mc,
     &  full_sum0, full_sum1, full_sum2,
     &  charge, charge2,
     &  plop, plopt,
     &  plopt2,
     &  dv_pl,
     &  bmin,
     &  eps, tl, bi, h2diss,
     &  r_ion3, r_neutral,
     &  match_variable,
     &  fl, rhostar, pstar,
     &  sstar, ustar,
     &  f, wf, eta, n_e, rmue,
     &  h, hd, he, rl, h2plus,
     &  sum0, sum2,
     &  extrasum,
     &  xextrasum,
     &  morder, ifexchange_in,
     &  ifrad,
     &  mion_end,
     &  n_partial_aux, partial_ions,
     &  ifdv, ifcoulomb_mod, if_dc,
     &  ifsame_zero_abundances,
     &  ifexcited,
     &  inv_aux,
     &  inv_ion, max_index,
     &  partial_elements,
     &  n_partial_elements,
     &  ion_end,
     &  ifionized, if_pteh, if_mc, ifreducedmass,
     &  ifsame_abundances, ifmtrace,
     &  iatomic_number, ifpi_local,
     &  ifpl, ifmodified, ifh2, ifh2plus,
     &  izlo, izhi,
     &  nmin,
     &  nmin_max,
     &  nmin_species, nmax,
     &  ifelement, nion,
     &  nextrasum, index_ref,
     &  ifestimate, ifdvzero
