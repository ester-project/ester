C      Coulomb fitting factors used to fit various detailed EOS results.
C      These are the log10(gamma) limits of the DH region and the modified
C      OCP region.
      real*8 xdh10, xmocp10
      data xdh10, xmocp10/-0.4d0, 0.d0/
C      Pressure-ionization fitting factors used to fit various detailed
C      EOS results.
      integer*4 nelements_pi_fit
      parameter (nelements_pi_fit = 20)
      real*8 pi_fit_neutral_ln_original(nelements_pi_fit+2),
     &  pi_fit_neutral_ln(nelements_pi_fit+2),
     &  pi_fit_neutral_ln_saumon(nelements_pi_fit+2)
C       adjust these values for neutral-neutral pressure ionization
C       best fit to MDH results
      data pi_fit_neutral_ln_original/0.d0, 0.d0, 18*0.d0, 2*0.d0/
C       best fit to opal results
!  data pi_fit_neutral_ln/0.5d0, 0.5d0, 18*0.5d0, -50.d0/
C       attempt to fit Saumon H results
!  data pi_fit_neutral_ln/-0.050d0, 0.5d0, 18*0.5d0, -2.d0/
C       keep neutral metals same as original mdh for simplicity
!      data pi_fit_neutral_ln/-0.050d0, 0.5d0, 18*0.d0, -2.d0/
C       try neutrals identical to Saumon fit for opal fit, metals hydrogenic
C       like excited-state pi_fit factors
C      Force best fit, quad=0 results for Saumon
!      data pi_fit_neutral_ln/0.8d0, 0.d0, 18*0.d0, -2.7d0/
C      Same for quad=10 Saumon fit.
!      data pi_fit_neutral_ln/-0.4d0, -0.2d0, 18*0.d0, 2*-4.0d0/
C      Same for quad=20 Saumon fit.
! best pre eos2001 fit.
!      data pi_fit_neutral_ln/0.d0, -0.2d0, 18*0.d0, -3.6d0/
! eos2001 fit
!      data pi_fit_neutral_ln/-3.d0, -4.d0, 18*-5.d0, -5.d0, -5.d0/
! eos2001 fit (refined solar, Helium fit to start)
!      data pi_fit_neutral_ln/-0.4d0, -0.2d0, 18*0.d0, 2*-4.0d0/
!      data pi_fit_neutral_ln/-0.4d0, -4.d0, 18*0.d0, 2*-4.0d0/
C      Go back to overall eos2001 fit since it gives reasonable fit
C      for solar case.
!      data pi_fit_neutral_ln/-3.d0, -4.d0, 18*-5.d0, -5.d0, -5.d0/
C      Revert back to 1.1.0
      data pi_fit_neutral_ln/-0.4d0, -0.2d0, 18*0.d0, -4.0d0, -40.d0/
C       fit to Saumon H, and He tables, metals don't matter
C      Best fit, quad=0 results.
!      data pi_fit_neutral_ln_saumon/0.8d0, 0.d0, 18*0.d0, -2.7d0/
C      Best fit, quad = 10 results.  (Note effective neutral radii
C      have to be reduced from their quad=0 results.)
!      data pi_fit_neutral_ln_saumon/-0.4d0, -0.2d0, 18*0.d0, -4.0d0/
C      Best fit, quad = 20 results.  This comparison with Saumon went
C      to higher density where we had to change coefficients such
C      that H+ reduced and H2 reduced.
!      data pi_fit_neutral_ln_saumon/0.d0, -0.2d0, 18*0.d0, 2*-3.6d0/
C      Revert back to 1.1.0.
      data pi_fit_neutral_ln_saumon/-0.4d0, -0.2d0, 18*0.d0,
     &  -4.0d0, -40.d0/
      integer*4 nions_pi_fit
      parameter (nions_pi_fit = 295)
      real*8 pi_fit_ion_ln_original(nions_pi_fit+2),
     &  pi_fit_ion_ln(nions_pi_fit+2),
     &  pi_fit_ion_ln_saumon(nions_pi_fit+2)
C       best fit to MDH results
      data pi_fit_ion_ln_original/0.d0, 0.d0, 0.d0, 292*0.d0, 0.d0,
     &  0.d0/
C       best fit to opal results (H2 and H2+ constants not relevant)
!  data pi_fit_ion_ln/-1.7d0, 0.d0, 0.2d0, 292*-1.d0, -50.d0, -50d0/
C       attempt to get some rough agreement with other pressure ionization
C       modes and also Saumon high density H results.
!  data pi_fit_ion_ln/-2.13d0, 0.d0, 0.2d0, 292*-1.d0, -0.7d0, -1.d0/
C       try opal fit as close to Saumon fit as possible
C       (i.e. reduced H, He, low ions of metals, H2, H2+  ion interactions).
C       further modification to high metal ions to force full metal ionization
C       at high rho, T
C       further modification to get monotonic metals for LMS model, change
C       -1.8 to zero for first two ions
      data pi_fit_ion_ln/
!Attempt to get smoothest LMS fit when unconstrained by eos2001
!     &  -2.5d0, !H
!OLD Attempt to solar fit eos2001
!     &  -40.d0, !H
!20041211 Attempt to solar fit eos2001 (start with quad=10 result)
!     &  -1.8d0, !H
C      back to overall eos2001 fit/
!     &  -40.d0, !H
!     &  -2.5d0, -2.5d0, !He
!Attempt to get smoothest LMS fit when unconstrained by eos2001
!     &  -1.5d0, -2.5d0, !He
!Old Attempt to solar fit eos2001
!     &  -40.d0, -5.d0, !He
!20041211 Attempt to solar fit eos2001 (start with quad=10 result)
!     &  -1.8d0, -1.8d0, !He
!     &  -40.d0, -5.d0, !He
!     &  0.d0, 0.d0,  4*6.d0, !C
!     &  -1.5d0, -2.5d0,  4*3.d0, !C
!     &  -1.5d0, -2.5d0,  5*3.d0, !N
!     &  -1.5d0, -2.5d0,  6*3.d0, !O
!     &  -1.5d0, -2.5d0,  8*3.d0, !Ne
!     &  -1.5d0, -2.5d0,  9*3.d0, !Na
!     &  -1.5d0, -2.5d0, 10*3.d0, !Mg
!     &  -1.5d0, -2.5d0, 11*3.d0, !Al
!     &  -1.5d0, -2.5d0, 12*3.d0, !Si
!     &  -1.5d0, -2.5d0, 13*3.d0, !P
!     &  -1.5d0, -2.5d0, 14*3.d0, !S
!     &  -1.5d0, -2.5d0, 15*3.d0, !Cl
!     &  -1.5d0, -2.5d0, 16*3.d0, !A
!     &  -1.5d0, -2.5d0, 18*3.d0, !Ca
!     &  -1.5d0, -2.5d0, 20*3.d0, !Ti
!     &  -1.5d0, -2.5d0, 22*3.d0, !Cr
!     &  -1.5d0, -2.5d0, 23*3.d0, !Mn
!     &  -1.5d0, -2.5d0, 24*3.d0, !Fe
!     &  -1.5d0, -2.5d0, 26*3.d0, !Ni
!     &  2*-1.8d0/ !H2, H2+
! quad=20 result used by Cassisi student for thesis work
!     &  -4.5d0, -4.5d0/ !H2, H2+
! quad=20 result adjusted for monotonic H charge fraction (H2+ + H+),
! suppression of H2+ above 10^5 K, and avoidance of grada sign change
! all for 0.075 model locus.
!     &  -2.0d0, -2.5d0/ !H2, H2+
! and suppress H2 and H2+ some more now I have discontinuity index to worry
! about
! Best fit that minimizes discontinuity for pre eos2001
!     &  -4.0d0, -4.5d0/ !H2, H2+
! eos_2001 fit.  Start with Saumon
!     &  -40.d0, -40.d0/ !H2, H2+
C      Revert back to 1.1.0
     &  -1.8d0, !H
     &  -1.8d0, -1.8d0, !He
     &  0.d0, 0.d0,  4*6.d0, !C
     &  0.d0, 0.d0,  5*6.d0, !N
     &  0.d0, 0.d0,  6*6.d0, !O
     &  0.d0, 0.d0,  8*6.d0, !Ne
     &  0.d0, 0.d0,  9*6.d0, !Na
     &  0.d0, 0.d0, 10*6.d0, !Mg
     &  0.d0, 0.d0, 11*6.d0, !Al
     &  0.d0, 0.d0, 12*6.d0, !Si
     &  0.d0, 0.d0, 13*6.d0, !P
     &  0.d0, 0.d0, 14*6.d0, !S
     &  0.d0, 0.d0, 15*6.d0, !Cl
     &  0.d0, 0.d0, 16*6.d0, !A
     &  0.d0, 0.d0, 18*6.d0, !Ca
     &  0.d0, 0.d0, 20*6.d0, !Ti
     &  0.d0, 0.d0, 22*6.d0, !Cr
     &  0.d0, 0.d0, 23*6.d0, !Mn
     &  0.d0, 0.d0, 24*6.d0, !Fe
     &  0.d0, 0.d0, 26*6.d0, !Ni
     &  -1.8d0, -1.8d0/ !H2, H2+
C       best fit to saumon H and He table has negligible ion interactions
      data pi_fit_ion_ln_saumon/-40.d0, -40.d0, -40.d0, 294*-40.d0/

C      quad pressure-ionization data.
      real*8 quad_original, quad_modified
      parameter(quad_original = 10.d0)
C      all modified forms now_use_this term (adjusted for best
C      possible fit and convergence I could find for highest density Saumon H table.)
!      parameter(quad_modified = 20.d0)
! eos2001 fit (value doesn't matter for relatively low density fit done at
!      this time [rho_lim = -0.5], but revisit later when want more stability
!      at higher density.)
!      parameter(quad_modified = 20.d0)
C        Revert back to 1.1.0 values.
      parameter(quad_modified = 10.d0)

C      pressure-ionization data for excited states.
      real*8 pi_fitx_neutral_ln_original, pi_fitx_ion_ln_original
      real*8 pi_fitx_neutral_ln, pi_fitx_ion_ln
      real*8 pi_fitx_neutral_ln_saumon, pi_fitx_ion_ln_saumon
      data pi_fitx_neutral_ln_original, pi_fitx_ion_ln_original
     &  /2*0.d0/
! pre eos2001 fit to original opal.
!      data pi_fitx_neutral_ln, pi_fitx_ion_ln/0.d0,-1.8d0/
! old fit to eos2001, start with Saumon
!      data pi_fitx_neutral_ln, pi_fitx_ion_ln/-4.d0,-40.d0/
! new fit to eos2001, start with solar considerations.  It turns out the
!      second fudge factor is important (changes of 0.003 in ln P) for
!      solar conditions near log T = 5.
!      data pi_fitx_neutral_ln, pi_fitx_ion_ln/0.d0,-1.8d0/
!      data pi_fitx_neutral_ln, pi_fitx_ion_ln/-4.d0,-40.d0/
C      Revert back to 1.1.0.
      data pi_fitx_neutral_ln, pi_fitx_ion_ln/0.d0,-1.8d0/
C       values which fit Saumon H table  (Saumon
C       He table calculated with no excited states)
      data pi_fitx_neutral_ln_saumon, pi_fitx_ion_ln_saumon
     &  /0.d0,-40.d0/
