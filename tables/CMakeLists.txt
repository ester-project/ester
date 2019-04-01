add_library(cesam SHARED
    cesam/mod_variables.f
    cesam/mod_kind.f
    cesam/mod_numerique.f
    cesam/test_opa.f
    cesam/mod_donnees.f
    cesam/mod_nuc.f
    cesam/mod_opa.f
    cesam/sortie.f)

add_library(freeeos SHARED
    freeeos/aux_to_traditional.f
    freeeos/bfgs.f
    freeeos/coulomb_adjust.f
    freeeos/coulomb.f
    freeeos/deriv_product.f
    freeeos/deriv_ratio.f
    freeeos/effective_radius.f
    freeeos/effsum_calc.f
    freeeos/eos_bfgs.f
    freeeos/eos_calc.f
    freeeos/eos_cold_start.f
    freeeos/eos_free_calc.f
    freeeos/eos_jacobian.f
    freeeos/eos_sum_calc.f
    freeeos/eos_tqft.f
    freeeos/eos_tqsum_calc.f
    freeeos/eos_warm_step.f
    freeeos/exchange_coeff.f
    freeeos/exchange_gcpf.f
    freeeos/excitation_pi_end.f
    freeeos/excitation_pi.f
    freeeos/excitation_sum.f
    freeeos/exct_calc.f
    freeeos/fermi_dirac_coeff.f
    freeeos/fermi_dirac_ct.f
    freeeos/fermi_dirac_direct.f
    freeeos/fermi_dirac_exchange.f
    freeeos/fermi_dirac.f
    freeeos/fermi_dirac_minusnr_coeff.f
    freeeos/fermi_dirac_original_coeff.f
    freeeos/fermi_dirac_recursion.f
    freeeos/fjs_pi.f
    freeeos/fp12_calc.f
    freeeos/f_psi.f
    freeeos/free_eos_detailed.f
    freeeos/free_eos.f
    freeeos/free_eos_interpolated.f
    freeeos/free_non_ideal_calc.f
    freeeos/ionize.f
    freeeos/lnblnk.f
    freeeos/lnx_calc.f
    freeeos/Makefile.am
    freeeos/master_coulomb.f
    freeeos/mdh_pi.f
    freeeos/molecular_hydrogen.f
    freeeos/pl_prepare.f
    freeeos/plsum_approx.f
    freeeos/plsum.f
    freeeos/poly_sum.f
    freeeos/pteh_pi.f
    freeeos/pteh_theta.f
    freeeos/qmhd_calc.f
    freeeos/qryd_approx.f
    freeeos/qryd_calc.f
    freeeos/qstar_calc.f
    freeeos/quad4b.f
    freeeos/round_ln.f
    freeeos/tau_calc.f
    freeeos/version.f)

add_library(opint SHARED
    houdek/v9/lib/condux_ad.f
    houdek/v9/lib/dopinit.f
    houdek/v9/lib/dreadkx.f
    houdek/v9/lib/ratmat.f
    houdek/v9/lib/condux_cf.f
    houdek/v9/lib/dopintc.f
    houdek/v9/lib/dshep.f
    houdek/v9/lib/rbival.f
    houdek/v9/lib/dadjac.f
    houdek/v9/lib/dopintf.f
    houdek/v9/lib/duvip3d.f
    houdek/v9/lib/rbivpd.f
    houdek/v9/lib/dctang.f
    houdek/v9/lib/dopints.f
    houdek/v9/lib/duvip3p.f
    houdek/v9/lib/rde95.f
    houdek/v9/lib/dextrp.f
    houdek/v9/lib/dordgr.f
    houdek/v9/lib/gi.f
    houdek/v9/lib/rdi95.f
    houdek/v9/lib/dhunt.f
    houdek/v9/lib/dordgri.f
    houdek/v9/lib/gid.f
    houdek/v9/lib/readal.f
    houdek/v9/lib/dintrp.f
    houdek/v9/lib/dpderix.f
    houdek/v9/lib/inttwo.f
    houdek/v9/lib/rtperm.f
    houdek/v9/lib/dintrpi.f
    houdek/v9/lib/dpderiy.f
    houdek/v9/lib/itoh_coeff_cf.f
    houdek/v9/lib/sbilin.f
    houdek/v9/lib/dmaceps.f
    houdek/v9/lib/dpdmin.f
    houdek/v9/lib/itoh_ec_ad.f
    houdek/v9/lib/sfa_ad.f
    houdek/v9/lib/dmasub.f
    houdek/v9/lib/dpdste.f
    houdek/v9/lib/itoh_ec_cf.f
    houdek/v9/lib/sfa_cf.f
    houdek/v9/lib/dmasube.f
    houdek/v9/lib/dreade.f
    houdek/v9/lib/itoh_init_ad.f
    houdek/v9/lib/trdisa.f
    houdek/v9/lib/dmasubi.f
    houdek/v9/lib/dreadi.f
    houdek/v9/lib/itoh_init_cf.f
    houdek/v9/lib/trdisb.f
    houdek/v9/lib/dmaxmn.f
    houdek/v9/lib/dreadk.f
    houdek/v9/lib/rat2d.f
    houdek/v9/lib/dcspline.f)

add_custom_target(opal ALL
    COMMAND tar xf ${CMAKE_SOURCE_DIR}/tables/opal/eos_tables.tar.gz
    )
add_custom_target(houdek ALL
    COMMAND
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/opal95e.bin           > OPINTPATH_AX &&
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/alex94.bin            >> OPINTPATH_AX &&
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/opal95e+alex94-pd.bin >> OPINTPATH_AX &&
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/ival95.dat            >> OPINTPATH_AX)


# Opal tables
install(FILES
    ${CMAKE_SOURCE_DIR}/tables/opal/GN93hz
    DESTINATION share/ester/tables/opal)

# Cesam tables
install(FILES
    ${CMAKE_SOURCE_DIR}/tables/cesam/opa/opa_yveline.data
    DESTINATION share/ester/tables/cesam)

# Houdek tables
install(FILES
    ${CMAKE_SOURCE_DIR}/tables/houdek/v9/alex94/alex94.bin
    ${CMAKE_SOURCE_DIR}/tables/houdek/v9/opal95/ival95.dat
    ${CMAKE_SOURCE_DIR}/tables/houdek/v9/opal95/opal95e+alex94-pd.bin
    ${CMAKE_SOURCE_DIR}/tables/houdek/v9/opal95/opal95e.bin
    ${CMAKE_CURRENT_BINARY_DIR}/OPINTPATH_AX
    DESTINATION share/ester/tables/houdek/v9)


install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/eos_tables
    DESTINATION share/ester/tables/opal)

install(TARGETS cesam freeeos opint
    DESTINATION lib)
