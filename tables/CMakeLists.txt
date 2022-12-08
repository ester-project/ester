
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
# utils/private/utils_c_system.c

add_library(mesaUtilsC SHARED
    mesa/utils/private/utils_c_system.c)


    
add_library(mesaUtils SHARED
    mesa/utils/private/utils_system.f90
    mesa/utils/private/utils_openmp_stub.f90
    mesa/utils/private/utils_nan_sp.f90
    mesa/utils/private/utils_nan_qp.f90
    mesa/utils/private/utils_nan.f90
    mesa/utils/private/utils_nan_dp.f90
    mesa/utils/private/utils_idict.f90
    mesa/utils/private/utils_dict.f90
    mesa/utils/public/utils_lib.f90
    mesa/utils/public/utils_def.f90)
# 
# add_library(mesaUtilsPublic SHARED)
# 
# 

add_library(mesaConst SHARED
    mesa/const/public/const_lib.f90
    mesa/const/public/const_def.f90)
# 
# 
# 
# 
add_library(mesaMath SHARED
    mesa/math/private/math_pown.f90
    mesa/math/private/math_io.f90
    mesa/math/public/math_def.f90
#     mesa/math/public/math_lib_crmath.f90)
    mesa/math/public/math_lib_intrinsic.f90)


    
# add_library(mesaMtxL95DBL SHARED
# #     mesa/mtx/private/mtx_support.f90
#     mesa/mtx/private/my_lapack95.F90)
#     mesa/mtx/private/DGBSVX_wrapper.f90
#     mesa/mtx/private/bcyclic.f90
#     mesa/mtx/private/pre_conditioners.f90
#     mesa/mtx/public/mtx_def.f
#     mesa/mtx/public/mtx_lib.f90)


# add_library(mesaMtxL95QD SHARED
# #     mesa/mtx/private/mtx_support.f90
#     mesa/mtx/private/my_lapack95.F90)
    
add_library(mesaMtxf SHARED
#     mesa/mtx/private/mtx_support.f90
#     mesa/mtx/private/my_lapack95_quad.f90
#     mesa/mtx/private/my_lapack95_dble.f90
#     mesa/mtx/private/DGBSVX_wrapper.f90
#     mesa/mtx/private/bcyclic.f90
#     mesa/mtx/private/pre_conditioners.f90
    mesa/mtx/public/mtx_def.f)
#     mesa/mtx/public/mtx_lib.f90)   
    
add_library(mesaMtxf90 SHARED
    mesa/mtx/private/mtx_support.f90
    mesa/mtx/private/my_lapack95_quad.f90
    mesa/mtx/private/my_lapack95_dble.f90
    mesa/mtx/private/DGBSVX_wrapper.f90
    mesa/mtx/private/bcyclic.f90
    mesa/mtx/private/pre_conditioners.f90
#     mesa/mtx/public/mtx_def.f
    mesa/mtx/public/mtx_lib.f90)
    
# add_library(mesaMtx SHARED
#     mesa/mtx/private/mtx_support.f90
#     mesa/mtx/private/my_lapack95.F90
#     mesa/mtx/private/DGBSVX_wrapper.f90
# #     mesa/mtx/private/bcyclic.f90
#     mesa/mtx/private/pre_conditioners.f90
#     mesa/mtx/public/mtx_def.f
#     mesa/mtx/public/mtx_lib.f90)
    
    
add_library(mesaAutodiff SHARED
mesa/auto_diff/private/auto_diff_real_4var_order1_module.f90
mesa/auto_diff/private/auto_diff_real_2var_order3_module.f90
mesa/auto_diff/private/auto_diff_real_2var_order2_module.f90
mesa/auto_diff/private/auto_diff_real_2var_order1_module.f90
mesa/auto_diff/private/auto_diff_real_1var_order1_module.f90
mesa/auto_diff/private/auto_diff_real_15var_order1_module.f90
# mesa/auto_diff/private/auto_diff_quad_star_order1_module.f90.not_yet
# mesa/auto_diff/private/support_functions_quad.f90.not_yet
mesa/auto_diff/private/support_functions.f90
mesa/auto_diff/private/auto_diff_real_tdc_module.f90
mesa/auto_diff/private/auto_diff_real_star_order1_module.f90
mesa/auto_diff/public/auto_diff.f90)
    




add_library(mesaNumf90lib SHARED 
    mesa/num/private/mod_random.f90
    mesa/num/private/mod_qsort.f90
#     mesa/num/private/mod_integrate.f90
#     mesa/num/private/mod_simplex.f90
    mesa/num/private/mod_brent.f90)
#     mesa/num/public/accurate_sum_auto_diff_star_order1.f90
#     mesa/num/public/accurate_sum.f90   )


add_library(mesaNumflib SHARED 
    mesa/num/private/mod_rosenbrock.f
#     mesa/num/private/mod_root_x.f
    mesa/num/private/mod_root.f
    mesa/num/private/mod_newuoa.f
    mesa/num/private/mod_newton.f
    mesa/num/private/mod_dopri5.f
    mesa/num/private/mod_dop853.f
    mesa/num/private/mod_dc_decsol.f
    mesa/num/private/mod_bobyqa.f
    mesa/num/private/mod_solout_root.f
    mesa/num/public/num_def.f  )


add_library(mesaNumf90 SHARED 
#     mesa/num/private/mod_random.f90
#     mesa/num/private/mod_qsort.f90
    mesa/num/private/mod_simplex.f90
    mesa/num/private/mod_integrate.f90
#     mesa/num/private/mod_brent.f90
    mesa/num/public/accurate_sum_auto_diff_star_order1.f90
    mesa/num/public/accurate_sum.f90   )
add_library(mesaNumf SHARED 
#     mesa/num/private/mod_random.f90
#     mesa/num/private/mod_qsort.f90
#     mesa/num/private/mod_simplex.f90
    mesa/num/public/num_lib.f)

 
add_library(mesaHdf5io SHARED 
mesa/hdf5io/private/kinds_m.fypp
mesa/hdf5io/private/kinds_m.f90
mesa/hdf5io/private/hdf5io_m.fypp
mesa/hdf5io/private/hdf5io_m.f90
# mesa/hdf5io/private/hdf5io.inc
# mesa/hdf5io/private/orum.inc
mesa/hdf5io/private/cinter_m.fypp
mesa/hdf5io/private/cinter_m.f90
mesa/hdf5io/public/hdf5io_lib.f90)


# 
# 
add_library(mesaChem SHARED
    mesa/chem/private/nuclide_set_mod.f90
    mesa/chem/private/lodders_mod.f90
    mesa/chem/private/chem_isos_io.f90
    mesa/chem/public/chem_lib.f90
    mesa/chem/public/chem_def.f90)

add_library(mesaInterp1d SHARED
    mesa/interp_1d/private/interp_1d_pm_sg.f90
    mesa/interp_1d/private/interp_1d_pm.f90
    mesa/interp_1d/private/interp_1d_mp_sg.f90
    mesa/interp_1d/private/interp_1d_mp.f90
    mesa/interp_1d/private/interp_1d_misc_sg.f90
    mesa/interp_1d/private/interp_1d_misc.f90
    mesa/interp_1d/public/interp_1d_lib_sg.f90
    mesa/interp_1d/public/interp_1d_def.f90
    mesa/interp_1d/public/interp_1d_lib.f90)


add_library(mesaInterp2df SHARED
    mesa/interp_2d/private/renka790_sg.f
    mesa/interp_2d/private/renka790_db.f
    mesa/interp_2d/private/bipm_sg.f
    mesa/interp_2d/private/bipm_db.f
    mesa/interp_2d/private/bicub_sg.f
    mesa/interp_2d/private/bicub_db.f
    mesa/interp_2d/private/akima760_sg.f
    mesa/interp_2d/private/akima760_db.f)

add_library(mesaInterp2df90 SHARED
    mesa/interp_2d/public/interp_2d_lib_sg.f90
    mesa/interp_2d/public/interp_2d_lib_db.f90)

add_library(mesaKapf90lib SHARED
    mesa/kap/other/other_radiative_opacity.f90
    mesa/kap/other/other_elect_cond_opacity.f90
    mesa/kap/other/other_compton_opacity.f90
    mesa/kap/private/op_def.f90
    mesa/kap/private/load_kap.f90
    mesa/kap/private/load_co_kap.f90
    mesa/kap/private/kap_eval_support.f90
    mesa/kap/private/kap_eval_fixed.f90
    mesa/kap/private/kap_eval.f90
    mesa/kap/private/kap_eval_co.f90
    mesa/kap/private/kap_ctrls_io.f90
    mesa/kap/private/kapcn.f90
    mesa/kap/private/kap_aesopus.f90
    mesa/kap/private/condint.f90
#     mesa/kap/public/kap_lib.f90
    mesa/kap/public/kap_def.f90)
    
add_library(mesaKapf SHARED
    mesa/kap/private/op_radacc.f
    mesa/kap/private/op_osc.f
    mesa/kap/private/op_load.f
    mesa/kap/private/op_ev.f
    mesa/kap/private/op_eval.f
    mesa/kap/private/op_common.f)
add_library(mesaKapf90 SHARED
    mesa/kap/public/kap_lib.f90
    mesa/variables_ester_and_mesa.f90)

add_library(mesaEos SHARED
    mesa/eos/other/other_eos.f90
    mesa/eos/private/create_EXCOR7_table.f90
    mesa/eos/private/skye_thermodynamics.f90
    mesa/eos/private/skye.f90
    mesa/eos/private/skye_ideal.f90
    mesa/eos/private/skye_coulomb_solid.f90
    mesa/eos/private/skye_coulomb_liquid.f90
    mesa/eos/private/skye_coulomb.f90
    mesa/eos/private/pc_support.f90
    mesa/eos/private/pc_eos.f90
    mesa/eos/private/ion_offset.f90
    mesa/eos/private/ideal.f90
    mesa/eos/private/helm_polynomials.f90
    mesa/eos/private/helm.f90
    mesa/eos/private/helm_alloc.f90
    mesa/eos/private/gauss_fermi.f90
    mesa/eos/private/eospt_eval.f90
    mesa/eos/private/eospc_eval.f90
    mesa/eos/private/eos_initialize.f90
    mesa/eos/private/eos_helm_eval.f90
    mesa/eos/private/eosdt_support.f90
    mesa/eos/private/eosdt_load_tables.f90
    mesa/eos/private/eosdt_eval.f90
    mesa/eos/private/eosde_eval.f90
    mesa/eos/private/eos_ctrls_io.f90
    mesa/eos/private/eoscms_eval.f90
    mesa/eos/private/eos_blend.f90
    mesa/eos/private/create_FSCRliq8_table.f90
    mesa/eos/public/eos_lib.f90
    mesa/eos/public/eos_def.f90)




SET_TARGET_PROPERTIES(mesaUtils PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaConst PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaMath PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")


# SET_TARGET_PROPERTIES(mesaMtxL95DBL PROPERTIES COMPILE_FLAGS "-std=f2008 -cpp -DDBLE -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")
# 
# 
# SET_TARGET_PROPERTIES(mesaMtxL95QD PROPERTIES COMPILE_FLAGS "-std=f2008 -cpp -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaMtxf PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffixed-form -ffixed-line-length-none -x f77-cpp-input")

SET_TARGET_PROPERTIES(mesaMtxf90 PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaAutodiff PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaNumf90lib PROPERTIES COMPILE_FLAGS "-std=legacy -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination  -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaNumflib PROPERTIES COMPILE_FLAGS "-std=legacy -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -O3 -ftree-vectorize -c -ffixed-form -ffixed-line-length-none -x f77-cpp-input")

SET_TARGET_PROPERTIES(mesaNumf90 PROPERTIES COMPILE_FLAGS "-std=legacy -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan  -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination  -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaNumf PROPERTIES COMPILE_FLAGS "-std=legacy -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -O3 -ftree-vectorize -c -ffixed-form -ffixed-line-length-none -x f77-cpp-input")


SET_TARGET_PROPERTIES(mesaChem PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaInterp1d PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaInterp2df PROPERTIES COMPILE_FLAGS "-std=legacy  -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -O3 -ftree-vectorize -c -ffixed-form -ffixed-line-length-none -x f95-cpp-input ")

SET_TARGET_PROPERTIES(mesaInterp2df90 PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffixed-form -ffixed-line-length-none -x f95-cpp-input")

SET_TARGET_PROPERTIES(mesaHdf5io PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -lhdf5 -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")


SET_TARGET_PROPERTIES(mesaKapf90lib PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")


SET_TARGET_PROPERTIES(mesaKapf PROPERTIES COMPILE_FLAGS "-std=legacy -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -O3 -ftree-vectorize -c -ffixed-form -ffixed-line-length-none -x f77-cpp-input")



SET_TARGET_PROPERTIES(mesaKapf90 PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")


SET_TARGET_PROPERTIES(mesaEos PROPERTIES COMPILE_FLAGS "-std=f2008 -Wno-uninitialized -fno-range-check -fmax-errors=7 -fprotect-parens -fno-sign-zero -fbacktrace -ggdb -finit-real=snan -Wno-error=tabs -Wunused-value -W -Wno-compare-reals -Wno-unused-parameter -Wno-function-elimination -fimplicit-none -O3 -ftree-vectorize -c -ffree-form -ffree-line-length-none -x f95-cpp-input")

# target_compile_options( 
#      mesaUtils PRIVATE "-std=legacy")
# 
# target_compile_options(
#     mesaConst PRIVATE "-std=legacy")
# 
# target_compile_options(
#     mesaMath PRIVATE "-std=legacy")
# 
# target_compile_options(
#     mesaNum PRIVATE "-std=legacy")
# 
#     
# target_compile_options(
#     mesaChem PRIVATE "-std=legacy")
# 
#     
# target_compile_options(
#     mesaKap PRIVATE "-std=legacy")
# 
# target_include_directories(mesaNumf PRIVATE ${mesa/num/public})
# target_include_directories(mesaNumf90 PRIVATE ${mesa/num/public})
    
add_custom_target(opal ALL
    COMMAND tar xf ${CMAKE_SOURCE_DIR}/tables/opal/eos_tables.tar.gz
    )
add_custom_target(houdek ALL
    COMMAND
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/opal95e.bin           > OPINTPATH_AX &&
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/alex94.bin            >> OPINTPATH_AX &&
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/opal95e+alex94-pd.bin >> OPINTPATH_AX &&
    echo ${CMAKE_INSTALL_PREFIX}/share/ester/tables/houdek/v9/ival95.dat            >> OPINTPATH_AX)

# set_source_files_properties(
#   my_lapack95.F90 
#   PROPERTIES Fortran_PREPROCESS ON
# )
target_link_libraries(mesaUtils mesaConst mesaUtilsC)
target_link_libraries(mesaMath mesaConst mesaUtils)
target_link_libraries(mesaAutodiff mesaMath mesaConst mesaUtils)

target_link_libraries(mesaMtxf mesaAutodiff mesaConst )
target_link_libraries(mesaMtxf90 mesaMtxf mesaAutodiff mesaConst )


target_link_libraries(mesaNumf90lib mesaAutodiff mesaConst mesaMtxf90 mesaMtxf )

target_link_libraries(mesaNumflib mesaNumf90lib mesaMath mesaConst mesaUtils mesaMtxf90 mesaMtxf )

target_link_libraries(mesaNumf90 mesaNumf90lib mesaNumflib mesaMath mesaConst mesaUtils mesaMtxf90 mesaMtxf )
target_link_libraries(mesaNumf mesaNumf90 mesaNumflib mesaNumf90lib mesaMath mesaConst mesaUtils mesaMtxf90 mesaMtxf )


target_link_libraries(mesaChem mesaNumf90 mesaNumf mesaAutodiff mesaMath mesaConst mesaUtils)
target_link_libraries(mesaInterp1d mesaMath mesaConst mesaUtils)

target_link_libraries(mesaInterp2df mesaInterp1d mesaNumf mesaMath mesaConst mesaUtils)
target_link_libraries(mesaInterp2df90 mesaInterp2df mesaMath mesaConst mesaUtils)

target_link_libraries(mesaHdf5io ${HDF5_LIBRARIES})

target_link_libraries(mesaKapf90lib mesaInterp2df90 mesaEos mesaHdf5io mesaAutodiff mesaChem mesaNumf90 mesaNumf mesaMath mesaConst mesaUtils)

target_link_libraries(mesaKapf mesaKapf90lib mesaEos mesaChem mesaNumf90 mesaNumf mesaMath mesaConst mesaUtils)

target_link_libraries(mesaKapf90 mesaKapf mesaKapf90lib mesaEos mesaChem mesaNumf90 mesaNumf mesaMath mesaConst mesaUtils)

target_link_libraries(mesaEos mesaChem mesaNumf90 mesaNumf mesaMath mesaConst mesaUtils mesaInterp2df90 mesaInterp2df)

include_directories(mesa/num/public)
include_directories(mesa/num/private)
include_directories(mesa/eos/private)
include_directories(mesa/eos/defaults)
include_directories(mesa/utils/)
include_directories(mesa/kap/defaults)
include_directories(mesa/hdf5io/private)
include_directories(mesa/mtx/public)


# include_directories(/home/mkessar/mesasdk/include)
# link_directories(/home/mkessar/mesasdk/lib)

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

install(TARGETS cesam freeeos opint mesaKapf mesaKapf90 mesaKapf90lib mesaChem mesaNumf90 mesaNumf mesaMath mesaConst mesaUtils mesaInterp2df90 mesaHdf5io mesaAutodiff mesaEos mesaInterp1d mesaInterp2df mesaUtilsC mesaNumflib mesaNumf90lib mesaMtxf90 mesaMtxf
    DESTINATION lib)
