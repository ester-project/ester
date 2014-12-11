#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([ESTER], [1.0.108])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.9 no-dist-gzip dist-bzip2 tar-ustar -Wall subdir-objects foreign])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CC
AC_PROG_FC
AC_PROG_F77($FC)
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

AM_PATH_PYTHON

# Checks for libraries.

# Checks for header files.
have_pgplot="no"
AC_CHECK_HEADERS([cpgplot.h],
                 [AC_DEFINE([USE_PGPLOT], [1], [enable pgplot support])
                 AC_CHECK_LIB([cpgplot], [cpgebuf], [
                               LIBS="-lcpgplot $LIBS" have_pgplot="yes"],
                              [AC_ERROR(cpgplot library not found)])])

AC_CHECK_HEADERS([stddef.h stdlib.h string.h sys/time.h])
AC_CHECK_HEADERS([cblas.h lapack.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
# AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_LIB([m], [sqrt])

AC_ARG_WITH([hdf5-root],
            AS_HELP_STRING([--with-hdf5-root],
                           [Provide hdf5 root directory]),
                           hdf5root=$withval)

AC_ARG_ENABLE([hdf5],
              AS_HELP_STRING([--enable-hdf5],
                             [Enable HDF5 support]))
AS_IF([test "x$enable_hdf5" == "xyes"], [
       CFLAGS="$CFLAGS -I$hdf5root/include"
       CXXFLAGS="$CXXFLAGS -I$hdf5root/include"
       AC_CHECK_HEADERS([hdf5.h])
       AS_IF([test "x$ac_cv_header_hdf5_h" == "xno"], [
              AC_ERROR([Could not find HDF5 header
    Please privde the path to HDF5 with --with-hdf5-root])])
       LIBS="$LIBS -L$hdf5root/lib"
       AC_CHECK_LIB([hdf5], [H5Fopen])
       AS_IF([test "x$ac_cv_lib_hdf5_H5Fopen" == "xno"], [
              AC_ERROR([Could not find HDF5 library
    Please privde the path to HDF5 with --with-hdf5-root])])
       ], [enable_hdf5="no"])


# Search for mkl
AC_SEARCH_LIBS([cblas_dgemm], [mkl_intel_lp64],
               [have_mkl=yes
                AC_DEFINE([USE_MKL], [1],
                          [Use MKL lib])
                LIBS="$LIBS -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"],
               [have_mkl=no],
               [-lmkl_intel_thread -lmkl_core -liomp5 -lpthread])

AS_IF([test "x$have_mkl" == "xyes"],
      [BLAS_LIB="MKL"], [
       AC_SEARCH_LIBS([cblas_sgemm], [cblas blas],
                      [],
                      [AC_ERROR(cblas library not found)])
       AC_CHECK_LIB([blas], [dgemm_], [],
                    [AC_ERROR(blas library not found: consider setting LDFLAGS)])
       AC_CHECK_LIB([lapack], [dgetrf_], [],
                    [AC_ERROR(lapack library not found: consider setting LDFLAGS)])])

AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strtoul])

AC_CONFIG_FILES([Makefile
                 tables/Makefile
                 tables/cesam_nuc/Makefile
                 tables/houdek/v9/lib/Makefile
                 tables/opal/Makefile
                 src/Makefile
                 src/physics/Makefile
                 test/models/Makefile])

AC_LANG_PUSH([Fortran 77])
AC_MSG_CHECKING([whether $FC supports -fno-align-commons])
save_fflags=$FFLAGS
FFLAGS="-fno-align-commons"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
                  [AC_MSG_RESULT([yes])]
                  [fc_support_no_align_commons=yes],
                  [AC_MSG_RESULT([no])]
                  [fc_support_no_align_commons=no])
FFLAGS=$save_fflags
AC_LANG_POP([Fortran 77])

AS_IF([test "x$ac_cv_fc_compiler_gnu" = "xyes"], [
       AS_IF([test "x$fc_support_no_align_commons" = "xyes"],
             [FFLAGS="$FFLAGS -fno-align-commons"
              FFLAGS="$FFLAGS -fdefault-real-8"
              FFLAGS="$FFLAGS -fdefault-double-8"
              FFLAGS="$FFLAGS -ffixed-line-length-132"
              FFLAGS="$FFLAGS -O -c -w"
              FFLAGS="$FFLAGS -fd-lines-as-comments"],
              [AC_ERROR($FC does not support -fno-align-commons)])],
       [FFLAGS="$FFLAGS -r8 -nowarn -fpp -extend-source"])

AC_ARG_ENABLE([perf-log],
              AS_HELP_STRING([--enable-perf-log],
                             [Enable performance monitoring]), [
              perf_log=$enableval], [
              perf_log="no"])
AS_IF([test "x$enable_perf_log" == "xyes"], [
       AC_DEFINE([PERF_LOG], [1], [Should perf-log be enabled])
       perf_log="yes"])

AC_ARG_ENABLE([debug],
              AS_HELP_STRING([--enable-debug],
                             [Enable debug mode]), [
              AC_DEFINE([DEBUG], [1], [Eanble debug mode])
              CFLAGS="$CFLAGS -g"
              CXXFLAGS="$CXXFLAGS -g"])

AC_OUTPUT

cat<< EOF

-------------------------------------------------------------------------------
ester build configuration:

Fortran compiler: $FC
C compiler:       $CC
C++ compiler:     $CXX
pgplot enabled:   $have_pgplot
HDF5 enabled:     $enable_hdf5
found MKL:        $have_mkl
perf_log:         $perf_log
CXXFLAGS:         $CXXFLAGS
LIBS:             $LIBS
prefix:           $prefix
-------------------------------------------------------------------------------

EOF
