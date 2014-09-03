#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ESTER], [1.01])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])

AM_PROG_AR
LT_INIT

AM_INIT_AUTOMAKE([-Wall subdir-objects foreign])
AC_CONFIG_HEADERS([config.h])
# AC_CONFIG_SRCDIR([src/])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_FC
AC_PROG_F77($FC)
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
have_pgplot="no"
AC_CHECK_HEADERS([cpgplot.h],
                 [AC_DEFINE([USE_PGPLOT], [1], [enable pgplot support])
                 have_pgplot="yes"
                 AC_CHECK_LIB([cpgplot], [cpgebuf], [],
                              [AC_ERROR(cpgplot library not found)])])
AC_CHECK_HEADERS([stddef.h stdlib.h string.h sys/time.h])
AC_CHECK_HEADERS([cblas.h lapack.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
# Search for mkl
AC_SEARCH_LIBS([cblas_dgemm], [mkl_intel_lp64],
               [have_mkl=yes
                LIBS="$LIBS -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"],
               [have_mkl=no],
               [-lmkl_intel_thread -lmkl_core -liomp5 -lpthread])

AS_IF([test "x$have_mkl" == "xyes"],
      [BLAS_LIB="MKL"
       AC_MSG_NOTICE(blas library: mkl)],
      [AC_SEARCH_LIBS([cblas_sgemm], [cblas blas],
                      [],
                      [AC_ERROR(cblas library not found)])
      AC_CHECK_LIB([blas], [dgemm_], [],
                   [AC_ERROR(blas library not found: consider setting LDFLAGS)])
      AC_CHECK_LIB([lapack], [dgetrf_], [],
                   [AC_ERROR(lapack lib not found: consider setting LDFLAGS)])])

AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strtoul])

AC_CONFIG_FILES([Makefile
                 tables/Makefile
                 tables/cesam_nuc/Makefile
                 tables/houdek/v9/Makefile
                 tables/opal/Makefile
                 src/Makefile
                 src/main/ester])

AS_IF([test "x$ac_cv_fc_compiler_gnu" = "xyes"], [
       FFLAGS="$FFLAGS -fno-align-commons"
       FFLAGS="$FFLAGS -fdefault-real-8"
       FFLAGS="$FFLAGS -fdefault-double-8"
       FFLAGS="$FFLAGS -ffixed-line-length-132"
       FFLAGS="$FFLAGS -O -c -w"
       FFLAGS="$FFLAGS -fd-lines-as-comments"], [
       FFLAGS="$FFLAGS -r8 -nowarn -fpp -extend-source"])

AC_OUTPUT
AC_MSG_NOTICE()
AC_MSG_NOTICE(config:)
AC_MSG_NOTICE(FC:       $FC)
AC_MSG_NOTICE(CC:       $CC)
AC_MSG_NOTICE(LIBS:     $LIBS)
AC_MSG_NOTICE(PGPLOT:   $have_pgplot)

cat<< EOF

-------------------------------------------------------------------------------
ester build configuration:

Fortran compiler: $FC
C compiler:       $CC
C++ compiler:     $CXX
pgplot enabled:   $have_pgplot
found MKL:        $have_mkl
LIBS:             $LIBS
prefix:           $prefix
-------------------------------------------------------------------------------

EOF
