#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([ESTER], [1.0.108])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.9 dist-bzip2 dist-zip tar-ustar -Wall subdir-objects foreign])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_FC(gfortran44 gfortran ifort)
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
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([cpgplot.h],
                 [AC_DEFINE([USE_PGPLOT], [1], [enable pgplot support])
                 AC_CHECK_LIB([cpgplot], [cpgebuf],
                              [LIBS="-lcpgplot $LIBS" have_pgplot="yes"],
                              [AC_ERROR(cpgplot library not found)])])

AC_CHECK_HEADERS([stddef.h stdlib.h string.h sys/time.h])
AC_CHECK_HEADERS([cblas.h lapack.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
# AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_LIB([m], [sqrt])

AC_ARG_ENABLE([hdf5],
              AS_HELP_STRING([--enable-hdf5],
                             [Enable HDF5 support]))
AS_IF([test "x$enable_hdf5" == "xyes"],
      [CXX="h5c++"
       AC_MSG_CHECKING([if HDF5 compiler works])
       AC_COMPILE_IFELSE(
                      [AC_LANG_SOURCE([
                       int main(void) {}
                       ])],
                       [h5cxx_works=yes],
                       [h5cxx_works=no])
       AC_MSG_RESULT($h5cxx_works)
       AS_IF([test "x$h5cxx_works" == "xno"],
             [AC_ERROR(HDF5 compiler wrapper ($CXX) does not work)])
       AC_CHECK_HEADERS([H5Cpp.h])
       AC_DEFINE([USE_HDF5], [1], [Enable hdf5 support])
       AS_IF([test "x$ac_cv_header_hdf5_h" == "xno"], [
              AC_ERROR([Could not find HDF5 header hdf5.h])])
       ],
   [enable_hdf5="no"])

AM_CONDITIONAL(USE_HDF5, [test "x$enable_hdf5" == "xyes"])

AC_CHECK_FUNC(cblas_dgemm, [have_cblas=yes],  [have_cblas=no])
AC_CHECK_FUNC(dgemm_     , [have_blas=yes],   [have_blas=no])
AC_CHECK_FUNC(dgetrf_    , [have_lapack=yes], [have_lapack=no])

# Search for mkl
AC_MSG_CHECKING([for mkl lib])
AC_COMPILE_IFELSE(
                  [AC_LANG_SOURCE([
                   #include <mkl_cblas.h>
                   void main(void) {}
                   ])],
                   [have_mkl=yes],
                   [have_mkl=no])
AC_MSG_RESULT($have_mkl)

AS_IF([test "x$have_cblas" == "xno" -a "x$have_blas" == "xno" -a "x$have_lapack" == "xno"],
      [AC_SEARCH_LIBS([dgetrf_], [mkl_intel_lp64],
               [have_mkl=yes
                AC_DEFINE([USE_MKL], [1],
                          [Use MKL lib])
                LIBS="$LIBS -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"],
               [have_mkl=no],
               [-lmkl_intel_thread -lmkl_core -liomp5 -lpthread])])

AS_IF([test "x$have_mkl" == "xno"],
      [AS_IF([test "x$ac_cv_header_cblas_h" == "xno"],
             [AC_ERROR(could not find cblas.h)])
       AS_IF([test "x$have_cblas" == "xno"],
             [AC_SEARCH_LIBS([cblas_sgemm], [cblas blas],
                             [],
                             [AC_ERROR(cblas library not found)])])
       AS_IF([test "x$have_blas" == "xno"],
             [AC_CHECK_LIB([blas], [dgemm_], [],
                           [AC_ERROR(blas library not found)])])
       AS_IF([test "x$have_lapack" == "xno"],
             [AC_CHECK_LIB([lapack], [dgetrf_], [],
                           [AC_ERROR(blas library not found)])])])

AC_LANG_POP([C++])

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

AC_DEFINE([ESTER_DATADIR], ["@datadir@"], [Path to ester data directory])

AC_OUTPUT

cat<< EOF

-------------------------------------------------------------------------------
ester build configuration:

FC:         $FC
CXX:        $CXX
CXXFLAGS:   $CXXFLAGS
LIBS:       $LIBS
pgplot:     $have_pgplot
HDF5:       $enable_hdf5
perf_log:   $perf_log
prefix:     $prefix
-------------------------------------------------------------------------------

EOF
