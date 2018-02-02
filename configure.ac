#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
m4_include([version.m4])

AC_INIT([ESTER],
        [ESTER_VERSION],
        [http://ester-project.github.io/ester/],
        [ester])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([1.9 no-dist-gzip dist-bzip2 tar-ustar -Wall subdir-objects foreign])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])

AC_PREFIX_DEFAULT(${HOME}/local)

# Checks for programs.
AC_PROG_LIBTOOL
AC_PROG_CPP
AC_PROG_CXX([icpc g++ clang++])
AC_PROG_CC([icc gcc clang])
AC_PROG_FC([ifort gfortran44 gfortran])
AC_PROG_F77([$FC])
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_CHECK_PROG([ADDR2LINE], [addr2line], [yes], [no])
CXXFLAGS="$CXXFLAGS -std=c++11"

AS_IF([test "x$ADDR2LINE" == "xyes"],
      [AC_DEFINE([HAVE_ADDR2LINE], [1], [addr2line is installed])])

AM_PATH_PYTHON(2.3)
AX_PKG_SWIG([1.3.17], [have_swig="yes"], [have_swig="no"])

AM_CONDITIONAL(HAVE_SWIG, [test "x$have_swig" == "xyes"])

AC_ARG_ENABLE([plt],
              AS_HELP_STRING([--disable-plt],
                             [Enable/disable matplotlib support]))

AC_ARG_VAR([PYTHON_CONFIG], [Path to python-config])
AC_PATH_PROGS([PYTHON_CONFIG],
              [python$PYTHON_VERSION-config python-config],
              [no],
              [`dirname $PYTHON`])
AS_IF([test "x$PYTHON" != "x:"],
      [PYTHONINC=`$PYTHON -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"`
       CPPFLAGS="$CPPFLAGS -I$PYTHONINC"
       PYTHONLIBS=`$PYTHON_CONFIG --ldflags`
       LIBS+=" $PYTHONLIBS"
       AC_DEFINE([PYTHON_VERSION], ["@PYTHON_VERSION@"], [Version of python])])

AS_IF([test "x$enable_plt" != "xno"],
    [enable_plt="yes"
     AC_MSG_CHECKING([for matplotlib GTKAgg backend])
     have_gtk=`$PYTHON -c "import sys; \
             sys.path.insert(0, '/usr/lib/python$PYTHON_VERSION/dist-packages'); \
             import matplotlib; matplotlib.use('GTKAgg'); import matplotlib.pyplot; print('yes')" 2> /dev/null`

     AS_IF([test "x$have_gtk" == "xyes"],
           [AC_DEFINE([USE_GTK], [1], [Use GTKAgg matplotlib backend])
            AC_MSG_RESULT(yes)],
           [AC_MSG_RESULT(no)])

     AC_MSG_CHECKING([for numpy API version])
     AC_COMPILE_IFELSE(
                       [AC_LANG_SOURCE([
                        #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
                        #include <numpy/arrayobject.h>
                        void main(void) {
                        int i = NPY_ARRAY_FARRAY_RO;
                        }
                        ])],
                        [numpy_api=1.7],
                        [numpy_api=deprecated])
     AC_MSG_RESULT($numpy_api)]
     AC_DEFINE([ENABLE_PLT], [1], [Enable matplotlib grahpical output]),
    [enable_plt="no"])

AS_IF([test "x$numpy_api" == "xdeprecated"],
      [AC_DEFINE([USE_DEPRECATED_NUMPY], [1],
                 [Use deprecated (<1.7) numpy A API])])


# Checks for libraries.

# Checks for header files.
AC_LANG_PUSH([C++])

AC_CHECK_HEADERS([stddef.h stdlib.h string.h sys/time.h])
AC_CHECK_HEADERS([cblas.h lapack.h])
AC_CHECK_HEADERS([Python.h],
                 [pymodule=yes],
                 [pymodule=no])

AM_CONDITIONAL(BUILD_ESTER_PY, [test "x$pymodule" == "xyes"])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
# AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_LIB([m], [sqrt])

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
      [AC_SEARCH_LIBS([dgetrf_], [mkl_rt],
               [have_mkl=yes
                AC_DEFINE([USE_MKL], [1],
                          [Use MKL lib])
                LIBS="$LIBS -lmkl_rt"],
               [have_mkl=no])])

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

AC_ARG_ENABLE([hdf5],
              AS_HELP_STRING([--enable-hdf5],
                             [Enable HDF5 support]))

AX_LIB_HDF5([serial])
extension="out"
AS_IF([test "x$enable_hdf5" != "xno"],
      [AC_MSG_CHECKING([whether we can link against hdf5 libraries])
      AS_IF([test "x$HDF5_LIBS" != "x"],
            [LIBS="$LIBS $HDF5_LIBS -lhdf5_cpp"
             CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
             LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"])
       AC_LINK_IFELSE(
                      [AC_LANG_SOURCE([
                       #include <H5Cpp.h>
                       int main(void) {}
                       ])],
                       [h5cxx_works="yes"],
                       [h5cxx_works="no"])
       AC_MSG_RESULT($h5cxx_works)
       AC_MSG_CHECKING([whether h5cc works])
       AC_LANG_PUSH([C])
       AC_LINK_IFELSE(
                      [AC_LANG_SOURCE([
                       #include <hdf5.h>
                       int main(void) {}
                       ])],
                       [h5cc_works="yes"],
                       [h5cc_works="no"])
       AC_LANG_POP([C])
       AC_MSG_RESULT($h5cc_works)
       AS_IF([test "x$h5cxx_works" == "xno" -o "x$h5cc_works" == "xno"],
             [enable_hdf5="no"
              extension="out"],
             [enable_hdf5="yes"
              extension="h5"
              AC_DEFINE([USE_HDF5], [1], [Enable hdf5 support]) ])],
   [enable_hdf5="no"])

AM_CONDITIONAL(USE_HDF5, [test "x$enable_hdf5" == "xyes"])

AC_SUBST(extension)

AC_LANG_POP([C++])

AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strtoul])

AC_CONFIG_FILES([Makefile
                 tables/Makefile
                 tables/cesam/Makefile
                 tables/houdek/v9/Makefile
                 tables/houdek/v9/lib/Makefile
                 tables/opal/Makefile
                 src/Makefile
                 src/physics/Makefile
                 lib-example/Makefile
                 test/models/Makefile
                 python/Makefile])

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

AS_IF([test "x$ac_cv_fc_compiler_gnu" == "xyes"], [
       AS_IF([test "x$fc_support_no_align_commons" == "xyes"],
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

AC_MSG_CHECKING([whether we can run compiled binaries])
AC_LINK_IFELSE([AC_LANG_SOURCE([
                void main(void) {
                double a = 1.0, b = 1.0, c = 1.0;
                int one = 1;
                dgemm_("n", "n", &one, &one, &one, &a, &a, &one, &b, &one,
                       &a, &c, &one);
                if (c == 2.0)
                    exit(0);
                else
                    exit(1);
                }
                ])],
               [],
               [AC_MSG_RESULT(no)
                AC_ERROR(Could not link)])
AC_RUN_IFELSE([AC_LANG_SOURCE([
                void main(void) {
                double a = 1.0, b = 1.0, c = 1.0;
                int one = 1;
                dgemm_("n", "n", &one, &one, &one, &a, &a, &one, &b, &one,
                       &a, &c, &one);
                if (c == 2.0)
                    exit(0);
                else
                    exit(1);
                }
               ])],
              [can_run="yes"],
              [can_run="no"])
AC_MSG_RESULT($can_run)

AC_SUBST([CONFIGURE_DEPENDENCIES], ['$(top_srcdir)/version.m4'])

AC_OUTPUT

AS_IF([test "x$can_run" == "xno"],
      [AC_WARN(Could not run compiled binary: check your LD_LIBRARY_PATH)])

cat<< EOF

-------------------------------------------------------------------------------
ester-$VERSION build configuration:

prefix:     $prefix
FC:         $FC
CPPFLAGS:   $CPPFLAGS
CXX:        $CXX
CXXFLAGS:   $CXXFLAGS
LIBS:       $LIBS
HDF5:       $enable_hdf5
swig:       $have_swig
python:     $PYTHON_VERSION
py module:  $pymodule
matplotlib: $enable_plt
-------------------------------------------------------------------------------

EOF
