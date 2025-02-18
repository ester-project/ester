cmake_minimum_required(VERSION 3.14)
project(ester VERSION 1.2.0 LANGUAGES C CXX Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON) # test line 1


if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds not allowed.")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    string(CONCAT CMAKE_Fortran_FLAGS
        "-g -O3"
        " -fno-align-commons -fdefault-real-8 -fdefault-double-8"
        " -ffixed-line-length-none -O -w -fd-lines-as-comments")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    string(CONCAT CMAKE_Fortran_FLAGS
        "-g -O3"
        " -r8"
        " -nowarn -fpp -extend-source")
else()
    message(FATAL_ERROR "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

set(CMAKE_CXX_FLAGS "-g")

# Find libarries
find_package(BLAS REQUIRED)
if(BLAS_LIBRARIES MATCHES "mkl")
    message(STATUS "Using mkl BLAS implementation")
    add_compile_definitions(USE_MKL)
endif()
find_package(LAPACK REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(SWIG REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
# find_package(VTK 6.3)

if(NOT DEFINED CBLAS_LIBRARIES)
    message(WARNING "CBLAS_LIBRARIES not specified in command line, trying to find it automatically (this may fail)")
endif()

if(VTK_FOUND)
    add_compile_definitions(USE_VTK)
    include(${VTK_USE_FILE})
    set(OPTIONAL_LIBS "${OPTIONAL_LIBS};${VTK_LIBRARIES}")
else()
    add_compile_definitions(ENABLE_PLT)
endif()

###

# Add the new C++ source file
#add_library(read_hdf5_data STATIC src/physics/read_hdf5_data.cpp)

# Link HDF5 libraries to the new C++ library
#target_link_libraries(read_hdf5_data ${HDF5_LIBRARIES})

add_subdirectory(tables)
add_subdirectory(python)

include_directories(src/include)
include_directories(src/graphics)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

### LIBMESAC INTEGRATION ###

# Include and link directories for MESA
#include_directories($ENV{MESA_DIR}/include $ENV{MESASDK_ROOT}/include)
#link_directories($ENV{MESA_DIR}/lib $ENV{MESASDK_ROOT}/lib)

#include_directories($ENV{LIBMESAC_ROOT}/include)
#link_directories($ENV{LIBMESAC_ROOT}/lib) # does this need to be added soemhow? 

# Run the MESA shell commands and capture their output
execute_process(
    COMMAND mesasdk_lapack95_link
    OUTPUT_VARIABLE LAPACK95_LIB_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND mesasdk_crmath_link
    OUTPUT_VARIABLE CRMATH_LIB_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Additional MESA libraries

#set(MESA_LIBS "-L$ENV{MESA_DIR}/lib -lkap -lchem -lmath -latm -lconst -lutils -leos -linterp_1d -linterp_2d -lnum -lionization -lmtx -lauto_diff -lcrmath -lhdf5io -lrates -lnet -lmvec -fopenmp -lm -lquadmath")

# Combine the outputs into one variable for SDK libraries
#set(SDK_LIBS "${MESA_LIBS};${LAPACK95_LIB_FLAGS};${CRMATH_LIB_FLAGS}")

#message(LIBMESAC_ROOT="$ENV{LIBMESAC_ROOT}")

#set(LIBMESAC_LIBS "$ENV{LIBMESAC_ROOT}/lib/libkapcmesa.a;$ENV{LIBMESAC_ROOT}/lib/libchemcmesa.a;$ENV{LIBMESAC_ROOT}/lib/libutilscmesa.a;${SDK_LIBS}")

### -------------------- ### 

# Find HDF5 with C and C++ components
#find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Include HDF5 headers
#include_directories(${HDF5_INCLUDE_DIRS})

# Link against HDF5 libraries
#target_link_libraries(ester PRIVATE ${HDF5_LIBRARIES})

###

if(DEFINED PYTHON_NUMPY_INCLUDE_DIR)
    include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
else()
    message(WARNING "PYTHON_NUMPY_INCLUDE_DIR not set, this may fail if numpy has been installed locally, using pip for example")
endif()

add_compile_definitions(WITH_CMAKE)
add_compile_definitions(VERSION="${PROJECT_VERSION}")
add_compile_definitions(ESTER_DATADIR="${CMAKE_INSTALL_PREFIX}/share")
add_compile_definitions(PACKAGE_VERSION="${PROJECT_VERSION}")

add_library(ester SHARED
    src/global/global.cpp
    src/graphics/matplotlib.cpp
    src/solver/RKF_solver.cpp
    src/main/read_config.cpp
    src/solver/SDIRK_solver.cpp
    src/solver/solver.cpp
    src/solver/solver_full.cpp
    src/star/convection.cpp
    src/star/polytrope.cpp
    src/star/star1d_class.cpp
    src/mapping/mapping.cpp
    src/star/star1d_extra.cpp
    src/mapping/remapper.cpp
    src/star/star1d_solvers.cpp
    src/matrix/mat_math.cpp
    src/star/star2d_class.cpp
    src/matrix/matrix_block_diag.cpp
    src/star/star2d_extra.cpp
    src/matrix/matrix.cpp
    src/star/star2d_solvers.cpp
    src/matrix/matrix_map.cpp
    src/star/star_evol.cpp
    src/matrix/mat_spec.cpp
    src/star/star_map.cpp
    src/numdiff/diff_gl.cpp
    src/star/star_phys.cpp
    src/numdiff/diff_leg.cpp
    src/symbolic/rational.cpp
    src/symbolic/symbolic.cpp
    src/symbolic/sym.cpp
    src/symbolic/sym_expr.cpp
    src/symbolic/sym_tens.cpp
    src/symbolic/sym_vec.cpp
    src/utils/debug.cpp
    src/utils/iofile.cpp
    src/utils/parser.cpp
    src/utils/stack.cpp
    src/physics/opa_kramer.cpp
    src/physics/opa_opal.cpp
    src/physics/physics.cpp
    src/physics/atm_onelayer.cpp
    src/physics/composition.cpp
    src/physics/eos_freeeos.cpp
    src/physics/eos_ideal.cpp
    src/physics/eos_idealrad.cpp
    src/physics/eos_opal.cpp
    src/physics/nuc_cesam.cpp
    src/physics/nuc_simple.cpp
    src/physics/opa_cesam.cpp
    src/physics/opa_houdek.cpp
    src/physics/opa_cesam_for.F90
    src/physics/nuc_cesam_for.F
    src/physics/EOS5_xtrin.F
    src/physics/xztrin21.F
    src/physics/ZFS_interp_EOS5.F
    src/physics/opa_opmesa.F90
    src/physics/opa_opmesa.cpp
    src/physics/mesa_opacity_modules.F90
    src/physics/cubic_interpolator.f90
    src/physics/read_hdf5_data.cpp)
    #src/physics/opa_mesa_wrapper.cpp) # blocked for fixing other issues 
    
# Link libmesac libraries using --whole-archive to include all symbols
#set(LIBMESAC_LIBS "-Wl,--whole-archive $ENV{LIBMESAC_ROOT}/lib/libkapcmesa.a $ENV{LIBMESAC_ROOT}/lib/libchemcmesa.a $ENV{LIBMESAC_ROOT}/lib/libutilscmesa.a -#Wl,--no-whole-archive ${SDK_LIBS}") # test lines 2 

   
message(LIBMESAC_LIBS="${LIBMESAC_LIBS}")
target_link_libraries(ester cesam freeeos opint
    ${OPTIONAL_LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${CBLAS_LIBRARIES} ${PYTHON_LIBRARIES} ${HDF5_LIBRARIES} ${LIBMESAC_LIBS})

# Set properties for ester
#set_target_properties(ester PROPERTIES
#    LINK_SEARCH_START_STATIC ON
#    LINK_SEARCH_END_STATIC OFF
#) # test lines 3 
 

add_custom_target(ester_bin ALL
    COMMAND sed -e 's,[@]PYTHON[@],${PYTHON_EXECUTABLE},g' -e
    's,[@]prefix[@],${CMAKE_INSTALL_PREFIX},g' -e
    's,[@]datadir[@],${CMAKE_INSTALL_PREFIX}/share/ester,g' <
    ${CMAKE_SOURCE_DIR}/src/main/ester.in > ester)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ester
    DESTINATION bin)

install(FILES
    config/star.cfg
    config/1d_default.par
    config/2d_default.par
    DESTINATION share/ester/)

add_executable(polytrope1d src/main/polytrope1d.cpp)
target_link_libraries(polytrope1d ester)

add_executable(polytrope2d src/main/polytrope2d.cpp)
target_link_libraries(polytrope2d ester)

add_executable(star1d src/main/star1d.cpp)
target_link_libraries(star1d ester)

add_executable(star2d src/main/star2d.cpp)
target_link_libraries(star2d ester)

add_executable(star_evol src/main/star_evol.cpp)
target_link_libraries(star_evol ester)

add_executable(version src/main/version.cpp)
target_link_libraries(version ester)

add_executable(gen_output src/main/gen_output.cpp)
target_link_libraries(gen_output ester)

add_executable(ester_info src/main/ester-info.cpp)
target_link_libraries(ester_info ester)

add_executable(ester_vtk src/main/vtk.cpp)
target_link_libraries(ester_vtk ester)

install(TARGETS ester polytrope1d polytrope2d star1d star2d star_evol version gen_output ester_info ester_vtk
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)
