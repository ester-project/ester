add_definitions(-DWITH_CMAKE)
add_definitions(-DENABLE_PLT)
include_directories(${HDF5_INCLUDE_DIRS})

add_custom_target(__init__.py ALL
    COMMAND sed -e 's,[@]PYTHON[@],${PYTHON_EXECUTABLE},g' ${CMAKE_SOURCE_DIR}/python/ester.in > __init__.py)

add_custom_target(cesam.py ALL
    COMMAND sed -e 's,[@]PYTHON[@],${PYTHON_EXECUTABLE},g' ${CMAKE_SOURCE_DIR}/python/cesam.in > cesam.py)

execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(prefix='', plat_specific=True))"
    OUTPUT_VARIABLE PYTHON_SITE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

include(${SWIG_USE_FILE})

set(CMAKE_SWIG_FLAGS "-I${CMAKE_SOURCE_DIR}/src/include;-I${CMAKE_SOURCE_DIR}/src/graphics;")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/src/include)
include_directories(${CMAKE_SOURCE_DIR}/src/graphics)

set_property(SOURCE ester.i PROPERTY CPLUSPLUS ON)
if(${CMAKE_VERSION} VERSION_LESS "3.8")
    swig_add_module(ester_wrap python ester.i)
else()
    swig_add_library(ester_wrap LANGUAGE python SOURCES ester.i)
endif()
target_link_libraries(${SWIG_MODULE_ester_wrap_REAL_NAME} ester cesam ${HDF5_CXX_LIBRARIES})

install(TARGETS ${SWIG_MODULE_ester_wrap_REAL_NAME}
    DESTINATION ${PYTHON_SITE}/ester)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
    ${CMAKE_CURRENT_BINARY_DIR}/ester_wrap.py
    ${CMAKE_CURRENT_BINARY_DIR}/cesam.py
    DESTINATION ${PYTHON_SITE}/ester)
