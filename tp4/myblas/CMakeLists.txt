###
#
#  @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @author Mathieu Faverge
#  @date 2019-11-13
#
###
cmake_minimum_required (VERSION 2.8)

# Disable in source directory build
set( CMAKE_DISABLE_IN_SOURCE_BUILD ON )

# use, i.e. don't skip the full RPATH for the build tree
set( CMAKE_SKIP_BUILD_RPATH FALSE )

# when building, don't use the install RPATH already
# (but later on when installing)
set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )

option(BUILD_SHARED_LIBS
  "Build shared libraries" ON)

if ( NOT DEFINED ALGONUM_DIR )
  set( ALGONUM_DIR ${CMAKE_SOURCE_DIR}/../install/ )
endif()

include_directories( ${ALGONUM_DIR}/include )
link_directories( ${ALGONUM_DIR}/lib )

# Define the list of source files for the myblas library
set( myblas_srcs
  mygemm.c
  mygetrf.c
)

if( ENABLE_STARPU )
  include_directories( ${ALGONUM_DIR}/include )
  list( APPEND myblas_srcs
    codelet_dgemm.c
    codelet_dgetrf.c
    codelet_dtrsm.c
    starpu_dgemm.c
    starpu_dgetrf.c
    )
endif()

# Add the files of your own myblas library here
add_library( myblas
  ${myblas_srcs}
  )

if( ENABLE_STARPU )
  # Configuration with MKL
  set( LAPACKE_LIBRARY_DIRS_DEP
    $ENV{MKLROOT}/lib/intel64
    )
  set( LAPACKE_LIBRARIES_DEP
    mkl_gf_lp64
    mkl_sequential
    mkl_core
    )

  target_link_libraries( myblas
    algonum
    ${STARPU_LIBRARIES}
    ${LAPACKE_LIBRARIES_DEP}
    )
endif()

install(TARGETS myblas
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

