cmake_minimum_required(VERSION 3.22)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

project(calliope Fortran C)

find_package(PETACA REQUIRED)

# Set missing defaults for the NAG Fortran compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set(CMAKE_Fortran_FLAGS_DEBUG "-u -O0 -g -gline -C=all -nan")
  set(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -w=uda")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set(CMAKE_Fortran_LINKER_WRAPPER_FLAG "-Wl,")
  set(CMAKE_Fortran_LINKER_WRAPPER_FLAG_SEP ",")
endif()

# Required compiler options
add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-f2018>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
)

# Disable assertions in release builds
add_compile_definitions($<$<CONFIG:Release,RelWithDebInfo>:NDEBUG>)

# Workarounds for compiler bugs
add_compile_definitions(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:GNU_PR109846>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:INTEL_BUG20180222>
)

# Set the default installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(NOT WIN32 AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(CMAKE_INSTALL_PREFIX "/opt/${PROJECT_NAME}")
  endif()
endif()

set(CMAKE_INSTALL_RPATH $ORIGIN)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CTest)

add_subdirectory(src)
add_subdirectory(test)
