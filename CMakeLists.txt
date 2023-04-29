cmake_minimum_required(VERSION 3.20.2)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

project(Calliope Fortran C)

find_package(PETACA REQUIRED)

### Set the default Fortran compiler flags for NAG (CMake doesn't)
if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set(CMAKE_Fortran_FLAGS_DEBUG "-u -O0 -g -gline -C=all -nan")
  set(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -f2018 -w=uda")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")
endif()

# Disable assertions in release builds
add_compile_definitions($<$<CONFIG:Release,RelWithDebInfo>:NDEBUG>)

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fimplicit-none>
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>)

include(CTest)

add_subdirectory(src)
add_subdirectory(test)
