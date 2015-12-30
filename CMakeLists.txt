cmake_minimum_required(VERSION 2.8.12)

project(Calliope Fortran)

## # Create all the Fortran module files in a directory at the root build directory
## SET(Calliope_MODULE_DIR ${Calliope_BINARY_DIR}/mod_files/)
## SET(CMAKE_Fortran_MODULE_DIRECTORY ${Calliope_MODULE_DIR})

# Default is an debug build
SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type")

### Set the default Fortran compiler flags (CMake doesn't do this for Fortran)
IF(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  SET(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -w=uda -DNDEBUG")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-u -g90 -gline -C=all -nan")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -DNDEBUG")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-u -g")
ENDIF()

add_subdirectory(petaca)
add_subdirectory(src)
add_subdirectory(test)
