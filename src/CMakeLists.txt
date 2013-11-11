project(CalliopeLib)

include_directories(${Petaca_BINARY_DIR}/mod_files)

add_library(calliope kinds.F90 idaesol_type.F90 nka_type.F90 kinds.F90)
target_link_libraries(calliope petaca)

## # Create all the Fortran module files in a directory at the root build directory
## SET(Calliope_MODULE_DIR ${Calliope_BINARY_DIR}/mod_files/)
## SET(CMAKE_Fortran_MODULE_DIRECTORY ${Calliope_MODULE_DIR})
## 
## # Default is an debug build
## SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type")
## 
## ### Set the default Fortran compiler flags (CMake doesn't do this for Fortran)
## IF(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
##   SET(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -w=uda -DNDEBUG" CACHE STRING
##       "Flags used by the compiler during release builds.")
##   SET(CMAKE_Fortran_FLAGS_DEBUG   "-u -g90 -gline -C=all" CACHE STRING
##       "Flags used by the compiler during debug builds.")
## ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
##   SET(CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -DNDEBUG" CACHE STRING
##       "Flags used by the compiler during release builds.")
##   SET(CMAKE_Fortran_FLAGS_DEBUG   "-u -g" CACHE STRING
##       "Flags used by the compiler during debug builds.")
## ENDIF()
## 
## ### Hardwire the petaca library for now
## IF(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
##   SET(Petaca_INCLUDE_DIR /opt/petaca/nag/include)
##   SET(Petaca_LIBRARY_DIR 
## ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
## ENDIF()
## 
## SET(Calliope_SOURCE_FILES
##     bdf2_dae.F90 f90_assert.F90 kinds.F90 nodal_fem_1d.F90)
## 
## INCLUDE_DIRECTORIES(${Calliope_MODULE_DIR})
