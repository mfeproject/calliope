set(SOURCE
  block_linear_solver.F90
  btd_matrix_type.F90
  coord_grid_type.F90
  idaesol_type.F90
  index_func_type.F90
  mfe1_disc_type.F90
  mfe1_solver_type.F90
  mfe1_vector_type.F90
  mfe_env_type.F90
  mfe_idaesol_model_type.F90
  mfe_model_type.F90
  mfe_precon_type.F90
  mfe_sim_type.F90
  nka_type.F90
  sim_event_queue_type.F90
  simlog_type.F90
  state_history_type.F90
  time_step_sync_type.F90
  vector_class.F90
  string_utilities.F90
  mfe1.F90
  )

add_executable(mfe1 ${SOURCE})
target_link_libraries(mfe1 common)

set(BUILD_SHARED_LIBS YES)

add_library(common OBJECT pde_class.F90 mfe1_disc_core_type.F90)
set_target_properties(common PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(common petaca)

add_library(navier-stokes navier_stokes_type.F90)
target_link_libraries(navier-stokes common)

add_library(conv-diff conv_diff_type.F90)
target_link_libraries(conv-diff common)

add_library(burgers burgers_type.F90)
target_link_libraries(burgers common)

add_library(drift-diff drift_diff_type.F90)
target_link_libraries(drift-diff common)

## # Create all the Fortran module files in a directory at the root build directory
## SET(Calliope_MODULE_DIR ${Calliope_BINARY_DIR}/mod_files/)
## SET(CMAKE_Fortran_MODULE_DIRECTORY ${Calliope_MODULE_DIR})
##
## # Default is an debug build
## SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type")
##
