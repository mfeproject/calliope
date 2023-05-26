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

#add_library(common OBJECT pde_class.F90 mfe1_disc_core_type.F90)
#set_target_properties(common PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(common pde_class.F90 cell_data_type.F90 pde_utilities.F90)
target_link_libraries(common petaca)

# LIBRARIES IMPLEMENTING THE PDE CLASS
add_library(navier-stokes navier_stokes_type.F90)
target_link_libraries(navier-stokes common)

add_library(conv-diff conv_diff_type.F90)
target_link_libraries(conv-diff common)

add_library(burgers burgers_type.F90)
target_link_libraries(burgers common)

add_library(drift-diff drift_diff_type.F90)
target_link_libraries(drift-diff common)

set(PDELIBS common navier-stokes conv-diff burgers drift-diff)

install(TARGETS ${PDELIBS})

# User-facing executable is a shell script that adds an argument to the command
# line giving the absolute path to the PDE library directory. Actual executable
# goes into the lib directory.
install(PROGRAMS mfe1.sh
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    RENAME mfe1
)
install(TARGETS mfe1
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
