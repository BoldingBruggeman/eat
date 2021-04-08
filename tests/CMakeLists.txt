# Copyright (C) 2021 Bolding & Bruggeman

add_executable(eat_seamless_filter EXCLUDE_FROM_ALL seamless/eat_seamless_filter.F90)
target_link_libraries(eat_seamless_filter PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})
add_executable(eat_seamless_model EXCLUDE_FROM_ALL seamless/eat_seamless_model.F90)
target_link_libraries(eat_seamless_model PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})
add_executable(eat_seamless_obs EXCLUDE_FROM_ALL seamless/eat_seamless_obs.F90)
target_link_libraries(eat_seamless_obs PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})
add_custom_target(eat_seamless_test DEPENDS eat_seamless_filter DEPENDS eat_seamless_model DEPENDS eat_seamless_obs)
# mpiexec -np 1 eat_seamles_filter : -np 1 eat_seamless_obs :  -np 6 eat_seamless_model

add_executable(eat_2d_model EXCLUDE_FROM_ALL 2d/eat_2d_model.F90)
target_link_libraries(eat_2d_model PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})
add_executable(eat_2d_obs EXCLUDE_FROM_ALL 2d/eat_2d_obs.F90)
target_link_libraries(eat_2d_obs PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})
add_custom_target(eat_2d_test DEPENDS eat_2d_model DEPENDS eat_2d_obs)
# mpiexec -np 1 eat_pdaf : -np 1 eat_2d_obs :  -np 6 eat_2d_model

#add_custom_target(eat_all_tests eat_2d_test eat_seamless_test)
