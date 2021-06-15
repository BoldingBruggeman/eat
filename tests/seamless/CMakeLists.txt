# Copyright (C) 2021 Bolding & Bruggeman

add_executable(eat_seamless_filter EXCLUDE_FROM_ALL eat_seamless_filter.F90)
target_link_libraries(eat_seamless_filter PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})

add_executable(eat_seamless_model EXCLUDE_FROM_ALL eat_seamless_model.F90)
target_link_libraries(eat_seamless_model PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})

add_executable(eat_seamless_obs EXCLUDE_FROM_ALL eat_seamless_obs.F90)
target_link_libraries(eat_seamless_obs PRIVATE eat_config datetime ${MPI_Fortran_LIBRARIES})
