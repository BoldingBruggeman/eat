# Copyright (C) 2021 Bolding & Bruggeman

add_library(eat_config eat_config.F90 ${MPI_Fortran_EXTRA_SOURCES})

add_executable(eat_filter_pdaf eat_filter_pdaf.F90)
target_link_libraries(eat_filter_pdaf PRIVATE eat_config pdaf-d ${MPI_Fortran_LIBRARIES} )

add_executable(eat_obs eat_obs.F90)
target_link_libraries(eat_obs PRIVATE eat_config ${MPI_Fortran_LIBRARIES})
