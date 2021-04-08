# Copyright (C) 2021 Bolding & Bruggeman

add_library(eat_config eat_config.F90 ${MPI_Fortran_EXTRA_SOURCES})

add_executable(eat_pdaf_filter eat_pdaf_filter.F90)
target_link_libraries(eat_pdaf_filter PRIVATE eat_config pdaf-d ${MPI_Fortran_LIBRARIES} )

add_executable(eat_gotm_model eat_gotm_model.F90)
target_link_libraries(eat_gotm_model PRIVATE gotm datetime eat_config ${MPI_Fortran_LIBRARIES})
