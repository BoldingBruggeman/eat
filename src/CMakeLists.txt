# Copyright (C) 2021 Bolding & Bruggeman

add_library(eat_config eat_config.F90 ${MPI_Fortran_EXTRA_SOURCES})

add_executable(eat_pdaf eat_pdaf.F90)
target_link_libraries(eat_pdaf PRIVATE eat_config pdaf-d ${MPI_Fortran_LIBRARIES} )

add_executable(eat_model_gotm eat_model_gotm.F90)
target_link_libraries(eat_model_gotm PRIVATE gotm datetime eat_config ${MPI_Fortran_LIBRARIES})
