# Copyright (C) 2021 Bolding & Bruggeman

add_library(eat_config eat_config.F90 ${MPI_Fortran_EXTRA_SOURCES})

if(EAT_USE_PDAF)
add_definitions(-D_PDAF_)
add_executable(eat_pdaf eat_pdaf.F90)
target_link_libraries(eat_pdaf PRIVATE eat_config)
target_link_libraries(eat_pdaf PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

add_executable(eat_model eat_model.F90)
target_link_libraries(eat_model PRIVATE gotm datetime eat_config)
target_link_libraries(eat_model PRIVATE ${MPI_Fortran_LIBRARIES})

add_executable(eat_observations eat_observations.F90)
target_link_libraries(eat_observations PRIVATE eat_config)
target_link_libraries(eat_observations PRIVATE ${MPI_Fortran_LIBRARIES})

add_executable(eat_model_gotm eat_model_gotm.F90)
target_link_libraries(eat_model_gotm PRIVATE gotm datetime eat_config)
target_link_libraries(eat_model_gotm PRIVATE ${MPI_Fortran_LIBRARIES})

add_executable(eat_model_seamless eat_model_seamless.F90)
target_link_libraries(eat_model_seamless PRIVATE gotm datetime eat_config)
target_link_libraries(eat_model_seamless PRIVATE ${MPI_Fortran_LIBRARIES})

