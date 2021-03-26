# Copyright (C) 2021 Bolding & Bruggeman

add_library(mpi_config mpi_config.F90 ${MPI_Fortran_EXTRA_SOURCES})

if(EAT_USE_PDAF)
add_definitions(-D_PDAF_)
add_executable(eat_pdaf eat_pdaf.F90)
target_link_libraries(eat_pdaf PRIVATE mpi_config)
target_link_libraries(eat_pdaf PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

add_executable(eat_server eat_server.F90)
target_link_libraries(eat_server PRIVATE gotm mpi_config)
target_link_libraries(eat_server PRIVATE ${MPI_Fortran_LIBRARIES})

add_executable(eat_worker eat_worker.F90)
target_link_libraries(eat_worker PRIVATE gotm mpi_config)
target_link_libraries(eat_worker PRIVATE ${MPI_Fortran_LIBRARIES})

add_executable(eat_observations eat_observations.F90)
target_link_libraries(eat_observations PRIVATE mpi_config)
target_link_libraries(eat_observations PRIVATE ${MPI_Fortran_LIBRARIES})

