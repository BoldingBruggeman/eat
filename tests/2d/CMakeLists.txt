# Copyright (C) 2021 Bolding & Bruggeman

add_executable(eat_2d_model eat_2d_model.F90)
target_link_libraries(eat_2d_model PRIVATE eat_config)
target_link_libraries(eat_2d_model PRIVATE ${MPI_Fortran_LIBRARIES})
target_link_libraries(eat_2d_model PRIVATE ${LAPACK_LIBRARIES})
