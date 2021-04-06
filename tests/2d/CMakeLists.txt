# Copyright (C) 2021 Bolding & Bruggeman

add_executable(eat_model_2d eat_model_2d.F90)
target_link_libraries(eat_model_2d PRIVATE eat_config)
target_link_libraries(eat_model_2d PRIVATE ${MPI_Fortran_LIBRARIES})
