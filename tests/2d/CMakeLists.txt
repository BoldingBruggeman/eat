# Copyright (C) 2021 Bolding & Bruggeman

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   enable_language(C)
   set(BLA_VENDOR Intel10_64lp_seq)
endif()
find_package(LAPACK REQUIRED) # To get - dlarnv - and to keep configuration private in PDAF

add_library(eat_2d EXCLUDE_FROM_ALL eat_2d_data.F90)
target_link_libraries(eat_2d PRIVATE ${LAPACK_LIBRARIES})

add_executable(eat_2d_model EXCLUDE_FROM_ALL eat_2d_model.F90)
target_link_libraries(eat_2d_model PRIVATE eat_config datetime eat_2d ${MPI_Fortran_LIBRARIES})

add_executable(eat_2d_obs EXCLUDE_FROM_ALL eat_2d_obs.F90)
target_link_libraries(eat_2d_obs PRIVATE eat_config datetime eat_2d ${MPI_Fortran_LIBRARIES})
