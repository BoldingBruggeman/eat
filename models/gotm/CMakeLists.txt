add_executable(eat-gotm eat_model_gotm.F90)
target_link_libraries(eat-gotm PRIVATE gotm datetime eat_config ${MPI_Fortran_LIBRARIES})

add_executable(eat_obs_gotm eat_obs_gotm.F90)
target_link_libraries(eat_obs_gotm PRIVATE eat_config ${MPI_Fortran_LIBRARIES})

install(TARGETS eat-gotm eat_obs_gotm
   EXPORT gotmConfig
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
)

