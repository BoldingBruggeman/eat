#find_package(MPI REQUIRED)
#if(MPI_Fortran_INCLUDE_PATH)
    #include_directories("${MPI_Fortran_INCLUDE_PATH}")
#endif()
add_executable(eat_model_gotm eat_model_gotm.F90)
target_link_libraries(eat_model_gotm PRIVATE gotm datetime eat_config ${MPI_Fortran_LIBRARIES})

install(TARGETS eat_model_gotm
   EXPORT gotmConfig
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
)

