# Copyright (C) 2021 Bolding & Bruggeman

add_library(eat_config eat_config.F90 pdaf_wrapper.F90 ${MPI_Fortran_ADDITIONAL_SOURCES})
target_link_libraries(eat_config PUBLIC pdaf-d ${MPI_Fortran_LIBRARIES} )

add_executable(eat_filter_pdaf eat_filter_pdaf.F90)
target_link_libraries(eat_filter_pdaf PUBLIC eat_config ${MPI_Fortran_LIBRARIES} )

install(TARGETS eat_config
   EXPORT eatConfig
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
)
install(FILES
   $<TARGET_PROPERTY:eat_config,Fortran_MODULE_DIRECTORY>/\${BUILD_TYPE}/eat_config.mod
   $<TARGET_PROPERTY:eat_config,Fortran_MODULE_DIRECTORY>/\${BUILD_TYPE}/pdaf_wrapper.mod
   DESTINATION include
)
install(TARGETS eat_filter_pdaf RUNTIME DESTINATION bin)

#export(EXPORT eatConfig)
#install(EXPORT eatConfig DESTINATION cmake NAMESPACE eat::)
