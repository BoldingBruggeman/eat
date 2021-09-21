# Copyright (C) 2021 Bolding & Bruggeman

add_library(eat_config eat_config.F90 ${MPI_Fortran_EXTRA_SOURCES})

add_executable(eat_filter_pdaf eat_filter_pdaf.F90)
target_link_libraries(eat_filter_pdaf PUBLIC eat_config pdaf-d ${MPI_Fortran_LIBRARIES} )

install(TARGETS eat_config
   EXPORT eatConfig
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
)
install(FILES
   $<TARGET_PROPERTY:eat_config,Fortran_MODULE_DIRECTORY>/eat_config.mod
   DESTINATION include
)
install(TARGETS eat_filter_pdaf RUNTIME DESTINATION bin)

export(EXPORT eatConfig)
install(EXPORT eatConfig DESTINATION cmake NAMESPACE eat::)
