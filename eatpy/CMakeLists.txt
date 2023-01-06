# Copyright (C) 2021 Bolding & Bruggeman
cmake_minimum_required(VERSION 3.0)

# Create main project.
project(ceat_filter_pdaf VERSION 0.1.0 LANGUAGES Fortran C)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Find target Python interpreter.
if(Python3_EXECUTABLE)
  file(TO_CMAKE_PATH "${Python3_EXECUTABLE}" Python3_EXECUTABLE)
endif()
set(Python_FIND_STRATEGY LOCATION)
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development)

add_library(ceat_filter_pdaf ceat_filter_pdaf.F90)
target_link_libraries(ceat_filter_pdaf PRIVATE eat_config)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/_eat_filter_pdaf.c
  COMMAND ${Python3_EXECUTABLE} -m cython _eat_filter_pdaf.pyx
  MAIN_DEPENDENCY _eat_filter_pdaf.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

Python3_add_library(_eat_filter_pdaf MODULE
  _eat_filter_pdaf.c)
target_link_libraries(_eat_filter_pdaf PRIVATE ceat_filter_pdaf)

# Build wheel
add_custom_command(OUTPUT wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}" eatpy
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:_eat_filter_pdaf>" eatpy/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/../setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --force-reinstall --upgrade eatpy WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
