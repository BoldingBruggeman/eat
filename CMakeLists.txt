# Copyright (C) 2021 Bolding & Bruggeman

cmake_minimum_required(VERSION 3.0)

project(eat VERSION 0.1.0 LANGUAGES Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/")

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Specify default build type for single-build-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# To check version of configuration file
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
include_directories("${CMAKE_BINARY_DIR}/modules")

# ----------------------------
# User-configurable options
# ----------------------------

option(EAT_USE_MPI "Enable MPI parallelization" ON)
set(PDAF_USE_MPI ${EAT_USE_MPI} CACHE BOOL "Compile PDAF with MPI support" FORCE)

option(EAT_USE_PDAF "Include support for the PDAF library" ON)

option(EAT_USE_GOTM "Build the GOTM driver for EAT" ON)

# ----------------------------
# Dependencies
# ----------------------------

if(EAT_USE_MPI)
find_package(MPI REQUIRED)
if(MPI_Fortran_INCLUDE_PATH)
    include_directories("${MPI_Fortran_INCLUDE_PATH}")
endif()
endif(EAT_USE_MPI)

add_subdirectory("extern/datetime-fortran" EXCLUDE_FROM_ALL)

# As PDAF can't be added as a submodule we need to find out where it is located
# First we look in EAT_PDAF_BASE (as configured via the call to cmake) - next we look in .../extern/pdaf
if(EAT_USE_PDAF)
   find_path(PDAF_BASE src/PDAF-D_init.F90 HINTS ${EAT_PDAF_BASE} ${CMAKE_CURRENT_LIST_DIR}/extern/pdaf DOC "Path to PDAF source directory.")
   add_subdirectory(${PDAF_BASE} "extern/pdaf" EXCLUDE_FROM_ALL)
endif(EAT_USE_PDAF)

if(EAT_USE_GOTM)
#KB    find_package(gotm REQUIRED CONFIG HINTS "${GOTM_PREFIX}" NO_DEFAULT_PATH)
#KB    mark_as_advanced(CLEAR GOTM_PREFIX)
#KB    add_subdirectory("models/gotm" EXCLUDE_FROM_ALL)
#KB    set(GOTM_BUILD_LIBRARIES_ONLY ON)
    add_subdirectory("extern/gotm")
    add_subdirectory("models/gotm")
endif(EAT_USE_GOTM)

# ----------------------------
# EAT itself
# ----------------------------

# Set default installation prefix (done after dependencies to ensure ours takes precedence)
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/eat")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/eat")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/eat")
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install EAT in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install EAT in")
endif()

add_subdirectory("src")
add_subdirectory("tests")
