# Copyright (C) 2021 Bolding & Bruggeman

cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0091 NEW)

project(eat VERSION 0.1.0 LANGUAGES Fortran)

# In principle the code should be fixed - but here we provide a working fix
# Set -fallow-argument-mismatch for major version >= 10
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    string(REPLACE "." ";" VERSION_LIST ${CMAKE_Fortran_COMPILER_VERSION})
    list(GET VERSION_LIST 0 VERSION_MAJOR)
    if(${VERSION_MAJOR} GREATER "9")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
#KB - will fail when compiling C-code        add_definitions(-fallow-argument-mismatch)
    endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/")

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

if(MPI_msmpifec_LIBRARY)
  file(TO_CMAKE_PATH "${MPI_msmpifec_LIBRARY}" MPI_msmpifec_LIBRARY)
endif()

if(EAT_USE_MPI)
find_package(MPI REQUIRED)
if(MPI_Fortran_INCLUDE_PATH)
    include_directories("${MPI_Fortran_INCLUDE_PATH}")
endif()
endif(EAT_USE_MPI)

add_subdirectory("extern/datetime-fortran" EXCLUDE_FROM_ALL)
if(WIN32)
  set_property(TARGET strptime PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# As PDAF can't be added as a submodule we need to find out where it is located
# We look in EAT_PDAF_BASE (as configured via the call to cmake)
if(EAT_USE_PDAF)
   # PDAF naming scheme changed from V2.0 "-D" was dropped
   find_path(PDAF_BASE NAMES src/PDAF_init.F90 src/PDAF-D_init.F90 HINTS ${EAT_PDAF_BASE} extern/pdaf DOC "Path to PDAF source directory." NO_CMAKE_FIND_ROOT_PATH)
   message("${PDAF_BASE}")
   add_subdirectory(${PDAF_BASE} "extern/pdaf" EXCLUDE_FROM_ALL)
endif(EAT_USE_PDAF)

if(EAT_USE_GOTM)
    set(FABM_BASE "${CMAKE_CURRENT_LIST_DIR}/extern/fabm" CACHE PATH "FABM source directory")
    find_path(GOTM_BASE NAMES src/gotm/gotm.F90 HINTS ENV EAT_GOTM_BASE PATHS "extern/gotm" NO_DEFAULT_PATH)
    add_subdirectory("${GOTM_BASE}" extern/gotm EXCLUDE_FROM_ALL)
    message("GOTM folder: " "${GOTM_BASE}")
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
  set(DEFAULT_PREFIX "$ENV{HOME}/.local")
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install EAT in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install EAT in")
endif()

add_subdirectory("src")
add_subdirectory("eatpy")
add_subdirectory("tests")
