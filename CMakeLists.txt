# Copyright (C) 2021 Bolding & Bruggeman

cmake_minimum_required(VERSION 3.0)

if(NOT PROJECT_NAME STREQUAL eat)
   project(eat VERSION 0.1.0 LANGUAGES Fortran)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/")

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
   set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Specify default build type for single-build-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# To check version of configuration file
#KBadd_definitions(-D_CFG_VERSION_=${gotm_VERSION_MAJOR})

# ----------------------------
# User-configurable options
# ----------------------------

option(GOTM_USE_CVMIX "Include support for the CVMix turbulence closure models (github.com/CVMix)" OFF)

option(GOTM_USE_FABM "Include support for Framework for Aquatic Biogeochemical Models (fabm.net)" ON)

option(GOTM_USE_STIM "Include support for Simple Thermodynamic Ice Models - STIM" OFF)

# ----------------------------
# Dependencies
# ----------------------------

add_subdirectory("extern/gotm")

find_package(MPI REQUIRED)
if(MPI_Fortran_INCLUDE_PATH)
    include_directories("${MPI_Fortran_INCLUDE_PATH}")
endif()
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${MPI_Fortran_LIBRARIES})

#KBoption(GOTM_USE_NetCDF "Enable output in NetCDF format" ON)

#KBif(MSVC)
  #KBoption(GOTM_EMBED_VERSION "Embed GOTM version information" OFF)
#KBelse(MSVC)
  #KBoption(GOTM_EMBED_VERSION "Embed GOTM version information" ON)
#KBendif(MSVC)

#KBoption(GOTM_USE_SEAGRASS "Enable seagrass module" OFF)

# ----------------------------
# EAT itself
# ----------------------------

# Set default installation prefix (done after dependencies to ensure ours takes precedence)
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/gotm")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/gotm")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/gotm")
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GOTM in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GOTM in")
endif()

add_subdirectory("src")
add_subdirectory("tests")

#KBinclude(FeatureSummary)
#KBif(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
   #KB#feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES PACKAGES_FOUND)
   #KBfeature_summary(FILENAME ${CMAKE_CURRENT_BINARY_DIR}/features.log WHAT ALL)
#KBendif()
