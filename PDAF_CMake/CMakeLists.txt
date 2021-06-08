cmake_minimum_required(VERSION 3.12)

#KBset(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/Modules/")

# Create main project.
project(pdaf VERSION 1.16.0 LANGUAGES Fortran)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Not used now - KB
if(USE_GIT_VERSION_INFO)
   set(GIT_COMMIT_ID ${pdaf_VERSION})
   set(GIT_BRANCH_NAME unknown)
   configure_file("${PROJECT_SOURCE_DIR}/util/compilation.F90.in" compilation.F90)
   configure_file("${PROJECT_SOURCE_DIR}/util/pdaf_version.F90.in" pdaf_version.F90)
   option(PDAF_EMBED_VERSION "Embed PDAF version information" OFF)
   if(PDAF_EMBED_VERSION)
      add_custom_target(version
         ${CMAKE_COMMAND} -DINFILE=${PROJECT_SOURCE_DIR}/util/pdaf_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/pdaf_version.F90 -P "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Retrieving description of last PDAF commit..."
        SOURCES ${PROJECT_SOURCE_DIR}/util/pdaf_version.F90.in "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
        VERBATIM
      )
      set_property(TARGET version PROPERTY FOLDER pdaf)
   endif()
endif(USE_GIT_VERSION_INFO)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
endif()

#[[-
configure_file (
  "${PROJECT_SOURCE_DIR}/pdaf_configure.h.in"
  "${PROJECT_BINARY_DIR}/pdaf_configure.h"
  )
-]]

# Customize compiler flags
#[[--
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  else()
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -mkl")
  endif()
endif()

--]]

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#KBfind_package(FFTW3)
#set(LAPACK_LIBRARIES -llapack)
#set(BLAS_LIBRARIES -lblas)
#set(FFT_LIBRARIES -lfftw3)
#KBset(MATH_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${FFT_LIBRARIES})

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

add_subdirectory(src)
#KBadd_subdirectory(testsuite/src/dummymodel_1D)
#KBadd_subdirectory(test)
