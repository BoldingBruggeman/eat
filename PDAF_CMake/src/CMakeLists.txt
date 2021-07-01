set(MOD_PDAF
    PDAF_timer.F90
    PDAF_memcount.F90
    PDAF-D_mod_filtermpi.F90
    PDAF-D_mod_filter.F90
   )

set(MOD_INTERFACE
    PDAF_interfaces_module.F90
   )

set(SRC_PDAF_GEN
    PDAF-D_init.F90
    PDAF-D_init_si.F90
    PDAF-D_init_filters.F90
    PDAF-D_alloc_filters.F90
    PDAF-D_print_info.F90
    PDAF-D_print_version.F90
    PDAF_options_filters.F90
    PDAF-D_get_state.F90
    PDAF-D_get_state_si.F90
    PDAF-D_incremental.F90
    PDAF-D_incremental_si.F90
    PDAF-D_set_forget.F90
    PDAF-D_set_forget_local.F90
    PDAF-D_add_increment.F90
    PDAF_generate_rndmat.F90
    PDAF_local_weights.F90
    PDAF_local_weight.F90
    PDAF-D_gather_ens.F90
    PDAF-D_get_memberid.F90
    PDAF-D_get_obsmemberid.F90
    PDAF-D_smoother_shift.F90
    PDAF-D_smoother.F90
    PDAF-D_smoother_local.F90
    PDAF-D_set_smootherens.F90
    PDAF-D_get_smootherens.F90
    PDAF_set_ens_pointer.F90
    PDAF-D_put_state_prepost.F90
    PDAF-D_put_state_prepost_si.F90
    PDAF-D_prepost.F90
    PDAF-D_prepost_si.F90
    PDAF_inflate_ens.F90
    PDAF_sampleens.F90
    PDAF_mvnormalize.F90
    PDAF_eofcovar.F90
    PDAF_diag_histogram.F90
    PDAF_diag_ensstats.F90
    PDAF_diag_effsample.F90
    PDAF-D_gather_dim_obs_f.F90
    PDAF-D_gather_obs_f.F90
    PDAF-D_gather_obs_f2.F90
    PDAF-D_gather_obs_f_flex.F90
    PDAF-D_gather_obs_f2_flex.F90
    PDAF-D_allreduce.F90
    PDAF-D_deallocate.F90
    PDAF-D_get_assim_flag.F90
    PDAF_get_localfilter.F90
    PDAF_inflate_weights.F90
    PDAFomi_put_state_global.F90
    PDAFomi_put_state_local.F90
    PDAFomi_assimilate_global.F90
    PDAFomi_assimilate_local.F90
    PDAF_force_analysis.F90
   )

option(PDAF_USE_SEIK "Include SEIK filter?" ON)
if(PDAF_USE_SEIK)
   add_definitions(-D_PDAF_SEIK_)
   set(SRC_SEIK
       PDAF-D_seik_init.F90
       PDAF-D_seik_alloc.F90
       PDAF_seik_options.F90
       PDAF-D_seik_memtime.F90
       PDAF-D_put_state_seik.F90
       PDAF-D_put_state_seik_si.F90
       PDAF-D_assimilate_seik.F90
       PDAF-D_assimilate_seik_si.F90
       PDAF-D_seik_update.F90
       PDAF-D_seik_analysis.F90
       PDAF-D_seik_resample.F90
       PDAF-D_seik_analysis_newT.F90
       PDAF-D_seik_resample_newT.F90
       PDAF-D_seik_analysis_trans.F90
       PDAF-D_seik_matrixT.F90
       PDAF-D_seik_uinv.F90
       PDAF_seik_omega.F90
       PDAF_seik_TtimesA.F90
      )
else()
   set(SRC_SEIK)
endif(PDAF_USE_SEIK)

option(PDAF_USE_LOCAL_SEIK "Include Local SEIK filter?" ON)
if(PDAF_USE_LOCAL_SEIK)
   add_definitions(-D_PDAF_LOCAL_SEIK_)
   set(SRC_LOCAL_SEIK
       PDAF-D_lseik_init.F90
       PDAF-D_lseik_alloc.F90
       PDAF_lseik_options.F90
       PDAF-D_lseik_memtime.F90
       PDAF-D_put_state_lseik.F90
       PDAF-D_put_state_lseik_si.F90
       PDAF-D_assimilate_lseik.F90
       PDAF-D_assimilate_lseik_si.F90
       PDAF-D_lseik_update.F90
       PDAF-D_lseik_analysis.F90
       PDAF-D_lseik_resample.F90
       PDAF-D_lseik_analysis_trans.F90
     )
else()
   set(SRC_LOCAL_SEIK)
endif(PDAF_USE_LOCAL_SEIK)

option(PDAF_USE_SEEK "Include SEEK filter?" ON)
if(PDAF_USE_SEEK)
   add_definitions(-D_PDAF_SEEK_)
   set(SRC_SEEK
       PDAF-D_seek_init.F90
       PDAF-D_seek_alloc.F90
       PDAF_seek_options.F90
       PDAF-D_seek_memtime.F90
       PDAF-D_put_state_seek.F90
       PDAF-D_put_state_seek_si.F90
       PDAF-D_assimilate_seek.F90
       PDAF-D_assimilate_seek_si.F90
       PDAF-D_seek_update.F90
       PDAF-D_seek_analysis.F90
       PDAF-D_seek_rediag.F90
      )
else()
   set(SRC_SEEK)
endif(PDAF_USE_SEEK)

option(PDAF_USE_ENKF "Include EnKF filter?" ON)
if(PDAF_USE_ENKF)
   add_definitions(-D_PDAF_ENKF_)
   set(SRC_ENKF
       PDAF-D_enkf_init.F90
       PDAF-D_enkf_alloc.F90
       PDAF_enkf_options.F90
       PDAF-D_enkf_memtime.F90
       PDAF-D_put_state_enkf.F90
       PDAF-D_put_state_enkf_si.F90
       PDAF-D_assimilate_enkf.F90
       PDAF-D_assimilate_enkf_si.F90
       PDAF-D_enkf_update.F90
       PDAF-D_enkf_obs_ensemble.F90
       PDAF-D_enkf_gather_resid.F90
       PDAF-D_enkf_analysis_rlm.F90
       PDAF-D_enkf_analysis_rsm.F90
       PDAF_enkf_omega.F90
       PDAF_enkf_Tleft.F90
       PDAF-D_smoother_enkf.F90
     )
else()
   set(SRC_ENKF)
endif(PDAF_USE_ENKF)

option(PDAF_USE_ETKF "Include ETKF filter?" ON)
if(PDAF_USE_ETKF)
   add_definitions(-D_PDAF_ETKF_)
   set(SRC_ETKF
       PDAF-D_etkf_init.F90
       PDAF-D_etkf_alloc.F90
       PDAF_etkf_options.F90
       PDAF-D_etkf_memtime.F90
       PDAF-D_put_state_etkf.F90
       PDAF-D_put_state_etkf_si.F90
       PDAF-D_assimilate_etkf.F90
       PDAF-D_assimilate_etkf_si.F90
       PDAF-D_etkf_update.F90
       PDAF-D_etkf_analysis.F90
       PDAF-D_etkf_analysis_T.F90
       PDAF-D_etkf_analysis_fixed.F90
       PDAF-D_etkf_Tright.F90
       PDAF_etkf_Tleft.F90
     )
else()
   set(SRC_ETKF)
endif(PDAF_USE_ETKF)

option(PDAF_USE_LETKF "Include LETKF filter?" ON)
if(PDAF_USE_LETKF)
   add_definitions(-D_PDAF_LETKF_)
   set(SRC_LETKF
       PDAF-D_letkf_init.F90
       PDAF-D_letkf_alloc.F90
       PDAF_letkf_options.F90
       PDAF-D_letkf_memtime.F90
       PDAF-D_put_state_letkf.F90
       PDAF-D_put_state_letkf_si.F90
       PDAF-D_assimilate_letkf.F90
       PDAF-D_assimilate_letkf_si.F90
       PDAF-D_letkf_update.F90
       PDAF-D_letkf_analysis.F90
       PDAF-D_letkf_analysis_T.F90
       PDAF-D_letkf_analysis_fixed.F90
     )
else()
   set(SRC_LETKF)
endif(PDAF_USE_LETKF)

option(PDAF_USE_ESTKF "Include ESTKF filter?" ON)
if(PDAF_USE_ESTKF)
   add_definitions(-D_PDAF_ESTKF_)
   set(SRC_ESTKF
       PDAF-D_estkf_init.F90
       PDAF-D_estkf_alloc.F90
       PDAF_estkf_options.F90
       PDAF-D_estkf_memtime.F90
       PDAF-D_put_state_estkf.F90
       PDAF-D_put_state_estkf_si.F90
       PDAF-D_assimilate_estkf.F90
       PDAF-D_assimilate_estkf_si.F90
       PDAF-D_estkf_update.F90
       PDAF-D_estkf_analysis.F90
       PDAF-D_estkf_analysis_fixed.F90
       PDAF-D_estkf_AOmega.F90
       PDAF-D_estkf_OmegaA.F90
     )
else()
   set(SRC_ESTKF)
endif(PDAF_USE_ESTKF)

option(PDAF_USE_LESTKF "Include LESTKF filter?" ON)
if(PDAF_USE_LESTKF)
   add_definitions(-D_PDAF_LESTKF_)
   set(SRC_LESTKF
       PDAF-D_lestkf_init.F90
       PDAF-D_lestkf_alloc.F90
       PDAF_lestkf_options.F90
       PDAF-D_lestkf_memtime.F90
       PDAF-D_put_state_lestkf.F90
       PDAF-D_put_state_lestkf_si.F90
       PDAF-D_assimilate_lestkf.F90
       PDAF-D_assimilate_lestkf_si.F90
       PDAF-D_lestkf_update.F90
       PDAF-D_lestkf_analysis.F90
       PDAF-D_lestkf_analysis_fixed.F90
     )
else()
   set(SRC_LESTKF)
endif(PDAF_USE_LESTKF)

option(PDAF_USE_LENKF "Include LENKF filter?" ON)
if(PDAF_USE_LENKF)
   add_definitions(-D_PDAF_LENKF_)
   set(SRC_LENKF
       PDAF-D_lenkf_init.F90
       PDAF-D_lenkf_alloc.F90
       PDAF_lenkf_options.F90
       PDAF-D_lenkf_memtime.F90
       PDAF-D_put_state_lenkf.F90
       PDAF-D_put_state_lenkf_si.F90
       PDAFomi_put_state_lenkf.F90
       PDAF-D_assimilate_lenkf.F90
       PDAF-D_assimilate_lenkf_si.F90
       PDAFomi_assimilate_lenkf.F90
       PDAF-D_lenkf_update.F90
       PDAF-D_enkf_obs_ensemble.F90
       PDAF-D_enkf_gather_resid.F90
       PDAF-D_lenkf_analysis_rsm.F90
       PDAF_enkf_omega.F90
       PDAF_enkf_Tleft.F90
     )
else()
   set(SRC_LENKF)
endif(PDAF_USE_LENKF)

option(PDAF_USE_NETF "Include NETF filter?" ON)
if(PDAF_USE_NETF)
   add_definitions(-D_PDAF_NETF_)
   set(SRC_NETF
       PDAF-D_netf_init.F90
       PDAF-D_netf_alloc.F90
       PDAF_netf_options.F90
       PDAF-D_netf_memtime.F90
       PDAF-D_put_state_netf.F90
       PDAF-D_put_state_netf_si.F90
       PDAF-D_assimilate_netf.F90
       PDAF-D_assimilate_netf_si.F90
       PDAF-D_netf_update.F90
       PDAF-D_netf_smootherT.F90
       PDAF-D_netf_analysis.F90
       PDAF-D_smoother_netf.F90
      )
else()
   set(SRC_NETF)
endif(PDAF_USE_NETF)

option(PDAF_USE_LNETF "Include LNETF filter?" ON)
if(PDAF_USE_LNETF)
   add_definitions(-D_PDAF_LNETF_)
   set(SRC_LNETF
       PDAF-D_lnetf_init.F90
       PDAF-D_lnetf_alloc.F90
       PDAF_lnetf_options.F90
       PDAF-D_lnetf_memtime.F90
       PDAF-D_put_state_lnetf.F90
       PDAF-D_put_state_lnetf_si.F90
       PDAF-D_assimilate_lnetf.F90
       PDAF-D_assimilate_lnetf_si.F90
       PDAF-D_lnetf_update.F90
       PDAF-D_lnetf_analysis.F90
       PDAF-D_lnetf_smootherT.F90
       PDAF-D_smoother_lnetf.F90
     )
else()
   set(SRC_LNETF)
endif(PDAF_USE_LNETF)

# Specific PDAF-routines for
set(SRC_PF
    PDAF-D_pf_init.F90
    PDAF-D_pf_alloc.F90
    PDAF_pf_options.F90
    PDAF-D_pf_memtime.F90
    PDAF-D_put_state_pf.F90
    PDAF-D_put_state_pf_si.F90
    PDAF-D_assimilate_pf.F90
    PDAF-D_assimilate_pf_si.F90
    PDAF-D_pf_update.F90
    PDAF-D_pf_analysis.F90
    PDAF-D_pf_resampling.F90
    PDAF-D_pf_add_noise.F90
   )

# Specific PDAF-routines for generating observations
set(SRC_OBSGEN
    PDAF-D_genobs_init.F90
    PDAF-D_genobs_alloc.F90
    PDAF_genobs_options.F90
    PDAF-D_put_state_generate_obs.F90
    PDAF-D_put_state_generate_obs_si.F90
    PDAFomi_put_state_generate_obs.F90
    PDAF-D_generate_obs.F90
    PDAF-D_generate_obs_si.F90
    PDAFomi_generate_obs.F90
    PDAF-D_gen_obs.F90
   )

# Routines for PDAF-OMI
set(SRC_PDAFOMI
    PDAFomi_obs_f.F90
    PDAFomi_obs_l.F90
    PDAFomi_obs_op.F90
    PDAFomi.F90
    PDAFomi_callback.F90
   )

set(SRC_PDAF
    ${SRC_PDAF_GEN}
    ${SRC_SEIK}
    ${SRC_LOCAL_SEIK}
    ${SRC_SEEK}
    ${SRC_ENKF}
    ${SRC_ETKF}
    ${SRC_LETKF}
    ${SRC_ESTKF}
    ${SRC_LESTKF}
    ${SRC_LENKF}
    ${SRC_NETF}
    ${SRC_LNETF}
    ${SRC_PF}
    ${SRC_OBSGEN}
    ${SRC_PDAFOMI}
   )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
add_library(pdaf-d)
target_sources(pdaf-d PRIVATE
            ${MOD_PDAF}
            ${MOD_INTERFACE}
            ${SRC_PDAF}
           )
target_include_directories(pdaf-d
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)

# In principle the code should be fixed - but here we provide a working fix
# Set -fallow-argument-mismatch for major version >= 10
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
string(REPLACE "." ";" VERSION_LIST ${CMAKE_Fortran_COMPILER_VERSION})
list(GET VERSION_LIST 0 VERSION_MAJOR)
if(${VERSION_MAJOR} GREATER "9")
target_compile_options(pdaf-d PRIVATE "-fallow-argument-mismatch")
endif()
endif()

# This will only work in newer versions of CMake
#target_link_libraries(pdaf-d PRIVATE BLAS::BLAS)
#target_link_libraries(pdaf-d PRIVATE LAPACK::LAPACK)
# So we use
target_link_libraries(pdaf-d PRIVATE ${BLAS_LIBRARIES})
target_link_libraries(pdaf-d PRIVATE ${LAPACK_LIBRARIES})

if(PDAF_USE_MPI)
#   set(MPI_Fortran_COMPILER_FLAGS,-static_mpi)
   target_compile_options(pdaf-d PUBLIC -D_PDAF_USE_MPI_)
   target_include_directories(pdaf-d PUBLIC "${MPI_Fortran_INCLUDE_PATH}")
   target_link_libraries(pdaf-d PUBLIC ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS})
else()
   target_include_directories(pdaf-d PRIVATE "${PROJECT_SOURCE_DIR}/dummympi")
endif(PDAF_USE_MPI)

install(TARGETS pdaf-d
        EXPORT pdaf-d-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/
        DESTINATION include)

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/pdaf_configure.h.in" "${CMAKE_CURRENT_BINARY_DIR}/pdaf_configure.h")

add_executable(pdaf_configure pdaf_configure.F90)
target_include_directories(pdaf_configure PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(pdaf_configure PRIVATE pdaf-d)

#install(TARGETS pdaf_configure
        #RUNTIME DESTINATION bin
#)

