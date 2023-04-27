string(STRIP ${PDAF_VERSION} ver)
message("PDAF version: " ${ver})
message("PDAF basedir: " ${PDAF_BASE})
add_subdirectory(${ver})
