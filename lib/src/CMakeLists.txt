

ot_add_current_dir_to_include_dirs ()

ot_add_source_file ( SubsetInverseSampling.cxx )
ot_add_source_file ( SubsetInverseSamplingResult.cxx )

ot_install_header_file ( SubsetInverseSampling.hxx )
ot_install_header_file ( SubsetInverseSamplingResult.hxx )


include_directories ( ${INTERNAL_INCLUDE_DIRS} )

add_library ( otsubsetinverse ${SOURCEFILES} )
if ( BUILD_SHARED_LIBS )
  set_target_properties ( otsubsetinverse PROPERTIES COMPILE_DEFINITIONS "OTSUBSETINVERSE_DLL_EXPORTS" )
endif ()

if ( NOT DEFINED LIB_VERSION ) 
  set ( LIB_VERSION 0.0.0 )
endif ()
if ( NOT DEFINED LIB_SOVERSION ) 
  set ( LIB_SOVERSION 0 )
endif ()
set_target_properties ( otsubsetinverse PROPERTIES VERSION ${LIB_VERSION} )
set_target_properties ( otsubsetinverse PROPERTIES SOVERSION ${LIB_SOVERSION} )
target_link_libraries (otsubsetinverse ${OPENTURNS_LIBRARY})

# Add targets to the build-tree export set
export (TARGETS otsubsetinverse FILE ${PROJECT_BINARY_DIR}/OTSubsetInverse-Targets.cmake)

# Install the export set for use with the install-tree
install(EXPORT OTSubsetInverse-Targets
        DESTINATION "${OTSUBSETINVERSE_CONFIG_CMAKE_PATH}"
        COMPONENT Development)

install(TARGETS otsubsetinverse
        EXPORT  OTSubsetInverse-Targets 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION ${OTSUBSETINVERSE_LIBRARY_PATH}
        ARCHIVE DESTINATION ${OTSUBSETINVERSE_LIBRARY_PATH}
)

