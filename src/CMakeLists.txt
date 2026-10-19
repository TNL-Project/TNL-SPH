if( TNL-SPH_BUILD_BENCHMARKS )
   add_subdirectory( Benchmarks )
endif()

# Unified test suite (todo 1 keeps UnitTests alongside until the content is migrated in todo 2)
add_subdirectory( Tests )
add_subdirectory( UnitTests )
