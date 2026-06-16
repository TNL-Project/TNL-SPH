add_subdirectory( tools )

if( TNL-SPH_BUILD_BENCHMARKS )
   add_subdirectory( Benchmarks )
endif()

add_subdirectory( UnitTests )
