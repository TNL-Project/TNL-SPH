#add_subdirectory( tools )

if( TNL-SPH_BUILD_BENCHMARKS )
   add_subdirectory( Benchmarks )
endif()

if( TNL-SPH_BUILD_TESTS )
   add_subdirectory( UnitTests )
endif()
