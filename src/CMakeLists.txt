#add_subdirectory( Benchmarks )
#add_subdirectory( tools )
add_subdirectory( UnitTests )

if( TNL_BUILD_TESTS )
   add_subdirectory( UnitTests )
endif()
