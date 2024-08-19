#set( CPP_TESTS
#         ParticleSystemTest
#)

set( CUDA_TESTS
         DistributedParticleSystemTest
)

#foreach( target IN ITEMS ${CPP_TESTS} )
#   add_executable( ${target} ${target}.cpp )
#   target_compile_options( ${target} PUBLIC ${CXX_TESTS_FLAGS} )
#   target_link_libraries( ${target} PUBLIC TNL::TNL ${TESTS_LIBRARIES} )
#   target_link_options( ${target} PUBLIC ${TESTS_LINKER_FLAGS} )
#   add_test( ${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target}${CMAKE_EXECUTABLE_SUFFIX} )
#endforeach()

if( TNL_BUILD_CUDA )
   foreach( target IN ITEMS ${CUDA_TESTS} )
      add_executable( ${target} ${target}.cu )
      target_compile_options( ${target} PUBLIC ${CUDA_TESTS_FLAGS} )
      target_link_libraries( ${target} PUBLIC TNL::TNL ${TESTS_LIBRARIES} )
      target_link_options( ${target} PUBLIC ${TESTS_LINKER_FLAGS} )
      add_test( ${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target}${CMAKE_EXECUTABLE_SUFFIX} )

      # enable MPI
      find_package(MPI COMPONENTS CXX)
      if(MPI_FOUND)
          target_compile_definitions(${target} PUBLIC "-DHAVE_MPI")
          target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
      endif()
   endforeach()
endif()
