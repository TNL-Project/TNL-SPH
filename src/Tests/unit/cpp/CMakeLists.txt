# C++ unit tests (moved from src/UnitTests in todo 2; CTest registration
# verified in todo 26).
add_subdirectory( DistributedParticles )

add_executable( MultiresolutionNeighborLoopTest MultiresolutionNeighborLoopTest.cu )
target_link_libraries( MultiresolutionNeighborLoopTest PUBLIC TNL::TNL )
add_test( NAME MultiresolutionNeighborLoopTest COMMAND MultiresolutionNeighborLoopTest )
set_tests_properties( MultiresolutionNeighborLoopTest PROPERTIES LABELS "unit" )

add_executable( MultiresolutionBoundaryGhostRefreshTest MultiresolutionBoundaryGhostRefreshTest.cu )
target_link_libraries( MultiresolutionBoundaryGhostRefreshTest PUBLIC TNL::TNL )
add_test( NAME MultiresolutionBoundaryGhostRefreshTest COMMAND MultiresolutionBoundaryGhostRefreshTest )
set_tests_properties( MultiresolutionBoundaryGhostRefreshTest PROPERTIES LABELS "unit" )

add_executable( ParticleSetCapacityTest ParticleSetCapacityTest.cu )
target_link_libraries( ParticleSetCapacityTest PUBLIC TNL::TNL )
add_test( NAME ParticleSetCapacityTest COMMAND ParticleSetCapacityTest )
set_tests_properties( ParticleSetCapacityTest PROPERTIES LABELS "unit" )

add_executable( MultiresolutionBufferCapacityTest MultiresolutionBufferCapacityTest.cu )
target_link_libraries( MultiresolutionBufferCapacityTest PUBLIC TNL::TNL )
add_test( NAME MultiresolutionBufferCapacityTest COMMAND MultiresolutionBufferCapacityTest )
set_tests_properties( MultiresolutionBufferCapacityTest PROPERTIES LABELS "unit" )

add_executable( MultiresolutionSearchSchedulingTest MultiresolutionSearchSchedulingTest.cu )
target_link_libraries( MultiresolutionSearchSchedulingTest PUBLIC TNL::TNL )
add_test( NAME MultiresolutionSearchSchedulingTest COMMAND MultiresolutionSearchSchedulingTest )
set_tests_properties( MultiresolutionSearchSchedulingTest PROPERTIES LABELS "unit" )

add_executable( MultiresolutionGrowthSafetyTest MultiresolutionGrowthSafetyTest.cu )
target_link_libraries( MultiresolutionGrowthSafetyTest PUBLIC TNL::TNL )
add_test( NAME MultiresolutionGrowthSafetyTest COMMAND MultiresolutionGrowthSafetyTest )
set_tests_properties( MultiresolutionGrowthSafetyTest PROPERTIES LABELS "unit" )

add_executable( MultiresolutionMeasurementTest MultiresolutionMeasurementTest.cu )
target_link_libraries( MultiresolutionMeasurementTest PUBLIC TNL::TNL )
add_test( NAME MultiresolutionMeasurementTest COMMAND MultiresolutionMeasurementTest )
set_tests_properties( MultiresolutionMeasurementTest PROPERTIES LABELS "unit" )
