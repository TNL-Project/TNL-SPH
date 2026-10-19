# DistributedParticleSystemTest - C++ unit test of TNL's DistributedParticles
# with a 1D decomposition. Registered with CTest under the "unit" label
# (wired up in unified-test-suite todo 26).
#
# The test source includes <gtest/gtest.h> and <gmock/gmock.h> (via
# ../main_mpi.h) and uses MPI::Comm unconditionally, so both Google Test and
# MPI are required to build it. The project does not fetch gtest, so the
# system installation is used; skip the target with a warning if either
# dependency is missing instead of failing the whole configure.
find_package( GTest )
find_package( MPI COMPONENTS CXX )
if( NOT GTest_FOUND OR NOT MPI_FOUND )
   message( WARNING "GTest or MPI not found - DistributedParticleSystemTest will not be built." )
   return()
endif()

add_executable( DistributedParticleSystemTest DistributedParticleSystemTest.cu )
target_link_libraries( DistributedParticleSystemTest PUBLIC TNL::TNL )
# The test originates from TNL's own unit tests and includes TNL headers in the
# <Particles/...> style, i.e. relative to tnl-src/src/TNL (TNL_INCLUDE_DIRS is
# cached by tnl-src); TNL::TNL only exposes tnl-src/src (<TNL/Particles/...>).
# The include style was kept by the API port (todo 26), so add the dir as a
# SYSTEM include.
target_include_directories( DistributedParticleSystemTest SYSTEM PRIVATE ${TNL_INCLUDE_DIRS}/TNL )

# enable OpenMP (same idiom as the dummy fixtures)
find_package( OpenMP COMPONENTS CXX )
if( OPENMP_FOUND )
   target_compile_definitions( DistributedParticleSystemTest PUBLIC "-DHAVE_OPENMP" )
   target_link_libraries( DistributedParticleSystemTest PUBLIC OpenMP::OpenMP_CXX )
   # nvcc needs -Xcompiler
   target_compile_options( DistributedParticleSystemTest PUBLIC $<$<CUDA_COMPILER_ID:NVIDIA>: -Xcompiler=-fopenmp > )
endif()

# enable MPI (the test source uses MPI::Comm / TNL::MPI::GetRank unconditionally)
target_compile_definitions( DistributedParticleSystemTest PUBLIC "-DHAVE_MPI" )
target_link_libraries( DistributedParticleSystemTest PUBLIC MPI::MPI_CXX )

target_link_libraries( DistributedParticleSystemTest PUBLIC GTest::gtest GTest::gmock )

# The test asserts per-rank results for ranks 0, 1, 2 of a 3x1 subdomain
# decomposition (see SubdomainsParametersX3Y1 in DistributedParticleSystemTest.h),
# so it must run under exactly 3 MPI processes.
add_test( NAME DistributedParticleSystemTest
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:DistributedParticleSystemTest> )
set_tests_properties( DistributedParticleSystemTest PROPERTIES LABELS "unit" )
