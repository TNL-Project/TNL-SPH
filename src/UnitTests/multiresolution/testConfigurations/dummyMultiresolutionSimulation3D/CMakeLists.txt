add_executable(dummyMultiresolutionSimulation3D dummyMultiresolutionSimulation3D.cu)

target_link_libraries(dummyMultiresolutionSimulation3D PUBLIC TNL::TNL)

set(TARGETS
      dummyMultiresolutionSimulation3D
)

foreach(target IN ITEMS ${TARGETS})
    # enable OpenMP for the target
   find_package(OpenMP COMPONENTS CXX)
   if(OPENMP_FOUND)
        target_compile_definitions(${target} PUBLIC "-DHAVE_OPENMP")
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
        # nvcc needs -Xcompiler
        target_compile_options(${target} PUBLIC $<$<CUDA_COMPILER_ID:NVIDIA>: -Xcompiler=-fopenmp >)
   endif()
endforeach()
