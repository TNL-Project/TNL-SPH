# add all targets
add_executable(lidDrivenCavity2D_SHTC_cuda lidDrivenCavity2D_SHTC.cu)

# add TNL to all targets
target_link_libraries(lidDrivenCavity2D_SHTC_cuda PUBLIC TNL::TNL)

set(TARGETS
        lidDrivenCavity2D_SHTC_cuda
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
