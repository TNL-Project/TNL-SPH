# add all targets
add_executable(damBreak3D_WCSPH-DBC_mpi_cuda damBreak3D_WCSPH-DBC.cu)

# add TNL to all targets
target_link_libraries(damBreak3D_WCSPH-DBC_mpi_cuda PUBLIC TNL::TNL)

set(TARGETS
        damBreak3D_WCSPH-DBC_mpi_cuda
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

    # enable MPI
    find_package(MPI COMPONENTS CXX)
    if(MPI_FOUND)
        target_compile_definitions(${target} PUBLIC "-DHAVE_MPI")
        target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
    endif()
endforeach()
