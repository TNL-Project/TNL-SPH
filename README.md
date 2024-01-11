# TNL-SPH

__TNL-SPH__ is an implementation of the __Smoothed Particle Hydrodynamics Method__ using the  [Template Numerical Library](https://gitlab.com/tnl-project/tnl). This repository contains a general framework for writing SPH-based solvers and a several simple examples that show ho to adapt the code for particular problem.

TNL-SPH is a high-performance Smoothed Particle Hydrodynamics code for simulations of free surface flow. It is verified an validated on standard SPH benchmarks, comparing the provided results with available experimental data. The main features are:

- Modular architecture with pluggable components (formulations, schemes, diffusive terms, viscous terms, boundary conditions etc).
   - $\delta$-WCSPH, Boundary Integrals WCSPH and Riemann SPH formulations
   - Inlet and outlet boundary conditions.
- Optimized and efficient framework for general particle simulations allowing implementation of other particle methods.
- MultiGPU computations based on [CUDA-aware MPI][CUDA-aware MPI]
   - 1D domain decomposition using domain overlaps

[CUDA-aware MPI]: https://developer.nvidia.com/blog/introduction-cuda-aware-mpi/

## Getting started

1. Install [Git](https://git-scm.com/) and [Git LFS](https://git-lfs.com/).

2. Clone the repository, making sure that Git submodules are initialized:

       git clone --recurse-submodules <this_repo_url>

   If you cloned the repository without the `--recurse-submodules` option,
   the submodules can be initialized subsequently:

       git submodule update --init --recursive

3. Install the necessary tools and dependencies:

    - [CMake](https://cmake.org/) build system (version 3.24 or newer)
    - [CUDA](https://docs.nvidia.com/cuda/index.html) toolkit (version 11 or newer)
    - compatible host compiler (e.g. [GCC](https://gcc.gnu.org/) or
      [Clang](https://clang.llvm.org/))
    - [Python 3](https://www.python.org/) (including development header files)
    - Python modules [NumPy](https://numpy.org/) and [Python VTK](https://pypi.org/project/vtk/)
    - [zlib](https://www.zlib.net/) (available in most Linux distributions)
    - [tinyxml2](https://github.com/leethomason/tinyxml2)
    - (optional) MPI library – for distributed computing
      (tested with [OpenMPI](https://www.open-mpi.org/))

4. Configure the build using `cmake` in the root path of the Git repository:

       cmake -B build -S . <additional_configure_options...>

   This will use `build` in the current path as the build directory.
   The path for the `-S` option corresponds to the root path of the project.
   You may use additional options to configure the build:

   - `-DCMAKE_BUILD_TYPE=<type>` where `<type>` is one of `Debug`, `Release`
   - `-DSYSTEM_TNL` – to use TNL installed on the system instead of the Git submodule
   - `-DTNL-SPH_BUILD_TESTS` – to build unit tests

5. Build the targets using `cmake`:

       cmake --build build

6. Run the example problem:

       ./examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/run.py --device cuda

   This will use the default configuration prepared in the [template](examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/template)
   file. Use the `--help` option to see the options available in `init.py` and `run.py`.

## License

TNL-SPH is provided under the terms of the [MIT License](./LICENSE).
