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

## License

TNL-SPH is provided under the terms of the [MIT License](./LICENSE).
