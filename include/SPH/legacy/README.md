This directory quarantines the previous unitary solver generation (`SolverMultiSet.h/.hpp`) and its obsolete `MultiresolutionBuffer` implementation.
It is retained only so the `dummyMultiresolutionSimulation2D` and `dummyMultiresolutionSimulation3D` infrastructure fixtures continue to build.
The supported solver hierarchy lives in `include/SPH/solvers/`; new code must not depend on this legacy generation.
Do not extend these files or change their algorithms while they are quarantined.
Their eventual removal is a separate future decision.
