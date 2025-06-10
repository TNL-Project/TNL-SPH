# Included third-party tools

### DualSPHysics GenCase

In several examples, to discretize geometries using DBC and MDBC, GenCase by [DualSPHysics](https://dual.sphysics.org/) is used. Using option cmake option `-DENABLE_GENCAGE_DOWNLOAD=True` while building the poroject, the GenCase binary is downloaded from [the DualSPHysics github repository](https://github.com/DualSPHysics/DualSPHysics/tree/master/bin/linux) and placed to tis folder. For further details about GenCase, visit [DualSPHysics webpage](https://dual.sphysics.org/).
