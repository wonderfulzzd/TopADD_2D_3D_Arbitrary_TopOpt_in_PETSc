2D_3D_Arbitrary_TopOpt_in_PETSc
===============
A 2D/3D integrated topology optimization parallel-computing framework for arbitrary design domains
===============

The code solves topology optimization problems with arbitrary design domains in either 2D or 3D. By default, it solves a 2D compliance minimization problem as shown in the below picture:

![2D problem](pic/01.jpg?raw=true "2D Example")

Moreover, it can solve 3D problems as well as two other problems: compliant design and heat conduction problems:

![Different problem](pic/02.jpg?raw=true "2D Extruded and 3D")


## The Base-Code

This code is based on the version 2017 of the TopOpt_in_PETSc code. The
original code can be found here: https://github.com/topopt/TopOpt_in_PETSc 
(Aage, N., Andreassen, E., & Lazarov, B. S. (2015). Topology optimization using PETSc: An easy-to-use, fully parallel, open source topology optimization framework. Structural and Multidisciplinary Optimization, 51(3), 565–572. https://doi.org/10.1007/s00158-014-1157-0)


## The Code Overview

![Code overview](pic/03.jpg?raw=true "Code Overview")

The green are new added files and modules.


## Compilation and Execution

Compile following rules in makefile

To compile, e.g.:make topopt

To run, e.g.: mpiexec -np 4 ./topopt

To visulize, using Paraview

> **NOTE**: The code works with **PETSc version 3.9.0**


## 3D Printing

The optimized parts can be directly printed without any post-processing:

![3D printing](pic/04.jpg?raw=true "3D Printed Parts")

