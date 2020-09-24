#ifndef TOPOPT_H
#define TOPOPT_H

#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include "MMA.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <petsc/private/dmdaimpl.h>
#include <sstream>

#include "options.h" // framework options, zzd

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
 Updated: June 2019, Niels Aage
 Copyright (C) 2013-2019,

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
 */

/*
 *
 * Parameter container for the topology optimization problem
 *
 * min_x fx
 * s.t. gx_j <= 0, j=1..m
 *      xmin_i <= x_i <= xmax_i, i=1..n
 *
 * with filtering and a volume constraint
 *
 */

/*
 * Modified by Zhidong Brian Zhang in May 2020, University of Waterloo
 */

class TopOpt {

  public:
    // Constructor/Destructor
    TopOpt (PetscInt nconstraint);
    TopOpt ();
    ~TopOpt ();

    // Method to allocate MMA with/without restarting
    PetscErrorCode AllocateMMAwithRestart (PetscInt *itr, MMA **mma);
    PetscErrorCode WriteRestartFiles (PetscInt *itr, MMA *mma);

    // Physical domain variables
    PetscScalar xc[2 * DIM]; // Domain coordinates, zzd
    PetscScalar dx, dy, dz; // Element size
    PetscInt nxyz[DIM]; // Number of nodes in each direction, zzd
    PetscInt nlvls; // Number of multigrid levels
    PetscScalar nu; // Poisson's ratio
    /* NOTE: two meshes are needed such the both
     * nodal and element mesh share the partitioning
     */
    // Nodal mesh (basis for physics)
    DM da_nodes;
    // element mesh (basis for design)
    DM da_elem;

    // Optimization parameters
    PetscInt n; // Total number of design variables
    PetscInt nloc; // Local number of local nodes?
    PetscInt m; // Number of constraints
    PetscScalar fx; // Objective value
    PetscScalar fscale; // Scaling factor for objective
    PetscScalar *gx; // Array with constraint values
    PetscScalar Xmin; // Min. value of design variables
    PetscScalar Xmax; // Max. value of design variables

    PetscScalar movlim; // Max. change of design variables
    PetscScalar volfrac; // Volume fraction
    PetscScalar penal; // Penalization parameter
    PetscScalar Emin, Emax; // Modified SIMP, max and min E

    PetscInt maxItr; // Max iterations

    PetscScalar rmin; // filter radius
    PetscInt filter; // Filter type
    PetscBool projectionFilter; // Smooth heaviside projectionFilter
    PetscReal beta;
    PetscReal betaFinal;
    PetscReal eta;

    Vec x; // Design variables
    Vec xTilde; // Filtered field
    Vec xPhys; // Physical variables (filtered x)
    Vec dfdx; // Sensitivities of objective
    Vec xmin, xmax; // Vectors with max and min values of x
    Vec xold; // x from previous iteration
    Vec *dgdx; // Sensitivities of constraints (vector array)

    // Restart data for MMA:
    PetscBool restart, flip;
    std::string restdens_1, restdens_2;
    Vec xo1, xo2, U, L;

  private:
    // Allocate and set default values
    void Init ();
    PetscErrorCode SetUp ();

    PetscErrorCode SetUpMESH ();
    PetscErrorCode SetUpOPT ();

    // Restart filenames
    std::string filename00, filename00Itr, filename01, filename01Itr;

    // File existence
    inline PetscBool fexists (const std::string &filename) {
      std::ifstream ifile (filename.c_str ());
      if (ifile) {
        return PETSC_TRUE;
      }
      return PETSC_FALSE;
    }

    /**
     * Newly added class members
     */
  public:
    PetscScalar E; // Young's modulus
    // Optimization parameters
    PetscInt nnd; // Total number of non design variables
//    PetscInt nnd2; // The temporary test of the non-design variables
//    PetscInt nStructureEle;
//    PetscInt nSupportEle;
//    PetscScalar volfracSup; // Volume fraction of support structure
//    PetscScalar alpha1, alpha2, alpha3, alpha4; // Weight ratio of the gravity and ISM
//    PetscScalar fx2; // Objective value
//    PetscScalar fx3; // Objective value
//    PetscScalar fx4; // Objective value
//    PetscScalar fscale2;
//    PetscScalar fscale3;
//    PetscScalar fscale4;
    PetscInt numLoads; // number of load cases
    std::string *inputSTL_DES; // name of part file
    std::string *inputSTL_FIX; // name of fixture file
    std::string *inputSTL_LOD; // name of loading domain file
//    Vec xFinalDesign; // Final design distribution including structure and support
    Vec xPassive0; // the passive hole element index
    Vec xPassive1; // the passive fixture position element index
    Vec xPassive2; // the passive loading position element index
//    Vec xSupportConnection; // the connection solid element index for support opt
//    Vec dfdx2; // Sensitivities of objective
//    Vec dfdx3; // Sensitivities of objective
//    Vec dfdx4; // Sensitivities of objective
//    PetscScalar gapDis; // the gap distance between sup. and part.
    Vec nodeDensity; // node density
    Vec nodeAddingCounts; // node adding counts when summing node density from element density
};

#endif
