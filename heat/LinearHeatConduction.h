#ifndef __LINEARHEATCONDUCTION__
#define __LINEARHEATCONDUCTION__

#include <fstream>
#include <iostream>
#include <math.h>
#include <petsc.h>
#include <petsc/private/dmdaimpl.h>

#include "options.h" // framework options

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
 * Modified by Zhidong Brian Zhang in August 2020, University of Waterloo
 */

class LinearHeatConduction {

  public:
    // Constructor
    LinearHeatConduction (DM da_nodes, DM da_elem, PetscInt m,
        PetscInt numLoads,
        Vec xPassive0, Vec xPassive1, Vec xPassive2, Vec xPassive3);

    // Destructor
    ~LinearHeatConduction ();

    //  Compute objective and constraints and sensitivities at once: GOOD FOR
    //  SELF_ADJOINT PROBLEMS
    PetscErrorCode ComputeObjectiveConstraintsSensitivities (PetscScalar *fx,
        PetscScalar *gx, Vec dfdx, Vec *dgdx, Vec xPhys, PetscScalar Emin,
        PetscScalar Emax, PetscScalar penal, PetscScalar volfrac,
        Vec xPassive0, Vec xPassive1, Vec xPassive2, Vec xPassive3);

    // Restart writer
    PetscErrorCode WriteRestartFiles ();

    // Get pointer to the FE solution
    Vec GetStateField () {
      return (U);
    }

    // Get pointer to DMDA
    DM GetDM () {
      return (da_nodal);
    }

    // Logical mesh
    DM da_nodal; // Nodal mesh

  private:
    // Logical mesh
    PetscInt nn[DIM]; // Number of nodes in each direction
    PetscInt ne[DIM]; // Number of elements in each direction
    PetscScalar xc[2 * DIM]; // Domain coordinates

    // Linear algebra
    Mat K; // Global heat conduction matrix
    Vec U; // Temperature vector
    Vec *RHS; // Load vector
    Vec *N; // Dirichlet vector (used when imposing BCs)
#if DIM == 2
    static const PetscInt nedof = 4; // Number of elemental dofs
#elif DIM == 3
    static const PetscInt nedof = 8; // Number of elemental dofs
#endif

    PetscScalar KE[nedof * nedof]; // Element heat conductivity matrix

    // Solver
    KSP ksp; // Pointer to the KSP object i.e. the linear solver+prec
    PetscInt nlvls;

    // Loading conditions
    PetscInt numLoads; // # new; number of loading conditions

    // Element dimensions
    PetscScalar dx, dy, dz; // # new

    // Number of constraints
    PetscInt m; // # new

    // Set up the FE mesh and data structures
    PetscErrorCode SetUpLoadAndBC (DM da_nodes, DM da_elem, Vec xPassive0,
        Vec xPassive1, Vec xPassive2, Vec xPassive3);

    // Solve the FE problem
    PetscErrorCode SolveState (Vec xPhys, PetscScalar Emin, PetscScalar Emax,
        PetscScalar penal, PetscInt loadCondition);

    // Assemble the heat conductivity matrix
    PetscErrorCode AssembleConductivityMatrix (Vec xPhys, PetscScalar Emin,
        PetscScalar Emax, PetscScalar penal, PetscInt loadCondition);

    // Start the solver
    PetscErrorCode SetUpSolver ();

#if DIM == 2
    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_2D (DM dm, PetscInt *nel, PetscInt *nen,
        const PetscInt *e[]);

    // Methods used to assemble the element heat conductivity matrix
    PetscInt Quad4Isoparametric (PetscScalar *X, PetscScalar *Y,
        PetscInt redInt, PetscScalar *ke);
    void DifferentiatedShapeFunctions_2D (PetscScalar xi, PetscScalar eta,
        PetscScalar *dNdxi, PetscScalar *dNdeta);
    PetscScalar Inverse2M (PetscScalar J[][2], PetscScalar invJ[][2]);

#elif DIM == 3

    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_3D (DM dm, PetscInt *nel, PetscInt *nen,
        const PetscInt *e[]);

    // Methods used to assemble the element heat conductivity matrix
    PetscInt Hex8Isoparametric (PetscScalar *X, PetscScalar *Y, PetscScalar *Z,
        PetscInt redInt, PetscScalar *ke);
    void DifferentiatedShapeFunctions (PetscScalar xi, PetscScalar eta,
        PetscScalar zeta, PetscScalar *dNdxi, PetscScalar *dNdeta,
        PetscScalar *dNdzeta);
    PetscScalar Inverse3M (PetscScalar J[][3], PetscScalar invJ[][3]);
    #endif

    PetscScalar Dot (PetscScalar *v1, PetscScalar *v2, PetscInt l);

    // Restart
    PetscBool restart, flip;
    std::string filename00, filename01;

    // File existence
    inline PetscBool fexists (const std::string &filename) {
      std::ifstream ifile (filename.c_str ());
      if (ifile) {
        return PETSC_TRUE;
      }
      return PETSC_FALSE;
    }
};

#endif
