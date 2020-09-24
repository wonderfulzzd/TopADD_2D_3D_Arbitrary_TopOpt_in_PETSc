#ifndef __LINEARELASTICITY__
#define __LINEARELASTICITY__

#include <fstream>
#include <iostream>
#include <math.h>
#include <petsc.h>
#include <petsc/private/dmdaimpl.h>

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
 * Modified by Zhidong Brian Zhang in May 2020, University of Waterloo
 */

class LinearCompliant {

  public:
    // Constructor
    LinearCompliant (DM da_nodes, PetscInt numLoads, Vec xPassive0,
        Vec xPassive1, Vec xPassive2); //zzd

    // Destructor
    ~LinearCompliant ();

    // Compute objective and constraints and sensitivities at once: GOOD FOR
    // SELF_ADJOINT PROBLEMS
    PetscErrorCode ComputeObjectiveConstraintsSensitivities (PetscScalar *fx,
        PetscScalar *gx, Vec dfdx, Vec dgdx, Vec xPhys, PetscScalar Emin,
        PetscScalar Emax, PetscScalar penal, PetscScalar volfrac, Vec xPassive0,
        Vec xPassive1, Vec xPassive2);

    // Restart writer
    PetscErrorCode WriteRestartFiles ();

    // Get pointer to the FE solution
    Vec GetStateField () {
      return (U[0]);
    }

    // Get pointer to DMDA
    DM GetDM () {
      return (da_nodal);
    }

    // Logical mesh
    DM da_nodal; // Nodal mesh

  private:
    // Logical mesh
    PetscInt nn[DIM]; // Number of nodes in each direction, zzd
    PetscInt ne[DIM]; // Number of elements in each direction, zzd
    PetscScalar xc[2 * DIM]; // Domain coordinates, zzd

    // Linear algebra
    Mat K; // Global stiffness matrix
    Vec *U; // Displacement vector
    Vec RHS; // Load vector
    Vec N; // Dirichlet vector (used when imposing BCs)
#if DIM == 2
    static const PetscInt nedof = 8; // zzd Number of elemental dofs
#elif DIM == 3
    static const PetscInt nedof = 24; // zzd Number of elemental dofs
#endif
    PetscScalar KE[nedof * nedof]; // Element stiffness matrix, zzd

    // Solver
    KSP ksp; // Pointer to the KSP object i.e. the linear solver+prec
    PetscInt nlvls;
    PetscScalar nu; // Possions ratio

    // Element size
    PetscScalar dx, dy, dz;

    // Number of load domains
    PetscInt nl;

    // External spring information
    Vec Sv; // spring vector

    // Extract DMDA from input one and set up a new one for linear compliant
    PetscErrorCode SetUpDMDA (DM da_nodes); //zzd

    // Update the load and boundary conditions
    PetscErrorCode SetUpLoadAndBC (Vec xPassive0, Vec xPassive1,
        Vec xPassive2, PetscInt loadStep); //zzd

    // Solve the FE problem
    PetscErrorCode SolveState (Vec xPhys, PetscScalar Emin, PetscScalar Emax,
        PetscScalar penal, PetscInt loadStep);

    // Assemble the stiffness matrix
    PetscErrorCode AssembleStiffnessMatrix (Vec xPhys, PetscScalar Emin,
        PetscScalar Emax, PetscScalar penal);

    // Start the solver
    PetscErrorCode SetUpSolver ();

#if DIM == 2
    // Routine that doesn't change the element type upon repeated calls, zzd
    PetscErrorCode DMDAGetElements_2D (DM dm, PetscInt *nel, PetscInt *nen,
        const PetscInt *e[]);

    // Methods used to assemble the element stiffness matrix, zzd
    PetscInt Quad4Isoparametric (PetscScalar *X, PetscScalar *Y, PetscScalar nu,
        PetscInt redInt, PetscScalar *ke);
    void DifferentiatedShapeFunctions_2D (PetscScalar xi, PetscScalar eta,
        PetscScalar *dNdxi, PetscScalar *dNdeta); // zzd
    PetscScalar Inverse2M (PetscScalar J[][2], PetscScalar invJ[][2]); // zzd

#elif DIM == 3

    // Routine that doesn't change the element type upon repeated calls, zzd
    PetscErrorCode DMDAGetElements_3D (DM dm, PetscInt *nel, PetscInt *nen,
        const PetscInt *e[]);

    // Methods used to assemble the element stiffness matrix, zzd
    PetscInt Hex8Isoparametric (PetscScalar *X, PetscScalar *Y, PetscScalar *Z,
        PetscScalar nu, PetscInt redInt, PetscScalar *ke);
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
