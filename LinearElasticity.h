#ifndef __LINEARELASTICITY__
#define __LINEARELASTICITY__

#include <fstream>
#include <iostream>
#include <math.h>
#include <petsc.h>
#include <petsc/private/dmdaimpl.h>

#include "options.h" // # new; framework options

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

class LinearElasticity {

  public:
    // Constructor
    LinearElasticity (DM da_nodes, PetscInt numLoads, Vec xPassive0,
        Vec xPassive1, Vec xPassive2, Vec xPassive3); //zzd

    // Destructor
    ~LinearElasticity ();

    // Compute objective and constraints and sensitivities at once: GOOD FOR
    // SELF_ADJOINT PROBLEMS
    PetscErrorCode ComputeObjectiveConstraintsSensitivities (PetscScalar *fx,
        PetscScalar *gx, Vec dfdx, Vec dgdx, Vec xPhys, PetscScalar Emin,
        PetscScalar Emax, PetscScalar penal, PetscScalar volfrac, Vec xPassive0,
        Vec xPassive1, Vec xPassive2, Vec xPassive3); // # modified

    // Compute objective and constraints for the optimiation
    PetscErrorCode ComputeObjectiveConstraints (PetscScalar *fx,
        PetscScalar *gx, Vec xPhys, PetscScalar Emin, PetscScalar Emax,
        PetscScalar penal, PetscScalar volfrac, Vec xPassive0,
        Vec xPassive1, Vec xPassive2, Vec xPassive3); // # modified

    // Compute sensitivities
    PetscErrorCode ComputeSensitivities (Vec dfdx, Vec dgdx, Vec xPhys,
        PetscScalar Emin, PetscScalar Emax, PetscScalar penal,
        PetscScalar volfrac); // needs ....

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
    PetscInt nn[DIM]; // # modified; Number of nodes in each direction
    PetscInt ne[DIM]; // # modified; Number of elements in each direction
    PetscScalar xc[2 * DIM]; // # modified; Domain coordinates

    // Linear algebra
    Mat K; // Global stiffness matrix
    Vec U; // Displacement vector
    Vec RHS; // Load vector
    Vec N; // Dirichlet vector (used when imposing BCs)
#if DIM == 2  // # new
    static const PetscInt nedof = 8; // Number of elemental dofs
#elif DIM == 3  // # new
    static const PetscInt nedof = 24; // Number of elemental dofs
#endif
    PetscScalar KE[nedof * nedof]; // # new; Element stiffness matrix

    // Solver
    KSP ksp; // Pointer to the KSP object i.e. the linear solver+prec
    PetscInt nlvls;
    PetscScalar nu; // Possions ratio

    // Number of load domains
    PetscInt nl;  // # new

    // Element dimensions
    PetscScalar dx, dy, dz;  // # new

    // Set up the FE mesh and data structures
    PetscErrorCode SetUpNodalMesh (DM da_nodes); // # modified

    // Set up the FE mesh and data structures
    PetscErrorCode SetUpLoadAndBC (Vec xPassive0, Vec xPassive1, Vec xPassive2,
        Vec xPassive3); // # modified

    // Solve the FE problem
    PetscErrorCode SolveState (Vec xPhys, PetscScalar Emin, PetscScalar Emax,
        PetscScalar penal);

    // Assemble the stiffness matrix
    PetscErrorCode AssembleStiffnessMatrix (Vec xPhys, PetscScalar Emin,
        PetscScalar Emax, PetscScalar penal);

    // Start the solver
    PetscErrorCode SetUpSolver ();

#if DIM == 2    // # new
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
