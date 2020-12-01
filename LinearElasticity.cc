#include "LinearElasticity.h"

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
 */

/*
 * Modified by Zhidong Brian Zhang in May 2020, University of Waterloo
 */

LinearElasticity::LinearElasticity (DM da_nodes, PetscInt m, PetscInt numLoads,
    PetscInt numNodeLoadAddingCounts, PetscScalar nu, PetscScalar E,
    PetscScalar *loadVector, Vec xPassive0, Vec xPassive1, Vec xPassive2,
    Vec xPassive3) { // # modified
  // Set pointers to null
  K = NULL;
  U = NULL;
  RHS = NULL;
  N = NULL;
  ksp = NULL;
  da_nodal = NULL;

  // Parameters - to be changed on read of variables
  this->nu = nu; // # modified
  this->E = E; // # new
  nlvls = 4;
  PetscBool flg;
  PetscOptionsGetInt (NULL, NULL, "-nlvls", &nlvls, &flg);
  PetscOptionsGetReal (NULL, NULL, "-nu", &nu, &flg);

  this->m = m; // # new
  this->loadVector = NULL; // # new
  this->numLoads = numLoads; // # new; num of loads, save for internal uses
  this->numNodeLoadAddingCounts = numNodeLoadAddingCounts; // # new; num of node load adding counts
  this->loadVector = new PetscScalar[numLoads * DIM]; // # new; total body load (e.g. gravity accleration)

  for (PetscInt i = 0; i < this->numLoads * DIM; ++i) { // # new
    if (this->numNodeLoadAddingCounts != 0) { // # new
      this->loadVector[i] = loadVector[i] / this->numNodeLoadAddingCounts; // # new
    } else { // # new
      this->loadVector[i] = loadVector[i]; // # new
    } // # new
  } // # new

  RHS = new Vec[numLoads]; // # new
  N = new Vec[numLoads]; // # new

  // Setup sitffness matrix, load vector and bcs (Dirichlet) for the design
  // problem
  SetUpLoadAndBC (da_nodes, xPassive0, xPassive1, xPassive2, xPassive3); // # modified
}

LinearElasticity::~LinearElasticity () {
  // Deallocate
  VecDestroy (&(U)); // # modified
  VecDestroyVecs (numLoads, &(RHS)); // # modified
  VecDestroyVecs (numLoads, &(N)); // # modified
  MatDestroy (&(K));
  KSPDestroy (&(ksp));

  if (da_nodal != NULL) {
    DMDestroy (&(da_nodal));
  }
  if (loadVector != NULL) delete loadVector; // # new
}

PetscErrorCode LinearElasticity::SetUpLoadAndBC (DM da_nodes, Vec xPassive0,
    Vec xPassive1, Vec xPassive2, Vec xPassive3) {
  PetscErrorCode ierr = 0;

#if DIM == 2  // # new
  // Extract information from input DM and create one for the linear elasticity
  // number of nodal dofs: (u,v)
  PetscInt numnodaldof = 2;

  // Stencil width: each node connects to a box around it - linear elements
  PetscInt stencilwidth = 1;

  PetscScalar dx, dy;
  DMBoundaryType bx, by;
  DMDAStencilType stype;
  {
    // Extract information from the nodal mesh
    PetscInt M, N, md, nd;
    DMDAGetInfo (da_nodes, NULL, &M, &N, NULL, &md, &nd, NULL, NULL, NULL, &bx,
        &by, NULL, &stype);

    // Find the element size
    Vec lcoor;
    DMGetCoordinatesLocal (da_nodes, &lcoor);
    PetscScalar *lcoorp;
    VecGetArray (lcoor, &lcoorp);

    PetscInt nel, nen;
    const PetscInt *necon;
    DMDAGetElements_2D (da_nodes, &nel, &nen, &necon);

    // Use the first element to compute the dx, dy, dz
    dx = lcoorp[DIM * necon[0 * nen + 1] + 0]
         - lcoorp[DIM * necon[0 * nen + 0] + 0];
    dy = lcoorp[DIM * necon[0 * nen + 2] + 1]
         - lcoorp[DIM * necon[0 * nen + 1] + 1];
    VecRestoreArray (lcoor, &lcoorp);

    nn[0] = M;
    nn[1] = N;

    ne[0] = nn[0] - 1;
    ne[1] = nn[1] - 1;

    xc[0] = 0.0;
    xc[1] = ne[0] * dx;
    xc[2] = 0.0;
    xc[3] = ne[1] * dy;
  }

  // Create the nodal mesh
  DMDACreate2d (PETSC_COMM_WORLD, bx, by, stype, nn[0], nn[1], PETSC_DECIDE,
  PETSC_DECIDE, numnodaldof, stencilwidth, 0, 0, &(da_nodal));
  // Initialize
  DMSetFromOptions (da_nodal);
  DMSetUp (da_nodal);

  // Set the coordinates
  DMDASetUniformCoordinates (da_nodal, xc[0], xc[1], xc[2], xc[3], 0.0, 0.0);
  // Set the element type to Q1: Otherwise calls to GetElements will change to
  // P1 ! STILL DOESN*T WORK !!!!
  DMDASetElementType (da_nodal, DMDA_ELEMENT_Q1);

  // Allocate matrix and the RHS and Solution vector and Dirichlet vector
  ierr = DMCreateMatrix (da_nodal, &(K));
  CHKERRQ(ierr);
  ierr = DMCreateGlobalVector (da_nodal, &(U));
  CHKERRQ(ierr);
  VecDuplicateVecs (U, numLoads, &(RHS));
  VecDuplicateVecs (U, numLoads, &(N));

  // Set the local stiffness matrix
  PetscScalar X[4] = { 0.0, dx, dx, 0.0 };
  PetscScalar Y[4] = { 0.0, 0.0, dy, dy };

  // Compute the element stiffnes matrix - constant due to structured grid
  Quad4Isoparametric (X, Y, nu, false, KE);

  // Save the element size for other uses
  this->dx = dx;
  this->dy = dy;

  // Set the RHS and Dirichlet vector
  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) {
    VecSet (N[loadCondition], 1.0);
    VecSet (RHS[loadCondition], 0.0);
  }

  // Global coordinates and a pointer
  Vec lcoor; // borrowed ref - do not destroy!
  PetscScalar *lcoorp;

  // Get local coordinates in local node numbering including ghosts
  ierr = DMGetCoordinatesLocal (da_nodal, &lcoor);
  CHKERRQ(ierr);
  VecGetArray (lcoor, &lcoorp);

  // Get local dof number
  PetscInt nn;
  VecGetSize (lcoor, &nn);

  // Compute epsilon parameter for finding points in space:
  PetscScalar epsi = PetscMin(dx * 0.05, dy * 0.05);
  // Passive design variable vector
  PetscScalar *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p;
  VecGetArray (xPassive0, &xPassive0p);
  VecGetArray (xPassive1, &xPassive1p);
  VecGetArray (xPassive2, &xPassive2p);
  VecGetArray (xPassive3, &xPassive3p);

  // Set the RHS and Dirichlet vector
  PetscScalar rhs_ele[8]; // local rhs
  PetscScalar n_ele[8]; // local n
  PetscInt edof[8];

  // Find the element size
  PetscInt nel, nen;
  const PetscInt *necon;
  DMDAGetElements_2D (da_nodes, &nel, &nen, &necon);

  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) {
    if (IMPORT_GEO == 0) {
      // Set the values:
      // In this case: N = the wall at x=xmin is fully clamped
      // RHS(z) = -0.001 at the middle point of the right boundary
      PetscScalar LoadIntensity = -0.079;
      for (PetscInt i = 0; i < nn; i++) {
        // Make a wall with all dofs clamped
        if (i % DIM == 0 && PetscAbsScalar (lcoorp[i] - xc[0]) < epsi) {
          VecSetValueLocal (N[loadCondition], i, 0.0, INSERT_VALUES);
          VecSetValueLocal (N[loadCondition], ++i, 0.0, INSERT_VALUES);
        }
        // Point load
        if (i % DIM == 0 && PetscAbsScalar (lcoorp[i] - xc[1]) < epsi
            && PetscAbsScalar (lcoorp[i + 1] - xc[2]) < epsi) {
          VecSetValueLocal (RHS[loadCondition], i + 1, LoadIntensity,
              INSERT_VALUES);
        }
      }
    } else { // # new
      // Set the values:
      // In this case:
      // xPassive1 indicates fix,
      // xPassive2 indicates loading.
      // Load and constraints
      for (PetscInt i = 0; i < nel; i++) {

        memset (rhs_ele, 0.0, sizeof(rhs_ele[0]) * 8);
        memset (n_ele, 0.0, sizeof(n_ele[0]) * 8);

        // Global dof in the RHS vector
        for (PetscInt l = 0; l < nen; l++) {
          for (PetscInt m = 0; m < DIM; m++) {
            edof[l * 2 + m] = 2 * necon[i * nen + l] + m; // dof in globe
          }
        }

        if (std::fmod ((xPassive2p[i] / std::pow (2.0, loadCondition)), 2) >= 1.0) {
          for (PetscInt j = 0; j < 8; j++) {
            rhs_ele[j] = loadVector[DIM * loadCondition + j % DIM];
          }
          ierr = VecSetValuesLocal (RHS[loadCondition], 8, edof, rhs_ele,
              ADD_VALUES);
          CHKERRQ(ierr);
        }

        if (std::fmod ((xPassive1p[i] / std::pow (2.0, loadCondition)), 2) >= 1.0) {
          for (PetscInt j = 0; j < 8; j++) {
            n_ele[j] = 0.0;
          }
          ierr = VecSetValuesLocal (N[loadCondition], 8, edof, n_ele,
              INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
    VecAssemblyBegin (N[loadCondition]); // # modified
    VecAssemblyEnd (N[loadCondition]); // # modified
    VecAssemblyBegin (RHS[loadCondition]); // # modified
    VecAssemblyEnd (RHS[loadCondition]); // # modified
  }

#elif DIM == 3
  // Extract information from input DM and create one for the linear elasticity
  // number of nodal dofs: (u,v,w)
  PetscInt numnodaldof = 3;

  // Stencil width: each node connects to a box around it - linear elements
  PetscInt stencilwidth = 1;

  PetscScalar dx, dy, dz;
  DMBoundaryType bx, by, bz;
  DMDAStencilType stype;
  {
    // Extract information from the nodal mesh
    PetscInt M, N, P, md, nd, pd;
    DMDAGetInfo (da_nodes, NULL, &M, &N, &P, &md, &nd, &pd, NULL, NULL, &bx,
        &by, &bz, &stype);

    // Find the element size
    Vec lcoor;
    DMGetCoordinatesLocal (da_nodes, &lcoor);
    PetscScalar *lcoorp;
    VecGetArray (lcoor, &lcoorp);

    PetscInt nel, nen;
    const PetscInt *necon;
    DMDAGetElements_3D (da_nodes, &nel, &nen, &necon);

    // Use the first element to compute the dx, dy, dz
    dx = lcoorp[3 * necon[0 * nen + 1] + 0]
         - lcoorp[3 * necon[0 * nen + 0] + 0];
    dy = lcoorp[3 * necon[0 * nen + 2] + 1]
         - lcoorp[3 * necon[0 * nen + 1] + 1];
    dz = lcoorp[3 * necon[0 * nen + 4] + 2]
         - lcoorp[3 * necon[0 * nen + 0] + 2];
    VecRestoreArray (lcoor, &lcoorp);

    nn[0] = M;
    nn[1] = N;
    nn[2] = P;

    ne[0] = nn[0] - 1;
    ne[1] = nn[1] - 1;
    ne[2] = nn[2] - 1;

    xc[0] = 0.0;
    xc[1] = ne[0] * dx;
    xc[2] = 0.0;
    xc[3] = ne[1] * dy;
    xc[4] = 0.0;
    xc[5] = ne[2] * dz;
  }

  // Create the nodal mesh
  DMDACreate3d (PETSC_COMM_WORLD, bx, by, bz, stype, nn[0], nn[1], nn[2],
  PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, numnodaldof, stencilwidth, 0, 0, 0,
      &(da_nodal));
  // Initialize
  DMSetFromOptions (da_nodal);
  DMSetUp (da_nodal);

  // Set the coordinates
  DMDASetUniformCoordinates (da_nodal, xc[0], xc[1], xc[2], xc[3], xc[4],
      xc[5]);
  // Set the element type to Q1: Otherwise calls to GetElements will change to
  // P1 ! STILL DOESN*T WORK !!!!
  DMDASetElementType (da_nodal, DMDA_ELEMENT_Q1);

  // Allocate matrix and the RHS and Solution vector and Dirichlet vector
  ierr = DMCreateMatrix (da_nodal, &(K));
  CHKERRQ(ierr);
  ierr = DMCreateGlobalVector (da_nodal, &(U));
  CHKERRQ(ierr);
  VecDuplicateVecs (U, numLoads, &(RHS)); // # modified
  VecDuplicateVecs (U, numLoads, &(N)); // # modified

  // Set the local stiffness matrix
  PetscScalar X[8] = { 0.0, dx, dx, 0.0, 0.0, dx, dx, 0.0 };
  PetscScalar Y[8] = { 0.0, 0.0, dy, dy, 0.0, 0.0, dy, dy };
  PetscScalar Z[8] = { 0.0, 0.0, 0.0, 0.0, dz, dz, dz, dz };

  // Compute the element stiffnes matrix - constant due to structured grid
  Hex8Isoparametric (X, Y, Z, nu, false, KE);

  // # new; Save the element size for other uses
  this->dx = dx;
  this->dy = dy;
  this->dz = dz;

  // Set the RHS and Dirichlet vector
  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) { // # new
    VecSet (N[loadCondition], 1.0); // # modified
    VecSet (RHS[loadCondition], 0.0); // # modified
  } // # new

  // Global coordinates and a pointer
  Vec lcoor; // borrowed ref - do not destroy!
  PetscScalar *lcoorp;

  // Get local coordinates in local node numbering including ghosts
  ierr = DMGetCoordinatesLocal (da_nodal, &lcoor);
  CHKERRQ(ierr);
  VecGetArray (lcoor, &lcoorp);

  // Get local dof number
  PetscInt nn;
  VecGetSize (lcoor, &nn);

  // Compute epsilon parameter for finding points in space:
  PetscScalar epsi = PetscMin(dx * 0.05, PetscMin(dy * 0.05, dz * 0.05));

  // # new; Passive design variable vector
  PetscScalar *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p;
  VecGetArray (xPassive0, &xPassive0p);
  VecGetArray (xPassive1, &xPassive1p);
  VecGetArray (xPassive2, &xPassive2p);
  VecGetArray (xPassive3, &xPassive3p);

  // # new; Set the local RHS and Dirichlet vector
  PetscScalar rhs_ele[24]; // local rhs
  PetscScalar n_ele[24]; // local n
  PetscInt edof[24];

  // # new; Find the element size
  PetscInt nel, nen;
  const PetscInt *necon;
  DMDAGetElements_3D (da_nodal, &nel, &nen, &necon);

  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) {
    if (IMPORT_GEO == 0) {
      // Set the values:
      // In this case: N = the wall at x=xmin is fully clamped
      //               RHS(z) = sin(pi*y/Ly) at x=xmax,z=zmin;
      // OR
      //               RHS(z) = -0.1 at x=xmax,z=zmin;
      PetscScalar LoadIntensity = -0.001;
      for (PetscInt i = 0; i < nn; i++) {
        // Make a wall with all dofs clamped
        if (i % 3 == 0 && PetscAbsScalar (lcoorp[i] - xc[0]) < epsi) {
          VecSetValueLocal (N[loadCondition], i, 0.0, INSERT_VALUES); // # modified
          VecSetValueLocal (N[loadCondition], i + 1, 0.0, INSERT_VALUES); // # modified
          VecSetValueLocal (N[loadCondition], i + 2, 0.0, INSERT_VALUES); // # modified
        }
        // Line load
        if (i % 3 == 0 && PetscAbsScalar (lcoorp[i] - xc[1]) < epsi
            && PetscAbsScalar (lcoorp[i + 2] - xc[4]) < epsi) {
          VecSetValueLocal (RHS[loadCondition], i + 2, LoadIntensity,
              INSERT_VALUES); // # modified
        }
      }

      VecAssemblyBegin (RHS[loadCondition]); // # new
      VecAssemblyEnd (RHS[loadCondition]); // # new

      for (PetscInt i = 0; i < nn; i++) {
        // Adjust the corners
        if (i % 3 == 0 && PetscAbsScalar (lcoorp[i] - xc[1]) < epsi
            && PetscAbsScalar (lcoorp[i + 1] - xc[2]) < epsi
            && PetscAbsScalar (lcoorp[i + 2] - xc[4] / 2.0) < epsi) {
          VecSetValueLocal (RHS[loadCondition], i + 2, LoadIntensity / 2.0,
              INSERT_VALUES); // # modified
        }
        if (i % 3 == 0 && PetscAbsScalar (lcoorp[i] - xc[1]) < epsi
            && PetscAbsScalar (lcoorp[i + 1] - xc[3]) < epsi
            && PetscAbsScalar (lcoorp[i + 2] - xc[4] / 2.0) < epsi) {
          VecSetValueLocal (RHS[loadCondition], i + 2, LoadIntensity / 2.0,
              INSERT_VALUES); // # modified
        }
      }
    } else { // # new
      // Set the values:
      // In this case:
      // xPassive1 indicates fix,
      // xPassive2 indicates loading.
      // Load and constraints
      for (PetscInt i = 0; i < nel; i++) {
        //std::fill_n(rhs_ele, 24, 0); // initialization
        memset (rhs_ele, 0.0, sizeof(rhs_ele[0]) * 24);
        memset (n_ele, 0.0, sizeof(n_ele[0]) * 24);

        // Global dof in the RHS vector
        for (PetscInt l = 0; l < nen; l++) {
          for (PetscInt m = 0; m < 3; m++) {
            edof[l * 3 + m] = 3 * necon[i * nen + l] + m; // dof in globe
          }
        }

        if (std::fmod ((xPassive2p[i] / std::pow (2.0, loadCondition)), 2) >= 1.0) {
          for (PetscInt j = 0; j < 24; j++) {
            rhs_ele[j] = loadVector[DIM * loadCondition + j % DIM];
          }
          ierr = VecSetValuesLocal (RHS[loadCondition], 24, edof, rhs_ele,
              ADD_VALUES);
          CHKERRQ(ierr);
        }

        if (std::fmod ((xPassive1p[i] / std::pow (2.0, loadCondition)), 2) >= 1.0) {
          for (PetscInt j = 0; j < 24; j++) {
            n_ele[j] = 0.0;
          }
          ierr = VecSetValuesLocal (N[loadCondition], 24, edof, n_ele,
              INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
    VecAssemblyBegin (N[loadCondition]); // # modified
    VecAssemblyEnd (N[loadCondition]); // # modified
    VecAssemblyBegin (RHS[loadCondition]); // # modified
    VecAssemblyEnd (RHS[loadCondition]); // # modified
  }

#endif

  // Restore vectors
  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) { // # new
    VecAssemblyBegin (N[loadCondition]); // # modified
    VecAssemblyEnd (N[loadCondition]); // # modified
    VecAssemblyBegin (RHS[loadCondition]); // # modified
    VecAssemblyEnd (RHS[loadCondition]); // # modified
  }
  VecRestoreArray (lcoor, &lcoorp);
  DMDARestoreElements (da_nodal, &nel, &nen, &necon); // # new
  VecRestoreArray (xPassive0, &xPassive0p); // # new
  VecRestoreArray (xPassive1, &xPassive1p); // # new
  VecRestoreArray (xPassive2, &xPassive2p); // # new
  VecRestoreArray (xPassive3, &xPassive3p); // # new

  return ierr;
}

PetscErrorCode LinearElasticity::SolveState (Vec xPhys, PetscScalar Emin,
    PetscScalar Emax, PetscScalar penal, PetscInt loadCondition) {

  PetscErrorCode ierr;

  double t1, t2;
  t1 = MPI_Wtime ();

  // Assemble the stiffness matrix
  ierr = AssembleStiffnessMatrix (xPhys, Emin, Emax, penal, loadCondition);
  CHKERRQ(ierr);

  // Setup the solver
  if (ksp == NULL) {
    ierr = SetUpSolver ();
    CHKERRQ(ierr);
  } else {
    ierr = KSPSetOperators (ksp, K, K);
    CHKERRQ(ierr);
    KSPSetUp (ksp);
  }

  // Solve
  ierr = KSPSolve (ksp, RHS[loadCondition], U);
  CHKERRQ(ierr);
  CHKERRQ(ierr);

  // DEBUG
  // Get iteration number and residual from KSP
  PetscInt niter;
  PetscScalar rnorm;
  KSPGetIterationNumber (ksp, &niter);
  KSPGetResidualNorm (ksp, &rnorm);
  PetscReal RHSnorm;
  ierr = VecNorm (RHS[loadCondition], NORM_2, &RHSnorm);
  CHKERRQ(ierr);
  rnorm = rnorm / RHSnorm;

  t2 = MPI_Wtime ();
  PetscPrintf (PETSC_COMM_WORLD,
      "State solver:  iter: %i, rerr.: %e, time: %f\n", niter, rnorm, t2 - t1);

  return ierr;
}

PetscErrorCode LinearElasticity::ComputeObjectiveConstraintsSensitivities (
    PetscScalar *fx, PetscScalar *gx, Vec dfdx, Vec *dgdx, Vec xPhys,
    PetscScalar Emin, PetscScalar Emax, PetscScalar penal, PetscScalar volfrac,
    Vec xPassive0, Vec xPassive1, Vec xPassive2, Vec xPassive3) { // # modified
  // Errorcode
  PetscErrorCode ierr;

  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) { // # new
    // Solve state eqs
    ierr = SolveState (xPhys, Emin, Emax, penal, loadCondition); // # modified
    CHKERRQ(ierr);

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt nel, nen;
    const PetscInt *necon;
#if DIM == 2    // # new
    ierr = DMDAGetElements_2D (da_nodal, &nel, &nen, &necon);
    CHKERRQ(ierr);
#elif DIM == 3
    ierr = DMDAGetElements_3D (da_nodal, &nel, &nen, &necon);
    CHKERRQ(ierr);
#endif
    // DMDAGetElements(da_nodes,&nel,&nen,&necon); // Still issue with elemtype
    // change !

    // Get pointer to the densities
    PetscScalar *xp, *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p; // # modified
    VecGetArray (xPhys, &xp);
    VecGetArray (xPassive0, &xPassive0p); // # new
    VecGetArray (xPassive1, &xPassive1p); // # new
    VecGetArray (xPassive2, &xPassive2p); // # new
    VecGetArray (xPassive3, &xPassive3p); // # new

    // Get Solution
    Vec Uloc;
    DMCreateLocalVector (da_nodal, &Uloc);
    DMGlobalToLocalBegin (da_nodal, U, INSERT_VALUES, Uloc);
    DMGlobalToLocalEnd (da_nodal, U, INSERT_VALUES, Uloc);

    // get pointer to local vector
    PetscScalar *up;
    VecGetArray (Uloc, &up);

    // Get dfdx
    PetscScalar *df;
    VecGetArray (dfdx, &df);

    // # new; Get dgdx
    PetscScalar **dg;
    for (PetscInt i = 0; i < m; ++i) {
      VecSet (dgdx[i], 0);
      gx[i] = 0;
    }
    VecGetArrays (dgdx, m, &dg);

    // Number of total elements and nonDesign domain elements
    PetscInt neltot = 0;
    PetscScalar nNonDesign = 0; // # new
    VecGetSize (xPhys, &neltot); // # modified

    // Edof array
    PetscInt edof[nedof]; // # modified

    fx[0] = 0.0;
    // # modified; Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
      // loop over element nodes
      if (xPassive0p[i] == 0 && xPassive1p[i] == 0 && xPassive2p[i] == 0
          && xPassive3p[i] == 0) {
        for (PetscInt j = 0; j < nen; j++) {
          // Get local dofs
          for (PetscInt k = 0; k < DIM; k++) {
            edof[j * DIM + k] = DIM * necon[i * nen + j] + k;
          }
        }
        // Use SIMP for stiffness interpolation
        PetscScalar uKu = 0.0;
        for (PetscInt k = 0; k < nedof; k++) {
          for (PetscInt h = 0; h < nedof; h++) {
            uKu += up[edof[k]] * KE[k * nedof + h] * up[edof[h]];
          }
        }
        // Add to objective
        fx[0] += (Emin + PetscPowScalar(xp[i], penal) * (Emax - Emin)) * uKu;
        // Set the Senstivity
        df[i] = -1.0 * penal * PetscPowScalar(xp[i], penal - 1) * (Emax - Emin)
                * uKu;
        // # new; Constraints
        for (PetscInt j = 0; j < m; ++j) {
          gx[j] += xp[i];
          dg[j][i] = 1;
        }
      } else if (xPassive0p[i] == 1) { // # new
        df[i] = 1.0E9; // # new
        nNonDesign += 1; // # new
      } else if (xPassive1p[i] != 0 || xPassive2p[i] != 0
                 || xPassive3p[i] != 0) { // # new
        df[i] = -1.0E9; // # new
        nNonDesign += 1; // # new
      }
    }

    // Allreduce fx[0]
    PetscScalar tmp = fx[0];
    fx[0] = 0.0;
    MPI_Allreduce(&tmp, &(fx[0]), 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);

    tmp = nNonDesign;
    nNonDesign = 0.0;
    MPI_Allreduce(&tmp, &(nNonDesign), 1, MPIU_SCALAR, MPI_SUM,
        PETSC_COMM_WORLD);

    // # modified; Allreduce gx
    for (PetscInt i = 0; i < m; ++i) {
      tmp = gx[i];
      gx[i] = 0.0;
      MPI_Allreduce(&tmp, &(gx[i]), 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);
      gx[i] = gx[i]
              / ((PetscScalar) neltot - nNonDesign)
              - volfrac; // # modified
      VecScale (dgdx[i],
          1.0 / ((PetscScalar) neltot - nNonDesign)); // # modified
    }

    VecRestoreArray (xPhys, &xp);
    VecRestoreArray (xPassive0, &xPassive0p); // # new
    VecRestoreArray (xPassive1, &xPassive1p); // # new
    VecRestoreArray (xPassive2, &xPassive2p); // # new
    VecRestoreArray (xPassive3, &xPassive3p); // # new
    VecRestoreArray (Uloc, &up);
    VecRestoreArray (dfdx, &df);
    VecRestoreArrays (dgdx, m, &dg);
    VecDestroy (&Uloc);

  } // # new

  return (ierr);
}

PetscErrorCode LinearElasticity::ComputeObjectiveConstraints (PetscScalar *fx,
    PetscScalar *gx, Vec xPhys, PetscScalar Emin, PetscScalar Emax,
    PetscScalar penal, PetscScalar volfrac, Vec xPassive0, Vec xPassive1,
    Vec xPassive2, Vec xPassive3) {

  // Error code
  PetscErrorCode ierr;

  for (PetscInt loadCondition = 0; loadCondition < numLoads; ++loadCondition) { // # new
    // Solve state eqs
    VecSet (U, 0.0);
    ierr = SolveState (xPhys, Emin, Emax, penal, loadCondition); // # modified
    CHKERRQ(ierr);

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt nel, nen;
    const PetscInt *necon;
#if DIM == 2   // # new
    ierr = DMDAGetElements_2D (da_nodal, &nel, &nen, &necon);
    CHKERRQ(ierr);
#elif DIM == 3
    ierr = DMDAGetElements_3D (da_nodal, &nel, &nen, &necon);
    CHKERRQ(ierr);
#endif

    // Get pointer to the densities
    PetscScalar *xp, *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p; // # modified
    VecGetArray (xPhys, &xp);
    VecGetArray (xPassive0, &xPassive0p); // # new
    VecGetArray (xPassive1, &xPassive1p); // # new
    VecGetArray (xPassive2, &xPassive2p); // # new
    VecGetArray (xPassive3, &xPassive3p); // # new

    // Get Solution
    Vec Uloc;
    DMCreateLocalVector (da_nodal, &Uloc);
    DMGlobalToLocalBegin (da_nodal, U, INSERT_VALUES, Uloc);
    DMGlobalToLocalEnd (da_nodal, U, INSERT_VALUES, Uloc);

    // get pointer to local vector
    PetscScalar *up;
    VecGetArray (Uloc, &up);

    // # new; Get gx
    for (PetscInt i = 0; i < m; ++i) {
      gx[i] = 0;
    }

    // Number of total elements and nonDesign domain elements
    PetscInt neltot = 0;
    PetscScalar nNonDesign = 0; // # new
    VecGetSize (xPhys, &neltot); // # modified

    // Edof array
    PetscInt edof[nedof]; // # modified

    fx[0] = 0.0;
    // # modified; Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
      // loop over element nodes
      if (xPassive0p[i] == 0 && xPassive1p[i] == 0 && xPassive2p[i] == 0
          && xPassive3p[i] == 0) {
        for (PetscInt j = 0; j < nen; j++) {
          // Get local dofs
          for (PetscInt k = 0; k < DIM; k++) {
            edof[j * DIM + k] = DIM * necon[i * nen + j] + k;
          }
        }
        // # modified; Use SIMP for stiffness interpolation
        PetscScalar uKu = 0.0;
        for (PetscInt k = 0; k < nedof; k++) {
          for (PetscInt h = 0; h < nedof; h++) {
            uKu += up[edof[k]] * KE[k * nedof + h] * up[edof[h]];
          }
        }
        // Add to objective
        fx[0] += (Emin + PetscPowScalar(xp[i], penal) * (Emax - Emin)) * uKu;
        // # new; Constraints
        for (PetscInt j = 0; j < m; ++j) {
          gx[j] += xp[i];
        }
      } else if (xPassive0p[i] == 1) { // # new
        nNonDesign += 1; // # new
      } else if (xPassive1p[i] >= 1 || xPassive2p[i] >= 1
                 || xPassive3p[i] == 1) { // # new
        nNonDesign += 1; // # new
      }

      // Allreduce fx[0]
      PetscScalar tmp = fx[0];
      fx[0] = 0.0;
      MPI_Allreduce(&tmp, &(fx[0]), 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);

      tmp = nNonDesign; // # new
      nNonDesign = 0.0; // # new
      MPI_Allreduce(&tmp, &(nNonDesign), 1, MPIU_SCALAR, MPI_SUM,
          PETSC_COMM_WORLD); // # new
      // # modified; Allreduce gx
      for (PetscInt i = 0; i < m; ++i) {
        tmp = gx[i];
        gx[i] = 0.0;
        MPI_Allreduce(&tmp, &(gx[i]), 1, MPIU_SCALAR, MPI_SUM,
            PETSC_COMM_WORLD);
        gx[i] = gx[i]
                / ((PetscScalar) neltot - nNonDesign)
                - volfrac; // # modified
      }
    }
    VecRestoreArray (xPhys, &xp);
    VecRestoreArray (xPassive0, &xPassive0p); // # new
    VecRestoreArray (xPassive1, &xPassive1p); // # new
    VecRestoreArray (xPassive2, &xPassive2p); // # new
    VecRestoreArray (xPassive3, &xPassive3p); // # new
    VecRestoreArray (Uloc, &up);
    VecDestroy (&Uloc);
  } //# new

  return (ierr);
}

PetscErrorCode
LinearElasticity::ComputeSensitivities (Vec dfdx, Vec *dgdx,
    Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar penal,
    PetscScalar volfrac, Vec xPassive0, Vec xPassive1,
    Vec xPassive2, Vec xPassive3) {

  PetscErrorCode ierr;

  // Get the FE mesh structure (from the nodal mesh)
  PetscInt nel, nen;
  const PetscInt *necon;
#if DIM == 2  // # new
  ierr = DMDAGetElements_2D (da_nodal, &nel, &nen, &necon);
#elif DIM == 3
  ierr = DMDAGetElements_3D (da_nodal, &nel, &nen, &necon);
#endif
  CHKERRQ(ierr);

  // Get pointer to the densities
  PetscScalar *xp, *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p; // # modified
  VecGetArray (xPhys, &xp);
  VecGetArray (xPassive0, &xPassive0p); // # new
  VecGetArray (xPassive1, &xPassive1p); // # new
  VecGetArray (xPassive2, &xPassive2p); // # new
  VecGetArray (xPassive3, &xPassive3p); // # new

  // Get Solution
  Vec Uloc;
  DMCreateLocalVector (da_nodal, &Uloc);
  DMGlobalToLocalBegin (da_nodal, U, INSERT_VALUES, Uloc);
  DMGlobalToLocalEnd (da_nodal, U, INSERT_VALUES, Uloc);

  // get pointer to local vector
  PetscScalar *up;
  VecGetArray (Uloc, &up);

  // Get dfdx
  PetscScalar *df;
  VecGetArray (dfdx, &df);

  // # new; Get dgdx
  // # new; Get dgdx
  PetscScalar **dg;
  for (PetscInt i = 0; i < m; ++i) {
    VecSet (dgdx[i], 0);
  }
  VecGetArrays (dgdx, m, &dg);

  // Number of total elements and nonDesign domain elements
  PetscInt neltot = 0;
  PetscScalar nNonDesign = 0; // # new
  VecGetSize (xPhys, &neltot); // # modified

  // Edof array
  PetscInt edof[nedof]; // # modified

  // # modified; Loop over elements
  for (PetscInt i = 0; i < nel; i++) {
    // loop over element nodes
    if (xPassive0p[i] == 0 && xPassive1p[i] == 0 && xPassive2p[i] == 0) {
      for (PetscInt j = 0; j < nen; j++) {
        // Get local dofs
        for (PetscInt k = 0; k < DIM; k++) {
          edof[j * DIM + k] = DIM * necon[i * nen + j] + k;
        }
      }
      // # modified; Use SIMP for stiffness interpolation
      PetscScalar uKu = 0.0;
      for (PetscInt k = 0; k < nedof; k++) {
        for (PetscInt h = 0; h < nedof; h++) {
          uKu += up[edof[k]] * KE[k * nedof + h] * up[edof[h]];
        }
      }
      // Set the Senstivity
      df[i] = -1.0 * penal * PetscPowScalar(xp[i], penal - 1) * (Emax - Emin)
              * uKu;
      // # new; Constraints
      for (PetscInt j = 0; j < m; ++j) {
        dg[j][i] = 1;
      }
    } else if (xPassive0p[i] == 1) { // # new
      df[i] = 1.0E9; // # new
      nNonDesign += 1; // # new
    } else if (xPassive1p[i] >= 1 || xPassive2p[i] >= 1
               || xPassive3p[i] == 1) { // # new
      df[i] = -1.0E9; // # new
      nNonDesign += 1; // # new
    }

    // Allreduce nNonDesign
    PetscScalar tmp = nNonDesign;
    nNonDesign = 0.0;
    MPI_Allreduce(&tmp, &(nNonDesign), 1, MPIU_SCALAR, MPI_SUM,
        PETSC_COMM_WORLD);

    // # modified; Allreduce gx
    for (PetscInt i = 0; i < m; ++i) {
      VecScale (dgdx[i],
          1.0 / ((PetscScalar) neltot - nNonDesign)); // # modified
    }

    VecRestoreArray (xPhys, &xp);
    VecRestoreArray (xPassive0, &xPassive0p); // # new
    VecRestoreArray (xPassive1, &xPassive1p); // # new
    VecRestoreArray (xPassive2, &xPassive2p); // # new
    VecRestoreArray (xPassive3, &xPassive3p); // # new
    VecRestoreArray (Uloc, &up);
    VecRestoreArray (dfdx, &df);
    VecRestoreArrays (dgdx, m, &dg);
    VecDestroy (&Uloc);

  } // # new

  return (ierr);
}

PetscErrorCode
LinearElasticity::WriteRestartFiles ()
{

  PetscErrorCode ierr = 0;

// Only dump data if correct allocater has been used
  if (!restart) {
    return -1;
  }

// Choose previous set of restart files
  if (flip) {
    flip = PETSC_FALSE;
  } else {
    flip = PETSC_TRUE;
  }

// Open viewers for writing
  PetscViewer view; // vectors
  if (!flip) {
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filename00.c_str (),
        FILE_MODE_WRITE, &view);
  } else if (flip) {
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filename01.c_str (),
        FILE_MODE_WRITE, &view);
  }

// Write vectors
  VecView (U, view);

// Clean up
  PetscViewerDestroy (&view);

  return ierr;
}

//##################################################################
//##################################################################
//##################################################################
// ######################## PRIVATE ################################
//##################################################################
//##################################################################

PetscErrorCode
LinearElasticity::AssembleStiffnessMatrix (Vec xPhys,
    PetscScalar Emin, PetscScalar Emax, PetscScalar penal,
    PetscInt loadCondition) {

  PetscErrorCode ierr;

// Get the FE mesh structure (from the nodal mesh)
  PetscInt nel, nen;
  const PetscInt *necon;
#if DIM == 2    // # new
  ierr = DMDAGetElements_2D (da_nodal, &nel, &nen, &necon);
  CHKERRQ(ierr);
#elif DIM == 3
  ierr = DMDAGetElements_3D (da_nodal, &nel, &nen, &necon);
  CHKERRQ(ierr);
#endif

// Get pointer to the densities
  PetscScalar *xp;
  VecGetArray (xPhys, &xp);

// Zero the matrix
  MatZeroEntries (K);

// # modified; Edof array
  PetscInt edof[nedof];
  PetscScalar ke[nedof * nedof];

// # modified; Loop over elements
  for (PetscInt i = 0; i < nel; i++) {
    // loop over element nodes
    for (PetscInt j = 0; j < nen; j++) {
      // Get local dofs
      for (PetscInt k = 0; k < DIM; k++) {
        edof[j * DIM + k] = DIM * necon[i * nen + j] + k;
      }
    }
    // Use SIMP for stiffness interpolation
    PetscScalar dens = Emin + PetscPowScalar(xp[i], penal) * (Emax - Emin);
    for (PetscInt k = 0; k < nedof * nedof; k++) {
      ke[k] = KE[k] * dens;
    }
    // Add values to the sparse matrix
    ierr = MatSetValuesLocal (K, nedof, edof, nedof, edof, ke, ADD_VALUES);
    CHKERRQ(ierr);
  }
  MatAssemblyBegin (K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (K, MAT_FINAL_ASSEMBLY);

// Impose the dirichlet conditions, i.e. K = N'*K*N - (N-I)
// 1.: K = N'*K*N
  MatDiagonalScale (K, N[loadCondition], N[loadCondition]); // # modified
// 2. Add ones, i.e. K = K + NI, NI = I - N
  Vec NI;
  VecDuplicate (N[loadCondition], &NI); // # modified
  VecSet (NI, 1.0);
  VecAXPY (NI, -1.0, N[loadCondition]); // # modified
  MatDiagonalSet (K, NI, ADD_VALUES);

// Zero out possible loads in the RHS that coincide
// with Dirichlet conditions
  VecPointwiseMult (RHS[loadCondition], RHS[loadCondition], N[loadCondition]); // # modified

  VecDestroy (&NI);
  VecRestoreArray (xPhys, &xp);
  DMDARestoreElements (da_nodal, &nel, &nen, &necon);

  return ierr;
}

PetscErrorCode
LinearElasticity::SetUpSolver ()
{

  PetscErrorCode ierr;

// CHECK FOR RESTART POINT
  restart = PETSC_TRUE;
  flip = PETSC_TRUE;
  PetscBool flg, onlyDesign;
  onlyDesign = PETSC_FALSE;
  char filenameChar[PETSC_MAX_PATH_LEN];
  PetscOptionsGetBool (NULL, NULL, "-restart", &restart, &flg);
  PetscOptionsGetBool (NULL, NULL, "-onlyLoadDesign", &onlyDesign, &flg); // DONT READ DESIGN IF THIS IS TRUE

// READ THE RESTART FILE INTO THE SOLUTION VECTOR(S)
  if (restart) {
    // THE FILES FOR WRITING RESTARTS
    std::string filenameWorkdir = "./";
    PetscOptionsGetString (NULL, NULL, "-workdir", filenameChar,
        sizeof(filenameChar), &flg);
    if (flg) {
      filenameWorkdir = "";
      filenameWorkdir.append (filenameChar);
    }
    filename00 = filenameWorkdir;
    filename01 = filenameWorkdir;
    filename00.append ("/RestartSol00.dat");
    filename01.append ("/RestartSol01.dat");

    // CHECK FOR SOLUTION AND READ TO STATE VECTOR(s)
    if (!onlyDesign) {
      // Where to read the restart point from
      std::string restartFileVec = ""; // NO RESTART FILE !!!!!
      // GET FILENAME
      PetscOptionsGetString (NULL, NULL, "-restartFileVecSol", filenameChar,
          sizeof(filenameChar), &flg);
      if (flg) {
        restartFileVec.append (filenameChar);
      }

      // PRINT TO SCREEN
      PetscPrintf (PETSC_COMM_WORLD,
          "# Restarting with solution (State Vector) from "
              "(-restartFileVecSol): %s \n", restartFileVec.c_str ());

      // Check if files exist:
      PetscBool vecFile = fexists (restartFileVec);
      if (!vecFile) {
        PetscPrintf (PETSC_COMM_WORLD, "File: %s NOT FOUND \n",
            restartFileVec.c_str ());
      }

      // READ
      if (vecFile) {
        PetscViewer view;
        // Open the data files
        ierr = PetscViewerBinaryOpen (PETSC_COMM_WORLD,
            restartFileVec.c_str (),
            FILE_MODE_READ, &view);

        VecLoad (U, view);

        PetscViewerDestroy (&view);
      }
    }
  }

  PC pc;

// The fine grid Krylov method
  KSPCreate (PETSC_COMM_WORLD, &(ksp));

// SET THE DEFAULT SOLVER PARAMETERS
// The fine grid solver settings
  PetscScalar rtol = 1.0e-5;
  PetscScalar atol = 1.0e-50;
  PetscScalar dtol = 1.0e5;
  PetscInt restart = 100;
  PetscInt maxitsGlobal = 200;

// Coarsegrid solver
  PetscScalar coarse_rtol = 1.0e-8;
  PetscScalar coarse_atol = 1.0e-50;
  PetscScalar coarse_dtol = 1e5;
  PetscInt coarse_maxits = 30;
  PetscInt coarse_restart = 30;

// Number of smoothening iterations per up/down smooth_sweeps
  PetscInt smooth_sweeps = 4;

// Set up the solver
  ierr = KSPSetType (ksp, KSPFGMRES); // KSPCG, KSPGMRES
  CHKERRQ(ierr);

  ierr = KSPGMRESSetRestart (ksp, restart);
  CHKERRQ(ierr);

  ierr = KSPSetTolerances (ksp, rtol, atol, dtol, maxitsGlobal);
  CHKERRQ(ierr);

  ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = KSPSetOperators (ksp, K, K);
  CHKERRQ(ierr);

// The preconditinoer
  KSPGetPC (ksp, &pc);
// Make PCMG the default solver
  PCSetType (pc, PCMG);

// Set solver from options
  KSPSetFromOptions (ksp);

// Get the prec again - check if it has changed
  KSPGetPC (ksp, &pc);

// Flag for pcmg pc
  PetscBool pcmg_flag = PETSC_TRUE;
  PetscObjectTypeCompare ((PetscObject) pc, PCMG, &pcmg_flag);

// Only if PCMG is used
  if (pcmg_flag) {

    // DMs for grid hierachy
    DM *da_list, *daclist;
    Mat R;

    PetscMalloc(sizeof(DM) * nlvls, &da_list);
    for (PetscInt k = 0; k < nlvls; k++)
      da_list[k] = NULL;
    PetscMalloc(sizeof(DM) * nlvls, &daclist);
    for (PetscInt k = 0; k < nlvls; k++)
      daclist[k] = NULL;

    // Set 0 to the finest level
    daclist[0] = da_nodal;

// Coordinates
#if DIM == 2  // # new
    PetscReal xmin = xc[0], xmax = xc[1], ymin = xc[2], ymax = xc[3];
#elif DIM == 3
    PetscReal xmin = xc[0], xmax = xc[1], ymin = xc[2], ymax = xc[3],
        zmin = xc[4], zmax = xc[5];
#endif
    // Set up the coarse meshes
    DMCoarsenHierarchy (da_nodal, nlvls - 1, &daclist[1]);
    for (PetscInt k = 0; k < nlvls; k++) {
      // NOTE: finest grid is nlevels - 1: PCMG MUST USE THIS ORDER ???
      da_list[k] = daclist[nlvls - 1 - k];
      // THIS SHOULD NOT BE NECESSARY
#if DIM == 2   // # new
      DMDASetUniformCoordinates (da_list[k], xmin, xmax, ymin, ymax, 0.0, 0.0);
#elif DIM == 3
      DMDASetUniformCoordinates (da_list[k], xmin, xmax, ymin, ymax, zmin,
          zmax);
#endif
    }
    // the PCMG specific options
    PCMGSetLevels (pc, nlvls, NULL);
    PCMGSetType (pc, PC_MG_MULTIPLICATIVE); // Default
    ierr = PCMGSetCycleType (pc, PC_MG_CYCLE_V);
    CHKERRQ(ierr);
    PCMGSetGalerkin (pc, PC_MG_GALERKIN_BOTH);
    for (PetscInt k = 1; k < nlvls; k++) {
      DMCreateInterpolation (da_list[k - 1], da_list[k], &R, NULL);
      PCMGSetInterpolation (pc, k, R);
      MatDestroy (&R);
    }

    // tidy up
    for (PetscInt k = 1; k < nlvls; k++) { // DO NOT DESTROY LEVEL 0
      DMDestroy (&daclist[k]);
    }
    PetscFree(da_list);
    PetscFree(daclist);

    // AVOID THE DEFAULT FOR THE MG PART
    {
      // SET the coarse grid solver:
      // i.e. get a pointer to the ksp and change its settings
      KSP cksp;
      PCMGGetCoarseSolve (pc, &cksp);
      // The solver
      ierr = KSPSetType (cksp, KSPGMRES); // KSPCG, KSPFGMRES
      ierr = KSPGMRESSetRestart (cksp, coarse_restart);
      // ierr = KSPSetType(cksp,KSPCG);

      ierr = KSPSetTolerances (cksp, coarse_rtol, coarse_atol, coarse_dtol,
          coarse_maxits);
      // The preconditioner
      PC cpc;
      KSPGetPC (cksp, &cpc);
      PCSetType (cpc, PCSOR); // PCGAMG, PCSOR, PCSPAI (NEEDS TO BE COMPILED), PCJACOBI

      // Set smoothers on all levels (except for coarse grid):
      for (PetscInt k = 1; k < nlvls; k++) {
        KSP dksp;
        PCMGGetSmoother (pc, k, &dksp);
        PC dpc;
        KSPGetPC (dksp, &dpc);
        ierr = KSPSetType (dksp,
        KSPGMRES); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
        ierr = KSPGMRESSetRestart (dksp, smooth_sweeps);
        // ierr = KSPSetType(dksp,KSPCHEBYSHEV);
        ierr = KSPSetTolerances (dksp, PETSC_DEFAULT, PETSC_DEFAULT,
        PETSC_DEFAULT, smooth_sweeps); // NOTE in the above maxitr=restart;
        PCSetType (dpc, PCSOR); // PCJACOBI, PCSOR for KSPCHEBYSHEV very good
      }
    }

// # new; The bleow commented code is about using the Cholesky direct solver
// for the coarse grid
//    // AVOID THE DEFAULT FOR THE MG PART
//    {
//      // SET the coarse grid solver:
//      // i.e. get a pointer to the ksp and change its settings
//      KSP cksp;
//      PCMGGetCoarseSolve (pc, &cksp);
//      // The solver
//      ierr = KSPSetType (cksp, KSPPREONLY); // KSPCG, KSPFGMRES
//      ierr = KSPGMRESSetRestart (cksp, coarse_restart);
//      // ierr = KSPSetType(cksp,KSPCG);
//
//      ierr = KSPSetTolerances (cksp, coarse_rtol, coarse_atol, coarse_dtol,
//          coarse_maxits);
//      // The preconditioner
//      PC cpc;
//      KSPGetPC (cksp, &cpc);
//      PCSetType (cpc, PCCHOLESKY); // PCGAMG, PCSOR, PCSPAI (NEEDS TO BE COMPILED), PCJACOBI
//
//      // Set smoothers on all levels (except for coarse grid):
//      for (PetscInt k = 1; k < nlvls; k++) {
//        KSP dksp;
//        PCMGGetSmoother (pc, k, &dksp);
//        PC dpc;
//        KSPGetPC (dksp, &dpc);
//        ierr = KSPSetType (dksp,
//        KSPCG); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
//        ierr = KSPGMRESSetRestart (dksp, smooth_sweeps);
//        ierr = KSPSetTolerances (dksp, PETSC_DEFAULT, PETSC_DEFAULT,
//        PETSC_DEFAULT, smooth_sweeps); // NOTE in the above maxitr=restart;
//        PCSetType (dpc, PCJACOBI); // PCJACOBI, PCSOR for KSPCHEBYSHEV very good
//      }
//    }
  }

// Write check to screen:
// Check the overall Krylov solver
  KSPType ksptype;
  KSPGetType (ksp, &ksptype);
  PCType pctype;
  PCGetType (pc, &pctype);
  PetscInt mmax;
  KSPGetTolerances (ksp, NULL, NULL, NULL, &mmax);
  PetscPrintf (PETSC_COMM_WORLD,
      "##############################################################\n");
  PetscPrintf (PETSC_COMM_WORLD,
      "################# Linear solver settings #####################\n");
  PetscPrintf (PETSC_COMM_WORLD,
      "# Main solver: %s, prec.: %s, maxiter.: %i \n", ksptype, pctype, mmax);

// Only if pcmg is used
  if (pcmg_flag) {
    // Check the smoothers and coarse grid solver:
    for (PetscInt k = 0; k < nlvls; k++) {
      KSP dksp;
      PC dpc;
      KSPType dksptype;
      PCMGGetSmoother (pc, k, &dksp);
      KSPGetType (dksp, &dksptype);
      KSPGetPC (dksp, &dpc);
      PCType dpctype;
      PCGetType (dpc, &dpctype);
      PetscInt mmax;
      KSPGetTolerances (dksp, NULL, NULL, NULL, &mmax);
      PetscPrintf (PETSC_COMM_WORLD,
          "# Level %i smoother: %s, prec.: %s, sweep: %i \n", k, dksptype,
          dpctype, mmax);
    }
  }
  PetscPrintf (PETSC_COMM_WORLD,
      "##############################################################\n");

  return (ierr);
}

#if DIM == 2   // # new
PetscErrorCode LinearElasticity::DMDAGetElements_2D (DM dm, PetscInt *nel,
    PetscInt *nen, const PetscInt *e[]) {
  PetscErrorCode ierr;
  DM_DA *da = (DM_DA*) dm->data;
  PetscInt i, xs, xe, Xs, Xe;
  PetscInt j, ys, ye, Ys, Ye;
  PetscInt cnt = 0, cell[4], ns = 1, nn = 4;
  PetscInt c;
  if (!da->e) {
    if (da->elementtype == DMDA_ELEMENT_Q1) {
      ns = 1;
      nn = 4;
    }
    ierr = DMDAGetCorners (dm, &xs, &ys, NULL, &xe, &ye, NULL);
    CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners (dm, &Xs, &Ys, NULL, &Xe, &Ye, NULL);
    CHKERRQ(ierr);
    xe += xs;
    Xe += Xs;
    if (xs != Xs) xs -= 1;
    ye += ys;
    Ye += Ys;
    if (ys != Ys) ys -= 1;

    da->ne = ns * (xe - xs - 1) * (ye - ys - 1);
    PetscMalloc((1 + nn * da->ne) * sizeof(PetscInt), &da->e);
    for (j = ys; j < ye - 1; j++) {
      for (i = xs; i < xe - 1; i++) {
        cell[0] = (i - Xs) + (j - Ys) * (Xe - Xs);
        cell[1] = (i - Xs + 1) + (j - Ys) * (Xe - Xs);
        cell[2] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs);
        cell[3] = (i - Xs) + (j - Ys + 1) * (Xe - Xs);
        if (da->elementtype == DMDA_ELEMENT_Q1) {
          for (c = 0; c < ns * nn; c++)
            da->e[cnt++] = cell[c];
        }
      }
    }
  }
  *nel = da->ne;
  *nen = nn;
  *e = da->e;
  return (0);
}

PetscInt LinearElasticity::Quad4Isoparametric (PetscScalar *X, PetscScalar *Y,
    PetscScalar nu, PetscInt redInt, PetscScalar *ke) {
  // QUA4_ISOPARAMETRIC - Computes QUA4 isoparametric element matrices
  // The element stiffness matrix is computed as:
  //
  //       ke = int(int(B^T*C*B,x),y)
  //
  // For an isoparameteric element this integral becomes:
  //
  //       ke = int(int(B^T*C*B*det(J),xi=-1..1),eta=-1..1)
  //
  // where B is the more complicated expression:
  // B = [dx*alpha1 + dy*alpha2]*N
  // where
  // dx = [invJ11 invJ12]*[dxi deta]
  // dy = [invJ21 invJ22]*[dxi deta]
  //
  // Remark: The elasticity modulus is left out in the below
  // computations, because we multiply with it afterwards (the aim is
  // topology optimization).
  // Furthermore, this is not the most efficient code, but it is readable.
  //
  /////////////////////////////////////////////////////////////////////////////////
  //////// INPUT:
  // X, Y = Vectors containing the coordinates of the four nodes
  //               (x1,y1,x2,y2,x3,y3,x4,y4). Where node 1 is in the
  //               lower left corner, and node 2 is the next node
  //               counterclockwise.
  // redInt   = Reduced integration option boolean (here an integer).
  //                  redInt == 0 (false): Full integration
  //                  redInt == 1 (true): Reduced integration
  // nu               = Poisson's ratio.
  //
  //////// OUTPUT:
  // ke  = Element stiffness matrix. Needs to be multiplied with elasticity
  // modulus
  //
  //   Written 2013 at
  //   Department of Mechanical Engineering
  //   Technical University of Denmark (DTU).
  //
  //   Modified 2020 by Zhidong Brian Zhang at
  //   Multi-scale additive manufacturing lab (MSAM)
  //   University of Waterloo
  /////////////////////////////////////////////////////////////////////////////////

  //// COMPUTE ELEMENT STIFFNESS MATRIX
  // Lame's parameters (with E=1.0):
  PetscScalar E = this->E; // # modified
  // Constitutive matrix, plane stress
  PetscScalar C[3][3] = { { E / (1.0 - nu * nu), E * nu / (1.0 - nu * nu), 0 },
                          { E * nu / (1.0 - nu * nu), E / (1.0 - nu * nu), 0 },
                          { 0, 0, E / (2.0 * (1.0 + nu)) } };
//  // Constitutive matrix, plane strain
//  PetscScalar C[3][3] = { { (1.0 - nu), nu, 0 },
//                          { nu, (1.0 - nu), 0 },
//                          { 0, 0, ((1.0 - 2.0 * nu) / 2.0) } };
//  for (PetscInt i = 0; i < 3 * 3; ++i) {
//    C[i / 3][i % 3] *= E / (1.0 + nu) / (1.0 - 2.0 * nu);
//  }

  // Gauss points (GP) and weigths
  // Two Gauss points in all directions (total of four)
  PetscScalar GP[2] = { -0.577350269189626, 0.577350269189626 };
  // Corresponding weights
  PetscScalar W[2] = { 1.0, 1.0 };
  // If reduced integration only use one GP
  if (redInt) {
    GP[0] = 0.0;
    W[0] = 2.0;
  }
  // Matrices that help when we gather the strain-displacement matrix:
  PetscScalar alpha1[3][2];
  PetscScalar alpha2[3][2];
  memset (alpha1, 0, sizeof(alpha1[0][0]) * 3 * 2); // zero out
  memset (alpha2, 0, sizeof(alpha2[0][0]) * 3 * 2); // zero out
  alpha1[0][0] = 1.0;
  alpha1[2][1] = 1.0;
  alpha2[1][1] = 1.0;
  alpha2[2][0] = 1.0;
  PetscScalar dNdxi[4];
  PetscScalar dNdeta[4];
  PetscScalar J[2][2];
  PetscScalar invJ[2][2];
  PetscScalar beta[3][2];
  PetscScalar B[3][8]; // Note: Small enough to be allocated on stack
  PetscScalar *dN;
  // Make sure the stiffness matrix is zeroed out:
  memset (ke, 0, sizeof(ke[0]) * 8 * 8);
  // Perform the numerical integration
  for (PetscInt ii = 0; ii < 2 - redInt; ii++) {
    for (PetscInt jj = 0; jj < 2 - redInt; jj++) {
      // Integration point
      PetscScalar xi = GP[ii];
      PetscScalar eta = GP[jj];
      // Differentiated shape functions
      DifferentiatedShapeFunctions_2D (xi, eta, dNdxi, dNdeta);
      // Jacobian
      J[0][0] = Dot (dNdxi, X, 4);
      J[0][1] = Dot (dNdxi, Y, 4);
      J[1][0] = Dot (dNdeta, X, 4);
      J[1][1] = Dot (dNdeta, Y, 4);
      // Inverse and determinant
      PetscScalar detJ = Inverse2M (J, invJ);
      // Weight factor at this point
      PetscScalar weight = W[ii] * W[jj] * detJ;
      // Strain-displacement matrix
      memset (B, 0, sizeof(B[0][0]) * 3 * 8); // zero out
      for (PetscInt ll = 0; ll < DIM; ll++) {
        // Add contributions from the different derivatives
        if (ll == 0) {
          dN = dNdxi;
        }
        if (ll == 1) {
          dN = dNdeta;
        }

        // Assemble strain operator
        for (PetscInt i = 0; i < 3; i++) {
          for (PetscInt j = 0; j < 2; j++) {
            beta[i][j] = invJ[0][ll] * alpha1[i][j]
                         + invJ[1][ll] * alpha2[i][j];
          }
        }
        // Add contributions to strain-displacement matrix
        for (PetscInt i = 0; i < 3; i++) {
          for (PetscInt j = 0; j < 8; j++) {
            B[i][j] = B[i][j] + beta[i][j % 2] * dN[j / 2];
          }
        }
      }
      // Finally, add to the element matrix
      for (PetscInt i = 0; i < 8; i++) {
        for (PetscInt j = 0; j < 8; j++) {
          for (PetscInt k = 0; k < 3; k++) {
            for (PetscInt l = 0; l < 3; l++) {

              ke[j + 8 * i] = ke[j + 8 * i]
                              + weight * (B[k][i] * C[k][l] * B[l][j]);
            }
          }
        }
      }
    }
  }
  return 0;
}

void LinearElasticity::DifferentiatedShapeFunctions_2D (PetscScalar xi,
    PetscScalar eta, PetscScalar *dNdxi, PetscScalar *dNdeta) {
  // differentiatedShapeFunctions - Computes differentiated shape functions
  // At the point given by (xi, eta, zeta).
  // With respect to xi:
  dNdxi[0] = -0.25 * (1.0 - eta);
  dNdxi[1] = 0.25 * (1.0 - eta);
  dNdxi[2] = 0.25 * (1.0 + eta);
  dNdxi[3] = -0.25 * (1.0 + eta);
  // With respect to eta:
  dNdeta[0] = -0.25 * (1.0 - xi);
  dNdeta[1] = -0.25 * (1.0 + xi);
  dNdeta[2] = 0.25 * (1.0 + xi);
  dNdeta[3] = 0.25 * (1.0 - xi);
}

PetscScalar LinearElasticity::Inverse2M (PetscScalar J[][2],
    PetscScalar invJ[][2]) {
  // inverse3M - Computes the inverse of a 3x3 matrix
  PetscScalar detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  invJ[0][0] = J[1][1] / detJ;
  invJ[0][1] = -J[0][1] / detJ;
  invJ[1][0] = -J[1][0] / detJ;
  invJ[1][1] = J[0][0] / detJ;
  return detJ;
}

#elif DIM == 3
PetscErrorCode
LinearElasticity::DMDAGetElements_3D (DM dm, PetscInt *nel,
    PetscInt *nen, const PetscInt *e[]) {
  PetscErrorCode ierr;
  DM_DA *da = (DM_DA*) dm->data;
  PetscInt i, xs, xe, Xs, Xe;
  PetscInt j, ys, ye, Ys, Ye;
  PetscInt k, zs, ze, Zs, Ze;
  PetscInt cnt = 0, cell[8], ns = 1, nn = 8;
  PetscInt c;
  if (!da->e) {
    if (da->elementtype == DMDA_ELEMENT_Q1) {
      ns = 1;
      nn = 8;
    }
    ierr = DMDAGetCorners (dm, &xs, &ys, &zs, &xe, &ye, &ze);
    CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners (dm, &Xs, &Ys, &Zs, &Xe, &Ye, &Ze);
    CHKERRQ(ierr);
    xe += xs;
    Xe += Xs;
    if (xs != Xs) xs -= 1;
    ye += ys;
    Ye += Ys;
    if (ys != Ys) ys -= 1;
    ze += zs;
    Ze += Zs;
    if (zs != Zs) zs -= 1;
    da->ne = ns * (xe - xs - 1) * (ye - ys - 1) * (ze - zs - 1);
    PetscMalloc((1 + nn * da->ne) * sizeof(PetscInt), &da->e);
    for (k = zs; k < ze - 1; k++) {
      for (j = ys; j < ye - 1; j++) {
        for (i = xs; i < xe - 1; i++) {
          cell[0] = (i - Xs) + (j - Ys) * (Xe - Xs)
                    + (k - Zs) * (Xe - Xs) * (Ye - Ys);
          cell[1] = (i - Xs + 1) + (j - Ys) * (Xe - Xs)
                    + (k - Zs) * (Xe - Xs) * (Ye - Ys);
          cell[2] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs)
                    + (k - Zs) * (Xe - Xs) * (Ye - Ys);
          cell[3] = (i - Xs) + (j - Ys + 1) * (Xe - Xs)
                    + (k - Zs) * (Xe - Xs) * (Ye - Ys);
          cell[4] = (i - Xs) + (j - Ys) * (Xe - Xs)
                    + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
          cell[5] = (i - Xs + 1) + (j - Ys) * (Xe - Xs)
                    + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
          cell[6] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs)
                    + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
          cell[7] = (i - Xs) + (j - Ys + 1) * (Xe - Xs)
                    + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
          if (da->elementtype == DMDA_ELEMENT_Q1) {
            for (c = 0; c < ns * nn; c++)
              da->e[cnt++] = cell[c];
          }
        }
      }
    }
  }
  *nel = da->ne;
  *nen = nn;
  *e = da->e;
  return (0);
}

PetscInt
LinearElasticity::Hex8Isoparametric (PetscScalar *X, PetscScalar *Y,
    PetscScalar *Z, PetscScalar nu, PetscInt redInt, PetscScalar *ke) {
// HEX8_ISOPARAMETRIC - Computes HEX8 isoparametric element matrices
// The element stiffness matrix is computed as:
//
//       ke = int(int(int(B^T*C*B,x),y),z)
//
// For an isoparameteric element this integral becomes:
//
//       ke = int(int(int(B^T*C*B*det(J),xi=-1..1),eta=-1..1),zeta=-1..1)
//
// where B is the more complicated expression:
// B = [dx*alpha1 + dy*alpha2 + dz*alpha3]*N
// where
// dx = [invJ11 invJ12 invJ13]*[dxi deta dzeta]
// dy = [invJ21 invJ22 invJ23]*[dxi deta dzeta]
// dy = [invJ31 invJ32 invJ33]*[dxi deta dzeta]
//
// Remark: The elasticity modulus is left out in the below
// computations, because we multiply with it afterwards (the aim is
// topology optimization).
// Furthermore, this is not the most efficient code, but it is readable.
//
/////////////////////////////////////////////////////////////////////////////////
//////// INPUT:
// X, Y, Z  = Vectors containing the coordinates of the eight nodes
//               (x1,y1,z1,x2,y2,z2,...,x8,y8,z8). Where node 1 is in the
//               lower left corner, and node 2 is the next node
//               counterclockwise (looking in the negative z-dir). Finish the
//               x-y-plane and then move in the positive z-dir.
// redInt   = Reduced integration option boolean (here an integer).
//           	redInt == 0 (false): Full integration
//           	redInt == 1 (true): Reduced integration
// nu 		= Poisson's ratio.
//
//////// OUTPUT:
// ke  = Element stiffness matrix. Needs to be multiplied with elasticity
// modulus
//
//   Written 2013 at
//   Department of Mechanical Engineering
//   Technical University of Denmark (DTU).
/////////////////////////////////////////////////////////////////////////////////

//// COMPUTE ELEMENT STIFFNESS MATRIX
// Lame's parameters (with E=1.0):
  PetscScalar E = this->E;  // # new
  PetscScalar lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // # modified
  PetscScalar mu = E / (2.0 * (1.0 + nu));  // # modified
// Constitutive matrix
  PetscScalar C[6][6] = { { lambda + 2.0 * mu, lambda, lambda, 0.0, 0.0, 0.0 }, { lambda, lambda
      + 2.0 * mu, lambda, 0.0, 0.0, 0.0 }, { lambda, lambda, lambda + 2.0 * mu, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, mu, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, mu, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, mu } };
// Gauss points (GP) and weigths
// Two Gauss points in all directions (total of eight)
  PetscScalar GP[2] = { -0.577350269189626, 0.577350269189626 };
// Corresponding weights
  PetscScalar W[2] = { 1.0, 1.0 };
// If reduced integration only use one GP
  if (redInt) {
    GP[0] = 0.0;
    W[0] = 2.0;
  }
// Matrices that help when we gather the strain-displacement matrix:
  PetscScalar alpha1[6][3];
  PetscScalar alpha2[6][3];
  PetscScalar alpha3[6][3];
  memset (alpha1, 0, sizeof(alpha1[0][0]) * 6 * 3); // zero out
  memset (alpha2, 0, sizeof(alpha2[0][0]) * 6 * 3); // zero out
  memset (alpha3, 0, sizeof(alpha3[0][0]) * 6 * 3); // zero out
  alpha1[0][0] = 1.0;
  alpha1[3][1] = 1.0;
  alpha1[5][2] = 1.0;
  alpha2[1][1] = 1.0;
  alpha2[3][0] = 1.0;
  alpha2[4][2] = 1.0;
  alpha3[2][2] = 1.0;
  alpha3[4][1] = 1.0;
  alpha3[5][0] = 1.0;
  PetscScalar dNdxi[8];
  PetscScalar dNdeta[8];
  PetscScalar dNdzeta[8];
  PetscScalar J[3][3];
  PetscScalar invJ[3][3];
  PetscScalar beta[6][3];
  PetscScalar B[6][24]; // Note: Small enough to be allocated on stack
  PetscScalar *dN;
// Make sure the stiffness matrix is zeroed out:
  memset (ke, 0, sizeof(ke[0]) * 24 * 24);
// Perform the numerical integration
  for (PetscInt ii = 0; ii < 2 - redInt; ii++) {
    for (PetscInt jj = 0; jj < 2 - redInt; jj++) {
      for (PetscInt kk = 0; kk < 2 - redInt; kk++) {
        // Integration point
        PetscScalar xi = GP[ii];
        PetscScalar eta = GP[jj];
        PetscScalar zeta = GP[kk];
        // Differentiated shape functions
        DifferentiatedShapeFunctions (xi, eta, zeta, dNdxi, dNdeta, dNdzeta);
        // Jacobian
        J[0][0] = Dot (dNdxi, X, 8);
        J[0][1] = Dot (dNdxi, Y, 8);
        J[0][2] = Dot (dNdxi, Z, 8);
        J[1][0] = Dot (dNdeta, X, 8);
        J[1][1] = Dot (dNdeta, Y, 8);
        J[1][2] = Dot (dNdeta, Z, 8);
        J[2][0] = Dot (dNdzeta, X, 8);
        J[2][1] = Dot (dNdzeta, Y, 8);
        J[2][2] = Dot (dNdzeta, Z, 8);
        // Inverse and determinant
        PetscScalar detJ = Inverse3M (J, invJ);
        // Weight factor at this point
        PetscScalar weight = W[ii] * W[jj] * W[kk] * detJ;
        // Strain-displacement matrix
        memset (B, 0, sizeof(B[0][0]) * 6 * 24); // zero out
        for (PetscInt ll = 0; ll < 3; ll++) {
          // Add contributions from the different derivatives
          if (ll == 0) {
            dN = dNdxi;
          }
          if (ll == 1) {
            dN = dNdeta;
          }
          if (ll == 2) {
            dN = dNdzeta;
          }
          // Assemble strain operator
          for (PetscInt i = 0; i < 6; i++) {
            for (PetscInt j = 0; j < 3; j++) {
              beta[i][j] = invJ[0][ll] * alpha1[i][j]
                           + invJ[1][ll] * alpha2[i][j]
                           + invJ[2][ll] * alpha3[i][j];
            }
          }
          // Add contributions to strain-displacement matrix
          for (PetscInt i = 0; i < 6; i++) {
            for (PetscInt j = 0; j < 24; j++) {
              B[i][j] = B[i][j] + beta[i][j % 3] * dN[j / 3];
            }
          }
        }
        // Finally, add to the element matrix
        for (PetscInt i = 0; i < 24; i++) {
          for (PetscInt j = 0; j < 24; j++) {
            for (PetscInt k = 0; k < 6; k++) {
              for (PetscInt l = 0; l < 6; l++) {

                ke[j + 24 * i] = ke[j + 24 * i]
                                 + weight * (B[k][i] * C[k][l] * B[l][j]);
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

void LinearElasticity::DifferentiatedShapeFunctions (PetscScalar xi,
    PetscScalar eta, PetscScalar zeta, PetscScalar *dNdxi, PetscScalar *dNdeta,
    PetscScalar *dNdzeta) {
// differentiatedShapeFunctions - Computes differentiated shape functions
// At the point given by (xi, eta, zeta).
// With respect to xi:
  dNdxi[0] = -0.125 * (1.0 - eta) * (1.0 - zeta);
  dNdxi[1] = 0.125 * (1.0 - eta) * (1.0 - zeta);
  dNdxi[2] = 0.125 * (1.0 + eta) * (1.0 - zeta);
  dNdxi[3] = -0.125 * (1.0 + eta) * (1.0 - zeta);
  dNdxi[4] = -0.125 * (1.0 - eta) * (1.0 + zeta);
  dNdxi[5] = 0.125 * (1.0 - eta) * (1.0 + zeta);
  dNdxi[6] = 0.125 * (1.0 + eta) * (1.0 + zeta);
  dNdxi[7] = -0.125 * (1.0 + eta) * (1.0 + zeta);
// With respect to eta:
  dNdeta[0] = -0.125 * (1.0 - xi) * (1.0 - zeta);
  dNdeta[1] = -0.125 * (1.0 + xi) * (1.0 - zeta);
  dNdeta[2] = 0.125 * (1.0 + xi) * (1.0 - zeta);
  dNdeta[3] = 0.125 * (1.0 - xi) * (1.0 - zeta);
  dNdeta[4] = -0.125 * (1.0 - xi) * (1.0 + zeta);
  dNdeta[5] = -0.125 * (1.0 + xi) * (1.0 + zeta);
  dNdeta[6] = 0.125 * (1.0 + xi) * (1.0 + zeta);
  dNdeta[7] = 0.125 * (1.0 - xi) * (1.0 + zeta);
// With respect to zeta:
  dNdzeta[0] = -0.125 * (1.0 - xi) * (1.0 - eta);
  dNdzeta[1] = -0.125 * (1.0 + xi) * (1.0 - eta);
  dNdzeta[2] = -0.125 * (1.0 + xi) * (1.0 + eta);
  dNdzeta[3] = -0.125 * (1.0 - xi) * (1.0 + eta);
  dNdzeta[4] = 0.125 * (1.0 - xi) * (1.0 - eta);
  dNdzeta[5] = 0.125 * (1.0 + xi) * (1.0 - eta);
  dNdzeta[6] = 0.125 * (1.0 + xi) * (1.0 + eta);
  dNdzeta[7] = 0.125 * (1.0 - xi) * (1.0 + eta);
}

PetscScalar
LinearElasticity::Inverse3M (PetscScalar J[][3],
    PetscScalar invJ[][3]) {
// inverse3M - Computes the inverse of a 3x3 matrix
  PetscScalar detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                     - J[0][1] * (J[1][0] * J[2][2] - J[2][0] * J[1][2])
                     + J[0][2] * (J[1][0] * J[2][1] - J[2][0] * J[1][1]);
  invJ[0][0] = (J[1][1] * J[2][2] - J[2][1] * J[1][2]) / detJ;
  invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) / detJ;
  invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / detJ;
  invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) / detJ;
  invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / detJ;
  invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) / detJ;
  invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / detJ;
  invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) / detJ;
  invJ[2][2] = (J[0][0] * J[1][1] - J[1][0] * J[0][1]) / detJ;
  return detJ;
}
#endif

PetscScalar
LinearElasticity::Dot (PetscScalar *v1, PetscScalar *v2,
    PetscInt l) { // # new
// Function that returns the dot product of v1 and v2,
// which must have the same length l
  PetscScalar result = 0.0;
  for (PetscInt i = 0; i < l; i++) {
    result = result + v1[i] * v2[i];
  }
  return result;
}
