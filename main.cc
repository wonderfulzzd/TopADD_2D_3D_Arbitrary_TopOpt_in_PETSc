#include "Filter.h"
#include "MMA.h"
#include "MPIIO.h"
#include "TopOpt.h"
#include "mpi.h"
#include <petsc.h>

#include "options.h" //zzd added
#include "timer.h" //zzd added

// Pre- and post-processing class
#include "PrePostProcess.h" //zzd added


// Choose the physical problem to be solved
#if PHYSICS == 0
#include "LinearElasticity.h"
#elif PHYSICS == 1
#include "LinearCompliant.h"   //zzd added
#elif PHYSICS == 2
#include "LinearHeatConduction.h" //zzd added
#endif

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

static char help[] = "2D/3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n"; // zzd

int main (int argc, char *argv[]) {

  // Error code for debugging
  PetscErrorCode ierr = 0;

  // Initialize PETSc / MPI and pass input arguments to PETSc
  PetscInitialize (&argc, &argv, PETSC_NULL, help);

  double ts, te; // time
  ts = MPI_Wtime ();

  // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
  TopOpt *opt = new TopOpt ();

  // STEP 2: Pre-processing to define the design domain by using passive element assigning method //zzd
  PrePostProcess *prepost = new PrePostProcess (opt); //zzd added
#if IMPORT_GEO == 1
  prepost->DesignDomainInitialization (opt);
#endif

  // STEP 3: THE PHYSICS
  // 0 - linear elasticity, 1 - linear heat conduction, 2 - compliant
#if PHYSICS == 0
  LinearElasticity *physics = new LinearElasticity (opt->da_nodes,
      opt->numLoads, opt->xPassive0, opt->xPassive1, opt->xPassive2);
#elif PHYSICS ==1
  LinearCompliant *physics = new LinearCompliant (opt->da_nodes, opt->numLoads,
      opt->xPassive0, opt->xPassive1, opt->xPassive2);
#elif PHYSICS == 2
  LinearHeatConduction *physics = new LinearHeatConduction (opt->da_nodes,
      opt->da_elem, opt->numLoads, opt->xPassive0, opt->xPassive1,
      opt->xPassive2);
#endif

  // STEP 4: THE FILTERING
  Filter *filter = new Filter (opt->da_nodes, opt->xPhys, opt->filter,
      opt->rmin, opt->xPassive0, opt->xPassive1, opt->xPassive2);

  // STEP 5: VISUALIZATION USING VTK
  MPIIO *output = new MPIIO (opt->da_nodes, 4, "ux, uy, uz, nodeDen", 6,
      "x, xTilde, xPhys, xPassive0, xPassive1, xPassive2"); // All point data must use 3 coordinates in VTK

  // STEP 6: THE OPTIMIZER MMA
  MMA *mma;
  PetscInt itr = 0;
  opt->AllocateMMAwithRestart (&itr, &mma); // allow for restart !
  // mma->SetAsymptotes(0.2, 0.65, 1.05);

  // STEP 7: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
  ierr = filter->FilterProject (opt->x, opt->xTilde, opt->xPhys,
      opt->projectionFilter, opt->beta, opt->eta);
  CHKERRQ(ierr);

  // STEP 8: OPTIMIZATION LOOP
  PetscScalar ch = 1.0;
  double t0, t1, t2, t3;
  t0 = MPI_Wtime ();

  while (itr < opt->maxItr && ch > 0.01) {
    // Update iteration counter
    itr++;

    // start timer
    t1 = MPI_Wtime ();
    t2 = MPI_Wtime ();

    // Compute (a) obj+const, (b) sens, (c) obj+const+sens
    ierr = physics->ComputeObjectiveConstraintsSensitivities (&(opt->fx),
        &(opt->gx[0]), opt->dfdx, opt->dgdx[0], opt->xPhys, opt->Emin,
        opt->Emax, opt->penal, opt->volfrac, opt->xPassive0, opt->xPassive1,
        opt->xPassive2);
    CHKERRQ(ierr); //zzd

    t3 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD, "tS: %f\n", t3 - t2);

    // Compute objective scale
    if (itr == 1) {
      opt->fscale = 10.0 / opt->fx; //zzd commented tmp
    }
    // Scale objectie and sens
    opt->fx = opt->fx * opt->fscale;
    VecScale (opt->dfdx, opt->fscale);

    t2 = MPI_Wtime ();

    // Filter sensitivities (chainrule)
    ierr = filter->Gradients (opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx,
        opt->projectionFilter, opt->beta, opt->eta);
    CHKERRQ(ierr);

    t3 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD, "tfilter: %f\n", t3 - t2);

    t2 = MPI_Wtime ();

    // Sets outer movelimits on design variables
    ierr = mma->SetOuterMovelimit (opt->Xmin, opt->Xmax, opt->movlim, opt->x,
        opt->xmin, opt->xmax);
    CHKERRQ(ierr);

    // Update design by MMA
    ierr = mma->Update (opt->x, opt->dfdx, opt->gx, opt->dgdx, opt->xmin,
        opt->xmax);
    CHKERRQ(ierr);

    // Inf norm on the design change
    ch = mma->DesignChange (opt->x, opt->xold);

    // Increase beta if needed
    PetscBool changeBeta = PETSC_FALSE;
    if (opt->projectionFilter) {
      changeBeta = filter->IncreaseBeta (&(opt->beta), opt->betaFinal,
          opt->gx[0], itr, ch);
    }

    t3 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD, "tOptization: %f\n", t3 - t2);

    // Filter design field
    ierr = filter->FilterProject (opt->x, opt->xTilde, opt->xPhys,
        opt->projectionFilter, opt->beta, opt->eta);
    CHKERRQ(ierr);

    // Discreteness measure
    PetscScalar mnd = filter->GetMND (opt->xPhys);

    // stop timer
    t2 = MPI_Wtime ();

    // Print to screen
    PetscPrintf (PETSC_COMM_WORLD,
        "It.: %i, True fx: %f, Scaled fx: %f, gx[0]: %f, ch.: %f, "
            "mnd.: %f, time: %f\n", itr, opt->fx / opt->fscale, opt->fx,
        opt->gx[0], ch, mnd, t2 - t1);

    // Write field data: first 10 iterations and then every 20th
    if (itr < 11 || itr % 20 == 0 || changeBeta) {
      prepost->UpdateNodeDensity (opt); // # new; Update node density
      output->WriteVTK (physics->da_nodal, physics->GetStateField (),
          opt->nodeDensity, opt->x, opt->xTilde, opt->xPhys, opt->xPassive0,
          opt->xPassive1, opt->xPassive2, itr); // # modified
    }

    // Dump data needed for restarting code at termination
    if (itr % 10 == 0) {
      opt->WriteRestartFiles (&itr, mma);
      physics->WriteRestartFiles ();
    }
  }
  // Write restart WriteRestartFiles
  opt->WriteRestartFiles (&itr, mma);
  physics->WriteRestartFiles ();

  // Dump final design
  output->WriteVTK (physics->da_nodal, physics->GetStateField (),
      opt->nodeDensity, opt->x, opt->xTilde, opt->xPhys, opt->xPassive0,
      opt->xPassive1, opt->xPassive2, itr); // # modified

  PetscPrintf (PETSC_COMM_WORLD, "tAll: %f\n", t3 - t0);

  // STEP 7: CLEAN UP AFTER YOURSELF
  delete mma;
  delete output;
  delete filter;
  delete opt;
  delete physics;
  delete prepost; // # new

  // Total computing time
  te = MPI_Wtime ();
  PetscPrintf (PETSC_COMM_WORLD, "Total computing time: %f min\n",
      (te - ts) / 60.0);

  // Finalize PETSc / MPI
  PetscFinalize ();
  return ierr;
}
