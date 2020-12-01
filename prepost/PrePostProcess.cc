// -------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the TopADD authors
//
// This file is part of the TopADD.
//
// The TopADD is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of TopADD.
//
// Author: Zhidong Brian Zhang
// Created on: Dec. 2018
// Updated on: Sep. 2020
//
// ---------------------------------------------------------------------

#include <PrePostProcess.h>

PrePostProcess::PrePostProcess (TopOpt *opt) {
  // Design domain dimensions
#if DIM == 2
  nx = opt->nxyz[0] - 1;
  ny = opt->nxyz[1] - 1;
  nz = 1;
  dx = opt->dx;
  dy = opt->dy;
  dz = opt->dx;
#elif DIM == 3
  nx = opt->nxyz[0] - 1;
  ny = opt->nxyz[1] - 1;
  nz = opt->nxyz[2] - 1;
  dx = opt->dx;
  dy = opt->dy;
  dz = opt->dz;
#endif
//  this->numLoads = 1; //opt->numLoads;
  this->numLoads = opt->numLoads;
  voxIndex = 0;
  occFIX.resize (numLoads);
  occLOD.resize (numLoads);
}

PrePostProcess::~PrePostProcess () {
}

PetscErrorCode PrePostProcess::DesignDomainInitialization (TopOpt *opt) {
  PetscErrorCode ierr = 0;

  double t1, t2;

  PetscPrintf (PETSC_COMM_WORLD,
      "################ Design domain initialization ################\n");

  // Import and voxelize
  t1 = MPI_Wtime ();
  ierr = ImportAndVoxelizeGeometry (opt);
  t2 = MPI_Wtime ();
  PetscPrintf (PETSC_COMM_WORLD,
      "# Importing and voxelizing totally took %f s\n", t2 - t1);

  // Assign passive element
  t1 = MPI_Wtime ();
  ierr = AssignPassiveElement (opt);
  t2 = MPI_Wtime ();
  PetscPrintf (PETSC_COMM_WORLD, "# Assigning passive element took %f s\n",
      t2 - t1);

  // Clean the occupancy data and free memory
  CleanUp ();

  // Calculate the node load adding total counts
  // This is for dividing the total force among all the loading nodes
  CalculateNodeLoadAddingTotalCounts (opt);

  PetscPrintf (PETSC_COMM_WORLD,
      "##############################################################\n");

  return ierr;
}

PetscErrorCode PrePostProcess::CalculateNodeLoadAddingTotalCounts (TopOpt *opt) {
  PetscErrorCode ierr = 0;

  // Get the FE mesh structure (from the nodal mesh)
  PetscInt nel, nen;
  const PetscInt *necon;
#if DIM == 2
  ierr = DMDAGetElements_2D (opt->da_nodes, &nel, &nen, &necon);
  CHKERRQ(ierr);
  PetscScalar eNodeAddingCounts[4] = { 1, 1, 1, 1 }; // elemental node adding counts
#elif DIM == 3
  ierr = DMDAGetElements_3D (opt->da_nodes, &nel, &nen, &necon);
  CHKERRQ(ierr);
  PetscScalar eNodeAddingCounts[8] = { 1, 1, 1, 1, 1, 1, 1, 1 }; // elemental node adding counts
#endif

  VecZeroEntries (opt->nodeAddingCounts); // zero off nodeAddingCounts vector
  // DMDAGetElements(da_nodes,&nel,&nen,&necon); // Still issue with elemtype
  // change !

  // Get pointer to the densities
  PetscScalar *xp, *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p;
  VecGetArray (opt->xPhys, &xp);
  VecGetArray (opt->xPassive0, &xPassive0p);
  VecGetArray (opt->xPassive1, &xPassive1p);
  VecGetArray (opt->xPassive2, &xPassive2p);
  VecGetArray (opt->xPassive3, &xPassive3p);

  // Edof array, new
  PetscInt edof[nen];

  // Loop over elements
  for (PetscInt i = 0; i < nel; i++) {
    // loop over element nodes
    if (xPassive2p[i] != 0) {
//      memset (eNodeDensity, 0.0, sizeof(eNodeDensity[0]) * nen);
      for (PetscInt j = 0; j < nen; j++) {
        // global numbering of each node
        if (PHYSICS == 0 || PHYSICS == 1) {
          edof[j] = DIM * necon[i * nen + j];
        } else if (PHYSICS == 2) {
          edof[j] = 1 * necon[i * nen + j]; // for heat conduction, only has 1 dof per node
        }
      }
      ierr = VecSetValuesLocal (opt->nodeAddingCounts, nen, edof,
          eNodeAddingCounts, ADD_VALUES);
      CHKERRQ(ierr);
    }
  }

  VecAssemblyBegin (opt->nodeAddingCounts);
  VecAssemblyEnd (opt->nodeAddingCounts);
  PetscScalar tmp;  // need to use this tmp variable because VecSum does not accept PetscInt
  VecSum(opt->nodeAddingCounts, &tmp);
  opt->numNodeLoadAddingCounts = static_cast<PetscInt> (tmp);
  VecAssemblyBegin (opt->nodeAddingCounts);
  VecAssemblyEnd (opt->nodeAddingCounts);
  DMDARestoreElements (opt->da_nodes, &nel, &nen, &necon);

  return ierr;
}

PetscErrorCode PrePostProcess::UpdateNodeDensity (TopOpt *opt) {
  PetscErrorCode ierr = 0;

  // Get the FE mesh structure (from the nodal mesh)
  PetscInt nel, nen;
  const PetscInt *necon;
#if DIM == 2
  ierr = DMDAGetElements_2D (opt->da_nodes, &nel, &nen, &necon);
  CHKERRQ(ierr);
  PetscScalar eNodeAddingCounts[4] = { 1, 1, 1, 1 }; // elemental node adding counts
  PetscScalar eNodeDensity[4]; // elemental node density
#elif DIM == 3
  ierr = DMDAGetElements_3D (opt->da_nodes, &nel, &nen, &necon);
  CHKERRQ(ierr);
  PetscScalar eNodeAddingCounts[8] = { 1, 1, 1, 1, 1, 1, 1, 1 }; // elemental node adding counts
  PetscScalar eNodeDensity[8]; // elemental node density
#endif

  VecZeroEntries (opt->nodeDensity); // zero off nodeDensity vector
  VecZeroEntries (opt->nodeAddingCounts); // zero off nodeAddingCounts vector
  // DMDAGetElements(da_nodes,&nel,&nen,&necon); // Still issue with elemtype
  // change !

  // Get pointer to the densities
  PetscScalar *xp, *xPassive0p, *xPassive1p, *xPassive2p, *xPassive3p;
  VecGetArray (opt->xPhys, &xp);
  VecGetArray (opt->xPassive0, &xPassive0p);
  VecGetArray (opt->xPassive1, &xPassive1p);
  VecGetArray (opt->xPassive2, &xPassive2p);
  VecGetArray (opt->xPassive3, &xPassive3p);

  // Edof array, new
  PetscInt edof[nen];

  // Loop over elements
  for (PetscInt i = 0; i < nel; i++) {
    // loop over element nodes
    if (xPassive0p[i] == 0 /*&& xPassive1p[i] == 0 && xPassive2p[i] == 0
     && xPassive3p[i] == 0*/) {
      memset (eNodeDensity, 0.0, sizeof(eNodeDensity[0]) * nen);
      for (PetscInt j = 0; j < nen; j++) {
        eNodeDensity[j] = xp[i];
        // global numbering of each node
        if (PHYSICS == 0 || PHYSICS == 1) {
          edof[j] = DIM * necon[i * nen + j];
        } else if (PHYSICS == 2) {
          edof[j] = 1 * necon[i * nen + j]; // for heat conduction, only has 1 dof per node
        }
      }
      ierr = VecSetValuesLocal (opt->nodeDensity, nen, edof, eNodeDensity,
          ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValuesLocal (opt->nodeAddingCounts, nen, edof,
          eNodeAddingCounts, ADD_VALUES);
      CHKERRQ(ierr);
    }
  }

  //  Calculate the average node density by using pointwise dividing
  VecAssemblyBegin (opt->nodeAddingCounts);
  VecAssemblyEnd (opt->nodeAddingCounts);
  VecAssemblyBegin (opt->nodeDensity);
  VecAssemblyEnd (opt->nodeDensity);
  VecPointwiseDivide (opt->nodeDensity, opt->nodeDensity,
      opt->nodeAddingCounts);
  VecAssemblyBegin (opt->nodeAddingCounts);
  VecAssemblyEnd (opt->nodeAddingCounts);
  DMDARestoreElements (opt->da_nodes, &nel, &nen, &necon);
  VecAssemblyBegin (opt->nodeDensity);
  VecAssemblyEnd (opt->nodeDensity);

  return ierr;
}

PetscErrorCode PrePostProcess::ImportAndVoxelizeGeometry (TopOpt *opt) {
  PetscErrorCode ierr = 0;

  double t1, t2;

  // Create voxelizer class
  StlVoxelizer *sv = new StlVoxelizer ();
  PetscPrintf (PETSC_COMM_WORLD, "# Start to voxelize the geometries\n");

  // Read and voxelize the geometries
  t1 = MPI_Wtime ();
  if (!opt->inputSTL_DES[0].empty ())
    ierr = sv->Read_file (opt->inputSTL_DES[0]);
  for (unsigned int loadCondition = 0; loadCondition < numLoads;
      ++loadCondition) {
    for (unsigned int backSearch = 0; backSearch <= loadCondition; ++backSearch) {
      if (!opt->inputSTL_FIX[loadCondition - backSearch].empty ()) {
        ierr = sv->Read_file (opt->inputSTL_FIX[loadCondition - backSearch]);
        break;
      }
    }
    for (unsigned int backSearch = 0; backSearch <= loadCondition; ++backSearch) {
      if (!opt->inputSTL_LOD[loadCondition - backSearch].empty ()) {
        ierr = sv->Read_file (opt->inputSTL_LOD[loadCondition - backSearch]);
        break;
      }
    }
  }
  if (!opt->inputSTL_SLD[0].empty ())
    ierr = sv->Read_file (opt->inputSTL_SLD[0]);
  t2 = MPI_Wtime ();
  PetscPrintf (PETSC_COMM_WORLD, "# Read STL files took: %f s\n", t2 - t1);

  t1 = MPI_Wtime ();
  sv->ScaleAndTranslate (nx, ny, nz, dx, dy, dz);
  t2 = MPI_Wtime ();
  PetscPrintf (PETSC_COMM_WORLD, "# Scale and translate took: %f s\n",
      t2 - t1);

  // Part domain voxelization
  if (!opt->inputSTL_DES[0].empty ()) {
    t1 = MPI_Wtime ();
    sv->Voxelize_surface (occDES, nx, ny, nz, dx, dy, dz);
    t2 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD,
        "# Voxelization of DES surface took: %f s\n",
        t2 - t1);
    t1 = MPI_Wtime ();
    sv->Voxelize_solid (occDES, nx, ny, nz);
    t2 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD, "# Voxelization of DES solid took: %f s\n",
        t2 - t1);
    PetscPrintf (PETSC_COMM_WORLD, "# Vexelized %s \n",
        opt->inputSTL_DES[0].c_str ());
  }

  for (unsigned int loadCondition = 0; loadCondition < numLoads;
      ++loadCondition) {
    // Fixture domain voxelization
    for (unsigned int backSearch = 0; backSearch <= loadCondition; ++backSearch) {
      if (!opt->inputSTL_FIX[loadCondition - backSearch].empty ()) {
        t1 = MPI_Wtime ();
        sv->Voxelize_surface (occFIX[loadCondition], nx, ny, nz, dx, dy, dz);
        t2 = MPI_Wtime ();
        PetscPrintf (PETSC_COMM_WORLD,
            "# Voxelization of FIX%d surface took: %f s\n", loadCondition,
            t2 - t1);
        t1 = MPI_Wtime ();
        sv->Voxelize_solid (occFIX[loadCondition], nx, ny, nz);
        t2 = MPI_Wtime ();
        PetscPrintf (PETSC_COMM_WORLD,
            "# Voxelization of FIX%d solid took: %f s\n", loadCondition,
            t2 - t1);
        PetscPrintf (PETSC_COMM_WORLD, "# Vexelized %s \n",
            opt->inputSTL_FIX[loadCondition - backSearch].c_str ());
        break;
      }
    }
    // Load domain voxelization
    for (unsigned int backSearch = 0; backSearch <= loadCondition; ++backSearch) {
      if (!opt->inputSTL_LOD[loadCondition - backSearch].empty ()) {
        t1 = MPI_Wtime ();
        sv->Voxelize_surface (occLOD[loadCondition], nx, ny, nz, dx, dy, dz);
        t2 = MPI_Wtime ();
        PetscPrintf (PETSC_COMM_WORLD,
            "# Voxelization of LOD%d surface took: %f s\n", loadCondition,
            t2 - t1);
        t1 = MPI_Wtime ();
        sv->Voxelize_solid (occLOD[loadCondition], nx, ny, nz);
        t2 = MPI_Wtime ();
        PetscPrintf (PETSC_COMM_WORLD,
            "# Voxelization of LOD%d solid took: %f s\n", loadCondition,
            t2 - t1);
        PetscPrintf (PETSC_COMM_WORLD, "# Vexelized %s \n",
            opt->inputSTL_LOD[loadCondition - backSearch].c_str ());
        break;
      }
    }
  }

  // Non-designable solid domain voxelization
  if (!opt->inputSTL_SLD[0].empty ()) {
    t1 = MPI_Wtime ();
    sv->Voxelize_surface (occSLD, nx, ny, nz, dx, dy, dz);
    t2 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD,
        "# Voxelization of SLD surface took: %f s\n",
        t2 - t1);
    t1 = MPI_Wtime ();
    sv->Voxelize_solid (occSLD, nx, ny, nz);
    t2 = MPI_Wtime ();
    PetscPrintf (PETSC_COMM_WORLD, "# Voxelization of SLD solid took: %f s\n",
        t2 - t1);
    PetscPrintf (PETSC_COMM_WORLD, "# Vexelized %s \n",
        opt->inputSTL_SLD[0].c_str ());
  }

  // Clean up the voxelizer class
  sv->CleanUp ();
  delete sv;

  return ierr;
}

// Passive element assignment
PetscErrorCode
PrePostProcess::AssignPassiveElement (TopOpt *opt)
    {
  PetscErrorCode ierr = 0;

  // Global coordinates of nodes and its pointer
  Vec lcoor;
  const PetscScalar *lcoorp;

  // Global coordinates of elements and its pointer
  Vec elcoor;
  const PetscScalar *elcoorp;

  // Get local coordinates in local node numbering including ghosts
  ierr = DMGetCoordinatesLocal (opt->da_nodes, &lcoor);
  CHKERRQ(ierr);
  VecGetArrayRead (lcoor, &lcoorp);

  // Get local coordinates in local element numbering including ghosts
  ierr = DMGetCoordinatesLocal (opt->da_elem, &elcoor);
  CHKERRQ(ierr);
  VecGetArrayRead (elcoor, &elcoorp);

#if DIM == 2
  // Get pointer to the densities and xpassive index vec
  PetscScalar **xp_2D, **xPassive0p_2D, **xPassive1p_2D, **xPassive2p_2D,
      **xPassive3p_2D;

  // Local vector of the global x and xPassive vectors
  Vec xloc, xPassive0loc, xPassive1loc, xPassive2loc, xPassive3loc;
  DMCreateLocalVector (opt->da_elem, &xloc);
  DMGlobalToLocalBegin (opt->da_elem, opt->x, INSERT_VALUES, xloc);
  DMGlobalToLocalEnd (opt->da_elem, opt->x, INSERT_VALUES, xloc);
  DMCreateLocalVector (opt->da_elem, &xPassive0loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive0, INSERT_VALUES,
      xPassive0loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive0, INSERT_VALUES,
      xPassive0loc);
  DMCreateLocalVector (opt->da_elem, &xPassive1loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive1, INSERT_VALUES,
      xPassive1loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive1, INSERT_VALUES,
      xPassive1loc);
  DMCreateLocalVector (opt->da_elem, &xPassive2loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive2, INSERT_VALUES,
      xPassive2loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive2, INSERT_VALUES,
      xPassive2loc);
  DMCreateLocalVector (opt->da_elem, &xPassive3loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive3, INSERT_VALUES,
      xPassive3loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive3, INSERT_VALUES,
      xPassive3loc);

  DMDAVecGetArray (opt->da_elem, xloc, &xp_2D);
  DMDAVecGetArray (opt->da_elem, xPassive0loc, &xPassive0p_2D);
  DMDAVecGetArray (opt->da_elem, xPassive1loc, &xPassive1p_2D);
  DMDAVecGetArray (opt->da_elem, xPassive2loc, &xPassive2p_2D);
  DMDAVecGetArray (opt->da_elem, xPassive3loc, &xPassive3p_2D);

  PetscInt xs, xe, Xs, Xe;
  PetscInt ys, ye, Ys, Ye;
  //    PetscInt start_x, start_y, start_z, m, n, p;
  DMDAGetCorners (opt->da_elem, &xs, &ys, NULL, &xe, &ye, NULL);
  DMDAGetGhostCorners (opt->da_elem, &Xs, &Ys, NULL, &Xe, &Ye, NULL);

  xe += xs;
  Xe += Xs;
  ye += ys;
  Ye += Ys;

  int occTmp; // occupancy temporary variable
  for (PetscInt i = Xs; i < Xe; i++) {
    for (PetscInt j = Ys; j < Ye; j++) {
      xp_2D[j][i] = 0.0;

      voxIndex = j * nx + i;
      if (!opt->inputSTL_DES[0].empty ()) {
        occTmp = occDES[voxIndex / BATCH];
        if ((occTmp >> (voxIndex % BATCH)) & 1) {
          xp_2D[j][i] = opt->volfrac;
          xPassive0p_2D[j][i] = 0;
          xPassive1p_2D[j][i] = 0;
          xPassive2p_2D[j][i] = 0;
          xPassive3p_2D[j][i] = 0;
        } else {
          xp_2D[j][i] = 0.0;
          xPassive0p_2D[j][i] = 1;
          xPassive1p_2D[j][i] = 0;
          xPassive2p_2D[j][i] = 0;
          xPassive3p_2D[j][i] = 0;
        }
      }

      for (unsigned int loadCondition = 0; loadCondition < numLoads;
          ++loadCondition) {
        for (unsigned int backSearch = 0; backSearch <= loadCondition; ++backSearch) {
          if (!opt->inputSTL_FIX[loadCondition - backSearch].empty ()) {
            occTmp = occFIX[loadCondition][voxIndex / BATCH];
            if ((occTmp >> (voxIndex % BATCH)) & 1) {
              xp_2D[j][i] = 1.0;
              xPassive0p_2D[j][i] = 0;
              xPassive1p_2D[j][i] += 1.0 * std::pow (2, 1.0 * loadCondition);
              xPassive2p_2D[j][i] = 0;
              xPassive3p_2D[j][i] = 0;
            }
            break;
          }
        }
        for (unsigned int backSearch = 0; backSearch <= loadCondition; ++backSearch) {
          if (!opt->inputSTL_LOD[loadCondition - backSearch].empty ()) {
            occTmp = occLOD[loadCondition][voxIndex / BATCH];
            if ((occTmp >> (voxIndex % BATCH)) & 1) {
              xp_2D[j][i] = 1.0;
              xPassive0p_2D[j][i] = 0;
              xPassive1p_2D[j][i] = 0;
              xPassive2p_2D[j][i] += 1.0 * std::pow (2, 1.0 * loadCondition);
              xPassive3p_2D[j][i] = 0;
            }
            break;
          }
        }
      }

      if (!opt->inputSTL_SLD[0].empty ()) {
        occTmp = occSLD[voxIndex / BATCH];
        if ((occTmp >> (voxIndex % BATCH)) & 1) {
          xp_2D[j][i] = 1.0;
          xPassive0p_2D[j][i] = 0;
          xPassive1p_2D[j][i] = 0;
          xPassive2p_2D[j][i] = 0;
          xPassive3p_2D[j][i] = 1;
        }
      }
    }
  }

// Restore the local x and xPassive vectors to their global vectors
  DMLocalToGlobalBegin (opt->da_elem, xloc, INSERT_VALUES, opt->x);
  DMLocalToGlobalEnd (opt->da_elem, xloc, INSERT_VALUES, opt->x);
  DMLocalToGlobalBegin (opt->da_elem, xPassive0loc, INSERT_VALUES,
      opt->xPassive0);
  DMLocalToGlobalEnd (opt->da_elem, xPassive0loc, INSERT_VALUES,
      opt->xPassive0);
  DMLocalToGlobalBegin (opt->da_elem, xPassive1loc, INSERT_VALUES,
      opt->xPassive1);
  DMLocalToGlobalEnd (opt->da_elem, xPassive1loc, INSERT_VALUES,
      opt->xPassive1);
  DMLocalToGlobalBegin (opt->da_elem, xPassive2loc, INSERT_VALUES,
      opt->xPassive2);
  DMLocalToGlobalEnd (opt->da_elem, xPassive2loc, INSERT_VALUES,
      opt->xPassive2);
  DMLocalToGlobalBegin (opt->da_elem, xPassive3loc, INSERT_VALUES,
      opt->xPassive3);
  DMLocalToGlobalEnd (opt->da_elem, xPassive3loc, INSERT_VALUES,
      opt->xPassive3);

  DMDAVecRestoreArray (opt->da_elem, xloc, &xp_2D);
  DMDAVecRestoreArray (opt->da_elem, xPassive0loc, &xPassive0p_2D);
  DMDAVecRestoreArray (opt->da_elem, xPassive1loc, &xPassive1p_2D);
  DMDAVecRestoreArray (opt->da_elem, xPassive2loc, &xPassive2p_2D);
  DMDAVecRestoreArray (opt->da_elem, xPassive3loc, &xPassive3p_2D);

#elif DIM ==3
  // Get pointer to the densities and xpassive index vec
  PetscScalar ***xp_3D, ***xPassive0p_3D, ***xPassive1p_3D, ***xPassive2p_3D,
      ***xPassive3p_3D;

  // Local vector of the global x and xPassive vectors
  Vec xloc, xPassive0loc, xPassive1loc, xPassive2loc, xPassive3loc;
  DMCreateLocalVector (opt->da_elem, &xloc);
  DMGlobalToLocalBegin (opt->da_elem, opt->x, INSERT_VALUES, xloc);
  DMGlobalToLocalEnd (opt->da_elem, opt->x, INSERT_VALUES, xloc);
  DMCreateLocalVector (opt->da_elem, &xPassive0loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive0, INSERT_VALUES,
      xPassive0loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive0, INSERT_VALUES,
      xPassive0loc);
  DMCreateLocalVector (opt->da_elem, &xPassive1loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive1, INSERT_VALUES,
      xPassive1loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive1, INSERT_VALUES,
      xPassive1loc);
  DMCreateLocalVector (opt->da_elem, &xPassive2loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive2, INSERT_VALUES,
      xPassive2loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive2, INSERT_VALUES,
      xPassive2loc);
  DMCreateLocalVector (opt->da_elem, &xPassive3loc);
  DMGlobalToLocalBegin (opt->da_elem, opt->xPassive3, INSERT_VALUES,
      xPassive3loc);
  DMGlobalToLocalEnd (opt->da_elem, opt->xPassive3, INSERT_VALUES,
      xPassive3loc);

  DMDAVecGetArray (opt->da_elem, xloc, &xp_3D);
  DMDAVecGetArray (opt->da_elem, xPassive0loc, &xPassive0p_3D);
  DMDAVecGetArray (opt->da_elem, xPassive1loc, &xPassive1p_3D);
  DMDAVecGetArray (opt->da_elem, xPassive2loc, &xPassive2p_3D);
  DMDAVecGetArray (opt->da_elem, xPassive3loc, &xPassive3p_3D);

  PetscInt xs, xe, Xs, Xe;
  PetscInt ys, ye, Ys, Ye;
  PetscInt zs, ze, Zs, Ze;
  // PetscInt start_x, start_y, start_z, m, n, p;
  DMDAGetCorners (opt->da_elem, &xs, &ys, &zs, &xe, &ye, &ze);
  DMDAGetGhostCorners (opt->da_elem, &Xs, &Ys, &Zs, &Xe, &Ye, &Ze);

  xe += xs;
  Xe += Xs;
  ye += ys;
  Ye += Ys;
  ze += zs;
  Ze += Zs;

  int occTmp; // occupancy temporary variable
  for (PetscInt k = Zs; k < Ze; k++) {
    for (PetscInt j = Ys; j < Ye; j++) {
      for (PetscInt i = Xs; i < Xe; i++) {

        xp_3D[k][j][i] = 0.0;
        voxIndex = k * nx * ny + j * nx + i;

        if (!opt->inputSTL_DES[0].empty ()) {
          occTmp = occDES[voxIndex / BATCH];
          if ((occTmp >> (voxIndex % BATCH)) & 1) {
            xp_3D[k][j][i] = opt->volfrac;
            xPassive0p_3D[k][j][i] = 0;
            xPassive1p_3D[k][j][i] = 0;
            xPassive2p_3D[k][j][i] = 0;
            xPassive3p_3D[k][j][i] = 0;
          } else {
            xp_3D[k][j][i] = 0.0;
            xPassive0p_3D[k][j][i] = 1;
            xPassive1p_3D[k][j][i] = 0;
            xPassive2p_3D[k][j][i] = 0;
            xPassive3p_3D[k][j][i] = 0;
          }
        }

        for (unsigned int loadCondition = 0; loadCondition < numLoads;
            ++loadCondition) {
          for (unsigned int backSearch = 0; backSearch <= loadCondition;
              ++backSearch) {
            if (!opt->inputSTL_FIX[loadCondition - backSearch].empty ()) {
              occTmp = occFIX[loadCondition][voxIndex / BATCH];
              if ((occTmp >> (voxIndex % BATCH)) & 1) {
                xp_3D[k][j][i] = 1.0;
                xPassive0p_3D[k][j][i] = 0;
                xPassive1p_3D[k][j][i] += 1.0
                                          * std::pow (2, 1.0 * loadCondition);
                xPassive2p_3D[k][j][i] = 0;
                xPassive3p_3D[k][j][i] = 0;
              }
              break;
            }
          }
          for (unsigned int backSearch = 0; backSearch <= loadCondition;
              ++backSearch) {
            if (!opt->inputSTL_LOD[loadCondition - backSearch].empty ()) {
              occTmp = occLOD[loadCondition][voxIndex / BATCH];
              if ((occTmp >> (voxIndex % BATCH)) & 1) {
                xp_3D[k][j][i] = 1.0;
                xPassive0p_3D[k][j][i] = 0;
                xPassive1p_3D[k][j][i] = 0;
                xPassive2p_3D[k][j][i] += 1.0
                                          * std::pow (2, 1.0 * loadCondition);
                xPassive3p_3D[k][j][i] = 0;
              }
              break;
            }
          }
        }

        if (!opt->inputSTL_SLD[0].empty ()) {
          occTmp = occSLD[voxIndex / BATCH];
          if ((occTmp >> (voxIndex % BATCH)) & 1) {
            xp_3D[k][j][i] = 1.0;
            xPassive0p_3D[k][j][i] = 0;
            xPassive1p_3D[k][j][i] = 0;
            xPassive2p_3D[k][j][i] = 0;
            xPassive3p_3D[k][j][i] = 1;
          }
        }
      }
    }
  }

// Restore the local x and xPassive vectors to their global vectors
  DMLocalToGlobalBegin (opt->da_elem, xloc, INSERT_VALUES, opt->x);
  DMLocalToGlobalEnd (opt->da_elem, xloc, INSERT_VALUES, opt->x);
  DMLocalToGlobalBegin (opt->da_elem, xPassive0loc, INSERT_VALUES,
      opt->xPassive0);
  DMLocalToGlobalEnd (opt->da_elem, xPassive0loc, INSERT_VALUES,
      opt->xPassive0);
  DMLocalToGlobalBegin (opt->da_elem, xPassive1loc, INSERT_VALUES,
      opt->xPassive1);
  DMLocalToGlobalEnd (opt->da_elem, xPassive1loc, INSERT_VALUES,
      opt->xPassive1);
  DMLocalToGlobalBegin (opt->da_elem, xPassive2loc, INSERT_VALUES,
      opt->xPassive2);
  DMLocalToGlobalEnd (opt->da_elem, xPassive2loc, INSERT_VALUES,
      opt->xPassive2);
  DMLocalToGlobalBegin (opt->da_elem, xPassive3loc, INSERT_VALUES,
      opt->xPassive3);
  DMLocalToGlobalEnd (opt->da_elem, xPassive3loc, INSERT_VALUES,
      opt->xPassive3);

  DMDAVecRestoreArray (opt->da_elem, xloc, &xp_3D);
  DMDAVecRestoreArray (opt->da_elem, xPassive0loc, &xPassive0p_3D);
  DMDAVecRestoreArray (opt->da_elem, xPassive1loc, &xPassive1p_3D);
  DMDAVecRestoreArray (opt->da_elem, xPassive2loc, &xPassive2p_3D);
  DMDAVecRestoreArray (opt->da_elem, xPassive3loc, &xPassive3p_3D);

#endif

  VecRestoreArrayRead (lcoor, &lcoorp);
  VecRestoreArrayRead (elcoor, &elcoorp);
  VecAssemblyBegin (opt->x);
  VecAssemblyEnd (opt->x);
  VecAssemblyBegin (opt->xPassive0);
  VecAssemblyEnd (opt->xPassive0);
  VecAssemblyBegin (opt->xPassive1);
  VecAssemblyEnd (opt->xPassive1);
  VecAssemblyBegin (opt->xPassive2);
  VecAssemblyEnd (opt->xPassive2);
  VecAssemblyBegin (opt->xPassive3);
  VecAssemblyEnd (opt->xPassive3);
  VecDestroy (&xloc);
  VecDestroy (&xPassive0loc);
  VecDestroy (&xPassive1loc);
  VecDestroy (&xPassive2loc);
  VecDestroy (&xPassive3loc);

  return ierr;
}

PetscErrorCode
PrePostProcess::CleanUp ()
{
  PetscErrorCode ierr = 0;

  std::vector<int> ().swap (occDES);
  for (unsigned int loadCondition = 0; loadCondition < numLoads;
      ++loadCondition) {
    std::vector<int> ().swap (occFIX[loadCondition]);
    std::vector<int> ().swap (occLOD[loadCondition]);
  }
  std::vector<int> ().swap (occSLD);

  ierr = occDES.size ();
  for (unsigned int loadCondition = 0; loadCondition < numLoads;
      ++loadCondition) {
    ierr += occFIX[loadCondition].size () + occLOD[loadCondition].size ();
  }
  ierr += occSLD.size ();

  return ierr;
}

#if DIM == 2
PetscErrorCode
PrePostProcess::DMDAGetElements_2D (DM dm, PetscInt *nel,
    PetscInt *nen, const PetscInt *e[]) { // new
  PetscErrorCode ierr = 0;
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

#elif DIM == 3
PetscErrorCode PrePostProcess::DMDAGetElements_3D (DM dm, PetscInt *nel,
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
#endif
