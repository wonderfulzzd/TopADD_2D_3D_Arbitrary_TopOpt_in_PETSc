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

#ifndef PrePostProcess_H_
#define PrePostProcess_H_

// Topology optimization class
#include <TopOpt.h>

#include "options.h"
// Stl voxelizer
#include <./vox/StlVoxelizer.h>

/**
 * class Pre- and post-processing class
 */
class PrePostProcess {
  public:

    /**
     * Constructor
     */
    PrePostProcess (TopOpt *opt);

    /**
     * Destructor
     */
    ~PrePostProcess ();

    /*
     * Voxel occupancy info
     */
    std::vector<int> occDES;
    std::vector<std::vector<int>> occFIX;
    std::vector<std::vector<int>> occLOD;
    std::vector<int> occSLD;

    /*
     * Number of loading conditions
     */
    unsigned int numLoads;

    /*
     * Voxel index, the voxel one-dimensional index
     */
    unsigned int voxIndex;

    /*
     *  Mesh parameters
     */
    unsigned int nx, ny, nz; // Voxel number in x, y, z
    float dx, dy, dz; // Voxel size in x, y, z

    /**
     * Design domain initialization
     * \param[in] pointer of the TopOpt class
     * \param[out]
     * \return PetscErrorCode
     */
    PetscErrorCode DesignDomainInitialization (TopOpt *opt);

    /**
     * Update the node density
     * \param[in] pointer of the TopOpt class
     * \param[out]
     * \return PetscErrorCode
     */
    PetscErrorCode UpdateNodeDensity (TopOpt *opt);

  private:

    /**
     * Import and voxelize the geometry
     * \param[in] vector of the voxelization information, 0 represents void, 1 represents solid
     * \param[out]
     * \return PetscErrorCode
     */
    PetscErrorCode ImportAndVoxelizeGeometry (TopOpt *opt);

    /**
     * Passive element assignment
     * \param[in] pointer of the TopOpt class
     * \param[out]
     * \return PetscErrorCode
     */
    PetscErrorCode AssignPassiveElement (TopOpt *opt);

    /*
     * Clean the occupancy data and free memory
     */
    PetscErrorCode CleanUp ();

    /**
     * Get the DMDA mesh info
     * \param[in] PETSc DM
     * \param[out] *nel, pointer to the number of element
     * \param[out] *nen, pointer to the number of nodes in each element
     * \param[out] *e[], pointer to the global indices of the elements' vertices
     * \return PetscErrorCode
     */
#if DIM == 2
    PetscErrorCode DMDAGetElements_2D (
        DM dm, PetscInt *nel,
        PetscInt *nen,
        const PetscInt *e[]); // new
#elif DIM == 3
    PetscErrorCode DMDAGetElements_3D (DM dm, PetscInt *nel, PetscInt *nen,
        const PetscInt *e[]);
#endif

};

#endif /* PrePostProcess_H_ */
