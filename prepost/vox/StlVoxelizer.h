//-------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the TopOptVox authors
//
// This file is part of the TopOptVox.
//
// The TopOptVox is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of TopOptVox.
//
// Author: Zhidong Brian Zhang
// Created on: April 6, 2020
//
// ---------------------------------------------------------------------

#ifndef STLVOXELIZER_H_
#define STLVOXELIZER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <sys/stat.h> // used in inquiry on file size
#include <cmath>

/*
 * BATCH size:
 * sizeof(int) * 8 digits in each occupancy vector item represents different
 * a specific number of voxels
 * For sizeof(int) == 4, each item in the occupancy vector represents
 * 4*8 = 32 voxels. This can help to save memory.
 *
 * Usage e.g.
 * occ[i / BATCH] |= (1 << (i % BATCH));
 * |= is or bitwise operation, 1 << (i % BATCH) is to move 1 to left with
 * digits of (i % BATCH). Thus, one value in occ vector can hold infomation of
 * 32 voxels
 */
#define BATCH 32

// Get the bitwise information of occupancy
#define GET_BIT_OCC(x, i) ((x>>(i%BATCH))&1)

// Throws an error with the given message
#define ERROR_THROW(msg) {std::ostringstream ss; ss << msg; throw(std::runtime_error(ss.str()));}

// Throws an error with the given message conditionally
#define ERROR_COND_THROW(cond, msg) if(cond) { std::ostringstream ss; ss << msg; throw(std::runtime_error(ss.str()));}

// Type names of stl file format
typedef enum {
  BinaryStl, AsciiStl
} StlType;

/**
 * 3D float vector
 */
typedef struct Vector3f {
    Vector3f () {
    }

    Vector3f (float a, float b, float c) {
      this->value[0] = a;
      this->value[1] = b;
      this->value[2] = c;
    }

    Vector3f (const Vector3f &vec) {
      this->value[0] = vec.value[0];
      this->value[1] = vec.value[1];
      this->value[2] = vec.value[2];
    }

    Vector3f& operator= (const Vector3f &vec) {
      this->value[0] = vec.value[0];
      this->value[1] = vec.value[1];
      this->value[2] = vec.value[2];
      return *this;
    }

    float value[3];
} Vector3f;

/**
 * 3D unsigned int vector
 */
typedef struct Vector3ui {
    Vector3ui () {
    }

    Vector3ui (unsigned int a, unsigned int b, unsigned int c) {
      this->value[0] = a;
      this->value[1] = b;
      this->value[2] = c;
    }

    Vector3ui (const Vector3ui &vec) {
      this->value[0] = vec.value[0];
      this->value[1] = vec.value[1];
      this->value[2] = vec.value[2];
    }

    Vector3ui& operator= (const Vector3ui &vec) {
      this->value[0] = vec.value[0];
      this->value[1] = vec.value[1];
      this->value[2] = vec.value[2];
      return *this;
    }

    unsigned int value[3];
} Vector3ui;

/**
 * class StlVoxelizer for voxelizing the mesh imported from stl file
 * Easy-access data structures
 */
class StlVoxelizer {
  public:

    /**
     * initializes the mesh from the stl file
     */
    StlVoxelizer ();
    StlVoxelizer (const char *filename);
    StlVoxelizer (const std::string &filename);

    // fills the mesh with the contents of the specified stl-file
    bool Read_file (const char *filename);
    bool Read_file (const std::string &filename);

    /**
     * Voxelize the given mesh into an occupancy grid.
     * \param[in] occupancy tensor
     * \param[in] mesh number nx, ny, nz
     * \param[in] mesh size dx, dy, dz
     * \param[out] surface voxelization - occupancy tensor
     * \return
     */
    void Voxelize_surface (std::vector<int> &occ, unsigned int nx, unsigned int ny, unsigned int nz, float dx, float dy, float dz);

    /**
     * Voxelize the solid domain based on the surface voxelization.
     * \param[in] occupancy tensor
     * \param[in] mesh number nx, ny, nz
     * \param[out] solid voxelization - occupancy tensor
     * \return
     */
    void Voxelize_solid (std::vector<int> &occ, unsigned int nx, unsigned int ny, unsigned int nz);

    /*
     * Bound adjust, scale and translate the vertices data
     * \param[out] background mesh info: element numbers and sizes
     * \return
     */
    void ScaleAndTranslate (unsigned int nx, unsigned int ny, unsigned int nz, float dx, float dy, float dz);

    /*
     * Clean the occupancy data and free memory
     */
    void CleanUp ();

  private:

    /**
     * Vertices vector with format as (x,y,z).
     */
    std::vector<Vector3f> vertices;

    /**
     * Faces Normals vector
     */
    std::vector<Vector3f> normals;

    /**
     * Triangles data
     */
    std::vector<Vector3ui> tris;

    /**
     * Solids data
     */
    std::vector<unsigned int> solidsRanges;
    unsigned int solidsNumber; // total solids number
    unsigned int solidsItr; // solids iterating number

    /*
     * Domain transform parameters
     */
    float bound[6]; // Input STL geometries bound
    float factor[3]; // Input geometries scaling
    float trans[3]; // Input geometries translation

    /*
     * Occupancy vector size
     */
    unsigned int occSize;

    /*
     * Voxel index, the voxel one-dimensional index
     */
    unsigned int voxIndex;

    /*
     * Occupancy temp vectors, surface, buffer
     */
    std::vector<int> occSUF, occBUF;

    /*
     * Array of elemental neighbour relationship in 2D/3D mesh
     */
    int E_2D[4][2] = { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } };
    int E_3D[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };

    /**
     * Check the type of stl file format
     * Return BinaryStl if it is in binary format
     * Return AsciiStl if it is in ascii format
     * \param[in] filepath path to the STL file
     * \param[out] StlType
     * \return StlType
     */
    StlType GetStlFileFormat (const char *filename);

    /**
     * Get stl file size
     * \param[in] file path to the STL file
     * \param[out] file size
     * \return file size
     */
    size_t GetStlFileSize (const char *filename);

    /**
     * Reading an binary stl file and returning the vertices x, y, z coordinates and
     * the face indices.
     * \param[in] filepath path to the STL file
     * \param[out] vertices vector
     * \param[out] normals vector
     * \param[out] tris vector (triangles)
     * \param[out] solid vector
     * \return success
     */
    bool ReadStlFile_BINARY (const char *filename, std::vector<Vector3f> &verticesOut, std::vector<
                                 Vector3f> &normalsOut, std::vector<Vector3ui> &trisOut, std::vector<
                                 unsigned int> &solidRangesOut);

    /**
     * Reading an ascii stl file and returning the vertices x, y, z coordinates and
     * the face indices.
     * \param[in] filepath path to the STL file
     * \param[out] vertices vector
     * \param[out] normals vector
     * \param[out] tris vector (triangles)
     * \param[out] solid vector
     * \return success
     */
    bool ReadStlFile_ASCII (const char *filename, std::vector<Vector3f> &verticesOut, std::vector<
                                Vector3f> &normalsOut, std::vector<Vector3ui> &trisOut, std::vector<
                                unsigned int> &solidRangesOut);

    /**
     * Fill the buffer domain of the imported CAD
     * \param[in] occupancy tensor
     * \param[in] buffer tensor
     * \param[in] mesh number nx, ny, nz
     * \param[out] buffer tensor
     * \return
     */
    void Fill_buffer (std::vector<int> &occ, unsigned int nx, unsigned int ny, unsigned int nz);

    /**
     * Voxelize the given mesh into an occupancy grid.
     * \param[in] occupancy tensor
     * \param[in] buffer tensor
     * \param[in] mesh number nx, ny, nz
     * \param[out] occupancy tensor
     * \return
     */
    void Fill_solid (std::vector<int> &occ, unsigned int nx, unsigned int ny, unsigned int nz);

    /**
     * Breath First Search (BFS) flood filling.
     * \param[in] occupancy tensor
     * \param[in] buffer tensor
     * \param[in] mesh number nx, ny, nz
     * \param[in] current voxel location, x, y, z
     * \param[out] occupancy tensor
     * \return
     */
    void BFS_flood_fill (std::vector<int> &occ, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int x0, unsigned y0, unsigned z0);

    /**
     * Translate the mesh.
     * \param[in] translation translation vector
     * \return
     */
    void Translate (const float *trans) {
      for (unsigned int v = 0; v < this->num_vertices (); ++v) {
        for (int i = 0; i < 3; ++i) {
          this->vertices[v].value[i] += trans[i];
        }
      }
    }

    /**
     * Scale the mesh.
     * \param[in] scale scale vector
     * \return
     */
    void Scale (const float *factor) {
      for (unsigned int v = 0; v < this->num_vertices (); ++v) {
        for (unsigned int i = 0; i < 3; ++i) {
          this->vertices[v].value[i] *= factor[i];
        }
      }
    }

    /**
     * Get bound of the mesh.
     * \param[in] translation translation vector
     * bd(0,2,4)-min corner, bd(1,3,5)-max corner
     * \return
     */
    void GetBound (float *bound) {
      // Initialization
      for (int i = 0; i < 6; ++i) {
        bound[i] = 0.0;
      }
      // Find bound
      for (unsigned int v = 0; v < this->num_vertices (); ++v) {
        bound[0] = std::min (this->vertices[v].value[0], bound[0]);
        bound[1] = std::max (this->vertices[v].value[0], bound[1]);
        bound[2] = std::min (this->vertices[v].value[1], bound[2]);
        bound[3] = std::max (this->vertices[v].value[1], bound[3]);
        bound[4] = std::min (this->vertices[v].value[2], bound[4]);
        bound[5] = std::max (this->vertices[v].value[2], bound[5]);
      }
    }

    /**
     * Compute triangle box intersection.
     * \param[in] min defining voxel
     * \param[in] max defining voxel
     * \param[in] v1 first vertex
     * \param[in] v2 second vertex
     * \param[in] v3 third vertex
     * \return success
     */
    inline bool Triangle_box_intersection (const Vector3f &min, Vector3f &max, const Vector3f &v1, const Vector3f &v2, const Vector3f &v3);

    /**
     * Get the number of vertices.
     * \return number of vertices
     */
    unsigned int num_vertices () {
      return this->vertices.size ();
    }

    /**
     * Get the number of faces.
     * \return number of faces
     */
    unsigned int num_faces () {
      return this->vertices.size () / 3;
    }
};

#endif // STLVOXELIZER_H_
