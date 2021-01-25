//-------------------------------------------------------------------
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
// Created on: Apr. 2020
// Updated on: Sep. 2020
//
// ---------------------------------------------------------------------

#include "StlVoxelizer.h"

// Point-triangle distance and ray-triangle intersection.
#include "box_triangle/aabb_triangle_overlap.h"
#include "box_triangle/aabb_triangle_overlap_remove_inflation.h"

StlVoxelizer::StlVoxelizer () {
  vertices.clear ();
  normals.clear ();
  tris.clear ();
  solidsRanges.clear ();
  solidsNumber = 0;
  solidsItr = 0;
  occSize = 0;
}

StlVoxelizer::StlVoxelizer (const char *filename) {
  vertices.clear ();
  normals.clear ();
  tris.clear ();
  solidsRanges.clear ();
  solidsNumber = 0;
  solidsItr = 0;
  occSize = 0;

  Read_file (filename);
}

StlVoxelizer::StlVoxelizer (const std::string &filename) {
  vertices.clear ();
  normals.clear ();
  tris.clear ();
  solidsRanges.clear ();
  solidsNumber = 0;
  solidsItr = 0;
  occSize = 0;

  Read_file (filename);
}

bool StlVoxelizer::Read_file (const char *filename) {
  bool res = false;

  StlType stltype = GetStlFileFormat (filename);
  if (stltype == BinaryStl) {
    res = ReadStlFile_BINARY (filename, vertices, normals, tris, solidsRanges);
  } else if (stltype == AsciiStl) {
    res = ReadStlFile_ASCII (filename, vertices, normals, tris, solidsRanges);
  }

  if (!res) {
    vertices.clear ();
    normals.clear ();
    tris.clear ();
    solidsRanges.clear ();
    solidsNumber = 0;
    throw("Error: in ReadStlFile...");
  }

  return res;
}

bool StlVoxelizer::Read_file (const std::string &filename) {
  return Read_file (filename.c_str ());
}

void StlVoxelizer::Voxelize_surface (std::vector<int> &occ, unsigned int nx,
    unsigned int ny, unsigned int nz, float dx, float dy, float dz) {

  unsigned int voxIndex, voxIndex1, voxIndex2;
  unsigned int occTmp, occTmp1, occTmp2;
  bool overlap, overlapInflation;
  float tolerance = 1E-4 * std::min (dx, std::min (dy, dz));
  // Initialize occupancy vector
  occSize = (nx * ny * nz - 1) / BATCH + 1; // occupancy vector size after batched
  occ.clear ();
  occ.resize (occSize);
  occSUF.clear ();
  occSUF.resize (occSize);
  occBUF.clear ();
  occBUF.resize (occSize);

  // triangular space range and vox index range variables
  Vector3f triMax = { 0.0, 0.0, 0.0 };
  Vector3f triMin = { 0.0, 0.0, 0.0 };
  Vector3ui voxMaxLocal = { 0, 0, 0 };
  Vector3ui voxMinLocal = { 0, 0, 0 };
  Vector3f voxSize = { dx, dy, dz };
  Vector3ui voxMaxGlobal = { nx, ny, nz };
  Vector3ui voxMinGlobal = { 0, 0, 0 };

  for (unsigned int l = solidsRanges[2 * solidsItr];
      l < solidsRanges[2 * solidsItr + 1]; ++l) { // loop over facet
    // Get the triangular space range and vox index range, so that the vox
    // range can be reduced. Thus, the voxelization can be accelerated.
    // triangular range
    triMax = vertices[3 * l];
    triMin = vertices[3 * l];
    voxMaxLocal = Vector3ui { 0, 0, 0 };
    voxMinLocal = Vector3ui { 0, 0, 0 };
    for (int dim = 0; dim < 3; ++dim) { // x, y, z
      for (int v = 0; v < 3; ++v) { // vertex 1, 2, 3
        triMax.value[dim] = std::max (triMax.value[dim],
            vertices[3 * l + v].value[dim]);
        triMin.value[dim] = std::min (triMin.value[dim],
            vertices[3 * l + v].value[dim]);
      }
    }

    // vox index range, having buffer of 2 voxels
    for (int dim = 0; dim < 3; ++dim) { // x, y, z
      unsigned int tmp1, tmp2;
      tmp1 = static_cast<unsigned int> (std::ceil (
          triMax.value[dim] / voxSize.value[dim]));
      tmp2 = static_cast<unsigned int> (std::floor (
          triMin.value[dim] / voxSize.value[dim]));
      voxMaxLocal.value[dim] = std::min (tmp1 + 1, voxMaxGlobal.value[dim]);
      voxMinLocal.value[dim] = std::max (tmp2 <= 1 ? 0 : tmp2 - 1,
          voxMinGlobal.value[dim]); // avoid “-1”
    }

    // loop over local voxels
    for (unsigned int k = voxMinLocal.value[2]; k < voxMaxLocal.value[2]; ++k) {
      for (unsigned int j = voxMinLocal.value[1]; j < voxMaxLocal.value[1];
          ++j) {
        for (unsigned int i = voxMinLocal.value[0]; i < voxMaxLocal.value[0];
            ++i) {
          voxIndex = k * ny * nx + j * nx + i;
          // voxel min and max bound
          Vector3f min = { dx * i, dy * j, dz * k };
          Vector3f max = { dx * (i + 1), dy * (j + 1), dz * (k + 1) };
          overlap = Triangle_box_intersection (min, max, vertices[3 * l],
              vertices[3 * l + 1], vertices[3 * l + 2]);
          if (overlap) {
            occSUF[voxIndex / BATCH] |= (1 << (voxIndex % BATCH));
          }
        }
      }
    }
  }

  // Excluding inflation when the triangle normals is parallel to one of the x,y,z axis
  // because when a triangle is coincident with one of the box surfaces
  // this box and its neighbor that sharing the same surface will be counted as solid
  // However, this will cause inflation of the voxelized geometry. So we need to
  // exclude one of these two neighboring boxes based on the normals of the triangle.
  for (unsigned int l = solidsRanges[2 * solidsItr];
      l < solidsRanges[2 * solidsItr + 1]; ++l) { // loop over facet
    // if triangle parallel to x, y, or z
    if (std::abs (normals[l].value[0]) >= 1.0 - tolerance
        || std::abs (normals[l].value[1]) >= 1.0 - tolerance
        || std::abs (normals[l].value[2]) >= 1.0 - tolerance) {
      // Get the triangular space range and vox index range, so that the vox
      // range can be reduced. Thus, the voxelization can be accelerated.
      // triangular range
      std::vector<int> dimseq { 1, 2, 3 }; // dim sequence
      std::vector<unsigned int> nxyzseq { nx, ny, nz }; // voxel maximum range
      std::vector<float> dxyzseq { dx, dy, dz }; // voxel maximum range
      for (int dim = 0; dim < 3; ++dim) { // x, y, z
        triMax = vertices[3 * l];
        triMin = vertices[3 * l];
        voxMaxLocal = Vector3ui { 0, 0, 0 };
        voxMinLocal = Vector3ui { 0, 0, 0 };
        for (int v = 0; v < 3; ++v) { // vertex 1, 2, 3
          triMax.value[dim] = std::max (triMax.value[dim],
              vertices[3 * l + v].value[dim]);
          triMin.value[dim] = std::min (triMin.value[dim],
              vertices[3 * l + v].value[dim]);
        }
        // vox index range
        unsigned int tmp1, tmp2;
        tmp1 = static_cast<unsigned int> (std::ceil (
            triMax.value[dim] / voxSize.value[dim]));
        tmp2 = static_cast<unsigned int> (std::floor (
            triMin.value[dim] / voxSize.value[dim]));
        voxMaxLocal.value[dim] = std::min (tmp1 + 1, voxMaxGlobal.value[dim]);
        voxMinLocal.value[dim] = std::max (tmp2 <= 1 ? 0 : tmp2 - 1,
            voxMinGlobal.value[dim]); // avoid “-1”
        if (std::abs (normals[l].value[dim]) >= 1.0 - tolerance) {
          // swap the index of the normal direction with the end in the sequence
          int tmp = dimseq.back ();
          dimseq.back () = dimseq[dim];
          dimseq[dim] = tmp;
          // swap the voxel maximum range of the normal direction with that of the end
          unsigned int tmp2 = nxyzseq.back ();
          nxyzseq.back () = nxyzseq[dim];
          nxyzseq[dim] = tmp2;
          // swap the voxel size of the normal direction with that of the end
          float tmp3 = dxyzseq.back ();
          dxyzseq.back () = dxyzseq[dim];
          dxyzseq[dim] = tmp3;
        }
      }

      // loop over local voxels
      // normalize the normals vector to 1
      int normal = normals[l].value[dimseq.back ()]
                   / std::abs (normals[l].value[dimseq.back ()]);
      for (unsigned int k = voxMinLocal.value[dimseq[2]];
          k < voxMaxLocal.value[dimseq[2]]; ++k) {
        for (unsigned int j = voxMinLocal.value[dimseq[1]];
            j < voxMaxLocal.value[dimseq[1]]; ++j) {
          for (unsigned int i = voxMinLocal.value[dimseq[0]];
              i < voxMaxLocal.value[dimseq[0]]; ++i) {
            unsigned int front = k, behind = k;
            bool checkInflation = 0;
            if (behind >= 1 && front + 1 <= nxyzseq.back ()) {
              front = front + normal;
              behind = behind - normal;
              voxIndex1 = front * nxyzseq[0] * nxyzseq[1] + j * nxyzseq[0]
                          + i; // voxel next to the current along the normal direction
              occTmp1 = (occSUF[voxIndex1 / BATCH] >> (voxIndex1 % BATCH))
                        & 1;
              voxIndex2 = behind * nxyzseq[0] * nxyzseq[1] + j * nxyzseq[0]
                          + i; // voxel next to the current along the normal direction
              occTmp2 = (occSUF[voxIndex1 / BATCH] >> (voxIndex2 % BATCH))
                        & 1;
              voxIndex = k * nxyzseq[0] * nxyzseq[1] + j * nxyzseq[0] + i;
              occTmp = (occSUF[voxIndex / BATCH] >> (voxIndex % BATCH)) & 1;
              // only when front is void, current and behind are solid
              checkInflation = ((occTmp1==0) & occTmp & occTmp2);
            }
            // voxel min and max bound
            Vector3f min = { dx * i, dy * j, dz * k };
            Vector3f max = { dx * (i + 1), dy * (j + 1), dz * (k + 1) };
            if (checkInflation) {
              overlapInflation = Triangle_box_intersection_remove_inflation (
                  min, max, vertices[3 * l], vertices[3 * l + 1],
                  vertices[3 * l + 2]);
              if (overlapInflation) {
                occSUF[voxIndex / BATCH] &= (~(1 << (voxIndex % BATCH))); // remove the inflated voxel
              }
            }
          }
        }
      }
    }
  }
  solidsItr++;
}

void StlVoxelizer::Voxelize_solid (std::vector<int> &occ, unsigned int nx,
    unsigned int ny, unsigned int nz) {
  Fill_buffer (nx, ny, nz);
  Fill_solid (occ, nx, ny, nz);
}

void StlVoxelizer::ScaleAndTranslate (unsigned int nx, unsigned int ny,
    unsigned int nz, float dx, float dy, float dz) {
// Get bound of the input models
  GetBound (bound);

  if ((bound[5] - bound[4]) == 0) { // z direction dimension 0, meaning it is 2D model, otherwise 3D
    // Scaling
    factor[0] = std::min (nx * dx / (bound[1] - bound[0]),
        ny * dy / (bound[3] - bound[2])); // scale the CAD model of the  part to the prescribe domain size
    factor[1] = factor[0];
    factor[2] = 1.0;
    Scale (factor);
    // Get bound of the input models
    GetBound (bound);
    // Translate the input models
    trans[0] = -bound[0];
    trans[1] = -bound[2];
    trans[2] = -bound[4];
    Translate (trans);
  } else {
    // Scaling
    factor[0] = std::min (nx * dx / (bound[1] - bound[0]),
        std::min (ny * dy / (bound[3] - bound[2]),
            nz * dz / (bound[5] - bound[4]))); // scale the CAD model of the  part to the prescribe domain size
    factor[1] = factor[0];
    factor[2] = factor[0];
    Scale (factor);
    // Get bound of the input models
    GetBound (bound);
    // Translate the input models
    trans[0] = -bound[0];
    trans[1] = -bound[2];
    trans[2] = -bound[4];
    Translate (trans);
  }
}

void StlVoxelizer::CleanUp () {
  std::vector<int> ().swap (occSUF);
  std::vector<int> ().swap (occBUF);
}

//##############################################################################
//#                               Private                                      #
//##############################################################################

StlType
StlVoxelizer::GetStlFileFormat (const char *filename) {
// Load the stl file
  std::ifstream in (filename);
  if (!in) {
    ERROR_THROW("Error during open the stl file...");
  }

  size_t fileSize = GetStlFileSize (filename); // the stl file size
  std::string line; // line string variable
  std::vector<char> charBuf; // char array buf stored in a vector
// Normals 3 float (3*4 bytes) + Vertices 3 x 3float (3*3*4 bytes) + AttributeCount 1 short (2 bytes) = 50 bytes
  const size_t eachTriangleSize = 3 * sizeof(float_t) + 3 * 3 * sizeof(float_t)
                                  + sizeof(uint16_t);

// ASCII
// 1. check file size, "solid \n" + "endsolid" should be at least = 15 bytes
  if (fileSize < 15) {
    std::string errStr;
    std::ostringstream oss;
    oss << "Stl file is not long enough (" << static_cast<unsigned int> (fileSize) << "bytes)... /n";
    errStr = oss.str ();
    in.close (); // close the file
    ERROR_THROW(errStr);
  }
// 2. check first 6 bytes whether is "solid "
  std::getline (in, line);
  if (line.compare (0, 6, "solid ") == 0) {
    // 3. then check whether ther is "endsolid"
    while (getline (in, line)) {
      if (line.compare (0, 8, "endsolid") == 0) {
        in.close ();
        return AsciiStl; // then it is ASCII
      }
    }
  }
  in.close ();

// BINARY
  in.open (filename, std::ios::binary); // reload the file stream as binary
// 1. check file size, 80 bytes (header size) + 4 bytes (number of triangular) should be at least = 84 bytes
  if (fileSize < 84) {
    std::string errStr;
    std::ostringstream oss;
    oss << "Stl file is not long enough (" << static_cast<unsigned int> (fileSize) << "bytes)... /n";
    errStr = oss.str ();
    in.close (); // close the file
    ERROR_THROW(errStr);
  }
// 2. check whether the total file size matches with the computed value
  in.seekg (80); // header is from 0 to 79, numTriangle starts from offset 80 to 83
  unsigned int numTriangles; // number of triangles
  in.read (reinterpret_cast<char*> (&numTriangles), 4);
  if (fileSize == (80 + 4 + eachTriangleSize * numTriangles)) {
    in.close ();
    return BinaryStl;
  }

// If not ASCII nor BINARY, throw an error
  ERROR_THROW("Error during open the stl file for both ASCII and BINARY ...");
}

size_t
StlVoxelizer::GetStlFileSize (const char *filename) {
  struct stat st;
  int statReturn = stat (filename, &st); // if succeed, statReturn = 0
  return (statReturn == 0) ? st.st_size : 0;
}

bool StlVoxelizer::ReadStlFile_BINARY (const char *filename, std::vector<
    Vector3f> &verticesOut, std::vector<
    Vector3f> &normalsOut, std::vector<
    Vector3ui> &trisOut, std::vector<
    unsigned int> &solidRangesOut) {
  std::ifstream in (filename, std::ios::binary); // reload the file stream as binary
  if (!in) {
    in.close ();
    ERROR_THROW("Error during open the stl file...");
  }

  char header[80];
  in.read (header, 80); // read the header
  if (!in)
  ERROR_THROW("Error during reading header in the stl file...");
  unsigned int numTriangles; // the number of triangles
  in.read (reinterpret_cast<char*> (&numTriangles), 4);
  if (!in)
  ERROR_THROW("Error during reading number of triangles in the stl file...");

// read the stl facet by facet
  unsigned int f = 0; // the index of the triangle data
  Vector3f nl, cr;
  Vector3ui tv;
  float buf[12];
  short att;
  solidRangesOut.push_back (trisOut.size ()); // save solid range to be able to deal with multiple solids: start
  while (f < numTriangles) {
    in.read (reinterpret_cast<char*> (buf), 12 * 4);
    if (!in)
    ERROR_THROW("Error during reading info. of a triangle in the stl file...");

    // 3 normals
    for (int i = 0; i < 3; ++i) {
      nl.value[i] = buf[i];
    }
    normalsOut.push_back (nl);

    // 3 vertices
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        cr.value[i] = buf[3 + 3 * j + i];
      }
      verticesOut.push_back (cr);
    }
    // 1 short of attribute
    in.read (reinterpret_cast<char*> (&att), 2);
    if (!in)
    ERROR_THROW("Error during reading attribute with 1 short size...");

    for (size_t i = 0; i < 3; ++i)
      tv.value[i] = (verticesOut.size () - 3) + i;
    trisOut.push_back (tv);
    f++;
  }
  solidRangesOut.push_back (trisOut.size ()); // save solid range to be able to deal with multiple solids: end
  solidsNumber++;

  return true;
}

bool StlVoxelizer::ReadStlFile_ASCII (const char *filename,
    std::vector<Vector3f> &verticesOut, std::vector<
        Vector3f> &normalsOut, std::vector<
        Vector3ui> &trisOut, std::vector<
        unsigned int> &solidRangesOut) {

  std::ifstream in (filename);
  if (!in) {
    in.close ();
    ERROR_THROW("Error during open the stl file...");
  }

  std::string line;
  std::vector<std::string> tokens;
  int lineCount = 1;
  size_t numFaceVertices = 0;

// read the file line-by-line
  while (std::getline (in, line)) {
    std::istringstream iss (line); // input string stream
    int tokenCount = 0; // number of tokens
    tokens.clear ();
    tokens.resize (tokenCount + 1);

    // read the line token-by-token
    while (iss >> tokens[tokenCount]) {
      ++tokenCount;
      tokens.resize (tokenCount + 1);
    }

    if (tokenCount > 0) {
      std::string &token = tokens[0];
      if (token.compare ("vertex") == 0) {
        // read the vertices position
        if (tokenCount < 4)
        ERROR_THROW("Error of vertex number in line " << lineCount);
        Vector3f cr;
        for (size_t i = 0; i < 3; ++i)
          cr.value[i] = static_cast<float> (std::stof (tokens[i + 1]));
        verticesOut.push_back (cr);
        ++numFaceVertices;
      } else if (token.compare ("facet") == 0) {
        // read the normal
        Vector3f nl;
        for (size_t i = 0; i < 3; ++i)
          nl.value[i] = static_cast<float> (std::stof (tokens[i + 2]));
        normalsOut.push_back (nl);
        numFaceVertices = 0;
      } else if (token.compare ("endfacet") == 0) {
        // triangular vertices sequential number
        if (numFaceVertices < 3)
        ERROR_THROW("Error of face vertex number in line " << lineCount);
        Vector3ui tv;
        for (size_t i = 0; i < 3; ++i)
          tv.value[i] = (verticesOut.size () - 3) + i;
        trisOut.push_back (tv);
      } else if (token.compare ("solid") == 0) {
        solidRangesOut.push_back (trisOut.size ());
      } else if (token.compare ("endsolid") == 0) {
        solidRangesOut.push_back (trisOut.size ());
        solidsNumber++;
      }
    }
    lineCount++;
  }

  return true;
}

void StlVoxelizer::Fill_buffer (unsigned int nx, unsigned int ny,
    unsigned int nz) {
  unsigned int x, y, z;
  unsigned int occTmp;
  unsigned int voxIndex;
// fill XY
  for (unsigned int k = 0; k < nz; ++k) {
    for (unsigned int j = 0; j < ny; ++j) {
      for (unsigned int i = 0; i < nx; ++i) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if voxel is surface, then break
          break;
        else {
          // Start flood fill from the voxel
          BFS_flood_fill_buffer (nx, ny, nz, i, j, k);
          // Search its neighbor voxels
          for (unsigned int e = 0; e < 4; ++e) {
            x = i + E_2D[e][0];
            y = j + E_2D[e][1];
            z = k;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_buffer (nx, ny, nz, x, y, z);
            }
          }
        }
      }

      for (unsigned int i = nx - 1; i != static_cast<unsigned int> (-1); --i) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if voxel is sur, then break
          break;
        else {
          // Start flood fill from the voxel
          BFS_flood_fill_buffer (nx, ny, nz, i, j, k);
          // Search its neighbor voxels
          for (unsigned int e = 0; e < 4; ++e) {
            x = i + E_2D[e][0];
            y = j + E_2D[e][1];
            z = k;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_buffer (nx, ny, nz, x, y, z);
            }
          }
        }
      }
    }
  }

// fill YZ
  for (unsigned int i = 0; i < nx; ++i) {
    for (unsigned int k = 0; k < nz; ++k) {
      for (unsigned int j = 0; j < ny; ++j) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if voxel is sur, then break
          break;
        else {
          // Start flood fill from the voxel
          BFS_flood_fill_buffer (nx, ny, nz, i, j, k);
          // Search its neighbor voxels
          for (unsigned int e = 0; e < 4; ++e) {
            y = j + E_2D[e][0];
            z = k + E_2D[e][1];
            x = i;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_buffer (nx, ny, nz, x, y, z);
            }
          }
        }
      }

      for (unsigned int j = ny - 1; j != static_cast<unsigned int> (-1); --j) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if voxel is sur, then break
          break;
        else {
          // Start flood fill from the voxel
          BFS_flood_fill_buffer (nx, ny, nz, i, j, k);
          // Search its neighbor voxels
          for (unsigned int e = 0; e < 4; ++e) {
            y = j + E_2D[e][0];
            z = k + E_2D[e][1];
            x = i;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_buffer (nx, ny, nz, x, y, z);
            }
          }
        }
      }
    }
  }

// fill ZX
  for (unsigned int j = 0; j < ny; ++j) {
    for (unsigned int i = 0; i < nx; ++i) {
      for (unsigned int k = 0; k < nz; ++k) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if voxel is sur, then break
          break;
        else {
          // Start flood fill from the voxel
          BFS_flood_fill_buffer (nx, ny, nz, i, j, k);
          // Search its neighbor voxels
          for (unsigned int e = 0; e < 4; ++e) {
            z = k + E_2D[e][0];
            x = i + E_2D[e][1];
            y = j;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_buffer (nx, ny, nz, x, y, z);
            }
          }
        }
      }

      for (unsigned int k = nz - 1; k != static_cast<unsigned int> (-1); --k) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if voxel is sur, then break
          break;
        else {
          // Start flood fill from the voxel
          BFS_flood_fill_buffer (nx, ny, nz, i, j, k);
          // Search its neighbor voxels
          for (unsigned int e = 0; e < 4; ++e) {
            z = k + E_2D[e][0];
            x = i + E_2D[e][1];
            y = j;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_buffer (nx, ny, nz, x, y, z);
            }
          }
        }
      }
    }
  }
}

void StlVoxelizer::Fill_solid (std::vector<int> &occ, unsigned int nx,
    unsigned int ny, unsigned int nz) {
  unsigned int x, y, z;
  unsigned int occTmp;
  unsigned int voxIndex;
// fill XY
  for (unsigned int k = 0; k < nz; ++k) {
    for (unsigned int j = 0; j < ny; ++j) {
      for (unsigned int i = 0; i < nx; ++i) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occ[voxIndex / BATCH] & occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1) // if occ is true and on surface, then break
          break;
        else {
          occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                    >> (voxIndex % BATCH))
                   & 1;
          if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, i, j, k);
          for (unsigned int e = 0; e < 4; ++e) {
            x = i + E_2D[e][0];
            y = j + E_2D[e][1];
            z = k;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, x, y, z);
            }
          }
        }
      }

      for (unsigned int i = nx - 1; i != static_cast<unsigned int> (-1); --i) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occ[voxIndex / BATCH] & occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1)
          break;
        else {
          occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                    >> (voxIndex % BATCH))
                   & 1;
          if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, i, j, k);
          for (unsigned int e = 0; e < 4; ++e) {
            x = i + E_2D[e][0];
            y = j + E_2D[e][1];
            z = k;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, x, y, z);
            }
          }
        }
      }
    }
  }

// fill YZ
  for (unsigned int i = 0; i < nx; ++i) {
    for (unsigned int k = 0; k < nz; ++k) {
      for (unsigned int j = 0; j < ny; ++j) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occ[voxIndex / BATCH] & occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1)
          break;
        else {
          occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                    >> (voxIndex % BATCH))
                   & 1;
          if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, i, j, k);
          for (unsigned int e = 0; e < 4; ++e) {
            y = j + E_2D[e][0];
            z = k + E_2D[e][1];
            x = i;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, x, y, z);
            }
          }
        }
      }

      for (unsigned int j = ny - 1; j != static_cast<unsigned int> (-1); --j) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occ[voxIndex / BATCH] & occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1)
          break;
        else {
          occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                    >> (voxIndex % BATCH))
                   & 1;
          if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, i, j, k);
          for (unsigned int e = 0; e < 4; ++e) {
            y = j + E_2D[e][0];
            z = k + E_2D[e][1];
            x = i;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, x, y, z);
            }
          }
        }
      }
    }
  }

// fill XZ
  for (unsigned int j = 0; j < ny; ++j) {
    for (unsigned int i = 0; i < nx; ++i) {
      for (unsigned int k = 0; k < nz; ++k) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occ[voxIndex / BATCH] & occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1)
          break;
        else {
          occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                    >> (voxIndex % BATCH))
                   & 1;
          if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, i, j, k);
          for (unsigned int e = 0; e < 4; ++e) {
            y = j + E_2D[e][0];
            z = k + E_2D[e][1];
            x = i;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, x, y, z);
            }
          }
        }
      }

      for (unsigned int k = nz - 1; k != static_cast<unsigned int> (-1); --k) {
        voxIndex = k * ny * nx + j * nx + i;
        occTmp = ((occ[voxIndex / BATCH] & occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 1)
          break;
        else {
          occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                    >> (voxIndex % BATCH))
                   & 1;
          if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, i, j, k);
          for (unsigned int e = 0; e < 4; ++e) {
            y = j + E_2D[e][0];
            z = k + E_2D[e][1];
            x = i;
            if (x >= 0 && x < nx && y >= 0 && y < ny) {
              voxIndex = k * ny * nx + y * nx + x;
              occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                        >> (voxIndex % BATCH))
                       & 1;
              if (occTmp == 0) BFS_flood_fill_solid (occ, nx, ny, nz, x, y, z);
            }
          }
        }
      }
    }
  }
}

void StlVoxelizer::BFS_flood_fill_buffer (unsigned int nx, unsigned int ny,
    unsigned int nz, unsigned int x0, unsigned int y0, unsigned int z0) {
  unsigned int voxIndex;
// Queue for recording breadth first search
  std::queue<Vector3ui> q;

// Initialize the seed voxel and mark it as solid
  Vector3ui voxTemp = { x0, y0, z0 };
  q.push (voxTemp);
  voxIndex = z0 * ny * nx + y0 * nx + x0;
  occBUF[voxIndex / BATCH] |= (1 << (voxIndex % BATCH));

  unsigned int x, y, z;
  unsigned int occTmp;
  while (!q.empty ()) {
    voxTemp = q.front ();
    q.pop ();
    x0 = voxTemp.value[0];
    y0 = voxTemp.value[1];
    z0 = voxTemp.value[2];

    for (unsigned int e = 0; e < 6; ++e) {
      x = x0 + E_3D[e][0];
      y = y0 + E_3D[e][1];
      z = z0 + E_3D[e][2];

      if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz) {
        voxIndex = z * ny * nx + y * nx + x;
        occTmp = ((occBUF[voxIndex / BATCH] | occSUF[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 0) { // not buffer nor occupied, then put it to queue
          voxTemp.value[0] = x;
          voxTemp.value[1] = y;
          voxTemp.value[2] = z;
          q.push (voxTemp);
          occBUF[voxIndex / BATCH] |= (1 << (voxIndex % BATCH)); // mark the voxel as solid
        }
      }
    }
  }
}

void StlVoxelizer::BFS_flood_fill_solid (std::vector<int> &occ, unsigned int nx,
    unsigned int ny, unsigned int nz, unsigned int x0, unsigned int y0,
    unsigned int z0) {
  unsigned int voxIndex;
// Queue for recording breadth first search
  std::queue<Vector3ui> q;

// Initialize the seed voxel and mark it as solid
  Vector3ui voxTemp = { x0, y0, z0 };
  q.push (voxTemp);
  voxIndex = z0 * ny * nx + y0 * nx + x0;
  occ[voxIndex / BATCH] |= (1 << (voxIndex % BATCH));

  unsigned int x, y, z;
  unsigned int occTmp;
  while (!q.empty ()) {
    voxTemp = q.front ();
    q.pop ();
    x0 = voxTemp.value[0];
    y0 = voxTemp.value[1];
    z0 = voxTemp.value[2];

    for (unsigned int e = 0; e < 6; ++e) {
      x = x0 + E_3D[e][0];
      y = y0 + E_3D[e][1];
      z = z0 + E_3D[e][2];

      if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz) {
        voxIndex = z * ny * nx + y * nx + x;
        occTmp = ((occBUF[voxIndex / BATCH] | occ[voxIndex / BATCH])
                  >> (voxIndex % BATCH))
                 & 1;
        if (occTmp == 0) { // not buffer nor occupied, then put it to queue
          voxTemp.value[0] = x;
          voxTemp.value[1] = y;
          voxTemp.value[2] = z;
          q.push (voxTemp);
          occ[voxIndex / BATCH] |= (1 << (voxIndex % BATCH)); // mark the voxel as solid
        }
      }
    }
  }
}

inline bool StlVoxelizer::Triangle_box_intersection (const Vector3f &min,
    Vector3f &max, const Vector3f &v1, const Vector3f &v2, const Vector3f &v3) {
  float half_size[3] = { (max.value[0] - min.value[0]) / (float) 2.0, (max.value[1]
      - min.value[1])
                                                                      / (float) 2.0, (max.value[2]
      - min.value[2])
                                                                                     / (float) 2.0 };
  float center[3] = { max.value[0] - half_size[0], max.value[1] - half_size[1], max.value[2]
      - half_size[2] };
  float vertices[3][3] = { { v1.value[0], v1.value[1], v1.value[2] }, { v2.value[0], v2.value[1], v2.value[2] }, { v3.value[0], v3.value[1], v3.value[2] } };
  return triBoxOverlap (center, half_size, vertices);
}

inline bool StlVoxelizer::Triangle_box_intersection_remove_inflation (
    const Vector3f &min,
    Vector3f &max, const Vector3f &v1, const Vector3f &v2, const Vector3f &v3) {
  bool removeInflationFlag = 0;
  float half_size[3] = { (max.value[0] - min.value[0]) / (float) 2.0, (max.value[1]
      - min.value[1])
                                                                      / (float) 2.0, (max.value[2]
      - min.value[2])
                                                                                     / (float) 2.0 };
  float center[3] = { max.value[0] - half_size[0], max.value[1] - half_size[1], max.value[2]
      - half_size[2] };
  float vertices[3][3] = { { v1.value[0], v1.value[1], v1.value[2] }, { v2.value[0], v2.value[1], v2.value[2] }, { v3.value[0], v3.value[1], v3.value[2] } };
// Judge to remove inflation
  removeInflationFlag = triBoxOverlapRemoveInflation (center, half_size,
      vertices);
  return removeInflationFlag;
}
