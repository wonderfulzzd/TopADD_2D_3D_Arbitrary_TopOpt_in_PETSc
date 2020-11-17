/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/

/*
 * Modified by Zhidong Brian Zhang in Nov 2020, University of Waterloo
 *
 * 1. Adapted to aabb_triangle_overlap_remove_inflation. The inflation of the
 *  voxelized geometry will be corrected by this file
 *
 */

#ifndef _AABB_TRIANGLE_OVERLAP_REMOVE_INFLATION_H_
#define _AABB_TRIANGLE_OVERLAP_REMOVE_INFLATION_H_

#include <math.h>
#include <stdio.h>

inline int planeBoxOverlapRemoveInflation (float normal[3], float vert[3],
    float maxbox[3]) // -NJMP-
    {
  int q;
  float vmin[3], vmax[3], v;
  float tolerance = 1E-4
      * std::min (2 * maxbox[0], std::min (2 * maxbox[1], 2 * maxbox[2])); // # new

  for (q = X; q <= Z; q++)
      {
    v = vert[q]; // -NJMP-
    if (normal[q] > 0.0f)
        {
      vmin[q] = -maxbox[q] - v; // -NJMP-
      vmax[q] = maxbox[q] - v; // -NJMP-
    }
    else
    {
      vmin[q] = maxbox[q] - v; // -NJMP-
      vmax[q] = -maxbox[q] - v; // -NJMP-
    }
  }
  if (std::abs (normal[0]) >= 1.0 - tolerance
      || std::abs (normal[1]) >= 1.0 - tolerance
      || std::abs (normal[2]) >= 1.0 - tolerance) {
    if (DOT (normal, vmax) < -2 * std::abs (DOT (normal, maxbox)) + tolerance
        || DOT (normal, vmax) > 2 * std::abs (DOT (normal, maxbox)) - tolerance)
      return 1;
    if (DOT (normal, vmin) > -tolerance
        || DOT (normal, vmin) < -2 * std::abs (DOT (normal, maxbox))
                                - tolerance)
      return 1;
  }
  return 0;
}

inline int triBoxOverlapRemoveInflation (float boxcenter[3],
    float boxhalfsize[3],
    float triverts[3][3])
    {

  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
  /*       we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
  /*       this gives 3x3=9 more tests */
  float v0[3], v1[3], v2[3];
//   float axis[3];
  float min, max, p0, p1, p2, rad, fex, fey, fez; // -NJMP- "d" local variable removed
  float normal[3], e0[3], e1[3], e2[3];

  float tolerance = 1E-4
      * std::min (2 * boxhalfsize[0],
          std::min (2 * boxhalfsize[1], 2 * boxhalfsize[2])); // # new

      /* This is the fastest branch on Sun */
  /* move everything so that the boxcenter is in (0,0,0) */
  SUB (v0, triverts[0], boxcenter);
  SUB (v1, triverts[1], boxcenter);
  SUB (v2, triverts[2], boxcenter);

  /* compute triangle edges */
  SUB (e0, v1, v0); /* tri edge 0 */
  SUB (e1, v2, v1); /* tri edge 1 */
  SUB (e2, v0, v2); /* tri edge 2 */

  /* Bullet 3:  */
  /*  test the 9 tests first (this was faster) */
  fex = fabsf (e0[X]);
  fey = fabsf (e0[Y]);
  fez = fabsf (e0[Z]);
  AXISTEST_X01 (e0[Z], e0[Y], fez, fey);
  AXISTEST_Y02 (e0[Z], e0[X], fez, fex);
  AXISTEST_Z12 (e0[Y], e0[X], fey, fex);

  fex = fabsf (e1[X]);
  fey = fabsf (e1[Y]);
  fez = fabsf (e1[Z]);
  AXISTEST_X01 (e1[Z], e1[Y], fez, fey);
  AXISTEST_Y02 (e1[Z], e1[X], fez, fex);
  AXISTEST_Z0 (e1[Y], e1[X], fey, fex);

  fex = fabsf (e2[X]);
  fey = fabsf (e2[Y]);
  fez = fabsf (e2[Z]);
  AXISTEST_X2 (e2[Z], e2[Y], fez, fey);
  AXISTEST_Y1 (e2[Z], e2[X], fez, fex);
  AXISTEST_Z12 (e2[Y], e2[X], fey, fex);

  /* Bullet 1: */
  /*  first test overlap in the {x,y,z}-directions */
  /*  find min, max of the triangle each direction, and test for overlap in */
  /*  that direction -- this is equivalent to testing a minimal AABB around */
  /*  the triangle against the AABB */

  /* test in X-direction */
  FINDMINMAX (v0[X], v1[X], v2[X], min, max);
  if (min > boxhalfsize[X] + tolerance || max < -boxhalfsize[X] - tolerance) // # modified
  return 0;

  /* test in Y-direction */
  FINDMINMAX (v0[Y], v1[Y], v2[Y], min, max);
  if (min > boxhalfsize[Y] + tolerance || max < -boxhalfsize[Y] - tolerance) // # modified
  return 0;

  /* test in Z-direction */
  FINDMINMAX (v0[Z], v1[Z], v2[Z], min, max);
  if (min > boxhalfsize[Z] + tolerance || max < -boxhalfsize[Z] - tolerance) // # modified
  return 0;

  /* Bullet 2: */
  /*  test if the box intersects the plane of the triangle */
  /*  compute plane equation of triangle: normal*x+d=0 */
  CROSS (normal, e0, e1);
// Normalize the normal vector
  float normleng = std::sqrt (
      normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]); // # new
  normal[0] = normal[0] / normleng; // # new
  normal[1] = normal[1] / normleng; // # new
  normal[2] = normal[2] / normleng; // # new
// -NJMP- (line removed here)
  if (!planeBoxOverlapRemoveInflation (normal, v0, boxhalfsize)) return 0; // -NJMP-

  return 1; /* box and triangle overlaps - remove inflation*/
}

#endif
