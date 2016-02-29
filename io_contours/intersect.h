// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_INTERSECT_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_INTERSECT_H__
#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include "contour_types.h"
#include <terrastream/common/tflow_types.h>
#include <terrastream/common/wlabel.h>

namespace terrastream{

  //Takes the input stream of triangles and creates the contour segments resulting from intersecting
  //the triangles with the contour planes with displacement equal to the granularity.
  //Each segment is marked with a sign, which indicates whether
  //contours lying on an edge of a triangle was generated as the bottom edge or top edge of the triangle.
  //Finally a map_info struct is returned, this is for later ridge removal use.
  map_info intersect(stream<triangle> &input,elev_t granularity,stream<signed_contour_segment> &output);

  //Does the same as the other intersect, but for every level t, edges for contours on levels
  //t-z_diff and t+z_diff are added.
  map_info intersect(stream<triangle> &input,elev_t granularity,float z_diff,
			 stream<signed_contour_segment> &output);
}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_INTERSECT_H__*/
