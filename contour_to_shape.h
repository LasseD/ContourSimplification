// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_CONTOUR_TO_SHAPE_H__
#define __TEST_CONTOUR_SIMPLIFICATION_CONTOUR_TO_SHAPE_H__
#include "io_contours/contour_types.h"
#include <tpie/stream.h>
#include <set>

using namespace terrastream;
using namespace tpie::ami;

/////////////////////////////////////////////////////////
///
///  Methods for outputting contours to the xfig (.fig) file format.
///
/////////////////////////////////////////////////////////
namespace shape {
  /////////////////////////////////////////////////////////
  ///
  ///  Writes the contours to shape files.
  ///  path_name is the name of the file to be written.
  ///  (If the file already exists, it is overwritten)
  ///  topo_stream contains the topology of the unsimplified stream
  ///  (The simplified stream is assumed to be a subset)
  ///	 unsimplified_stream is the unsimplified stream of contour segments.
  ///     These segments will be drawn black if info is null.
  ///	 simplified_stream is the simplified stream of contour segments.
  ///     These segments will be drawn blue if info is null or not at all if
  ///     simplified_stream is null.
  ///
  /////////////////////////////////////////////////////////
  void to_shape(char const* const file_suffix,
				stream<contour_point>& unsimplified_stream,
				stream<topology_edge> &topology,				   
				stream<contour_point>* const simplified_stream,
				float contour_interval, float e_z);
}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_CONTOUR_TO_SHAPE_H__*/
