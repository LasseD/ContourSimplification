// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_OUTER_CURVES_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_OUTER_CURVES_H__
#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include "contour_types.h"
#include <terrastream/common/tflow_types.h>

namespace terrastream{

	//This routine adds curves outside the convex triangulation.
	//Each segment in the output has their sign set to true.
	//The map_info object can be obtained by the intersect routine in intersect.h
	void add_outer_curves(stream<triangle> &tris,elev_t granularity,map_info &info,stream<labelling_signed_contour_segment> &output);

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_OUTER_CURVES_H__*/
