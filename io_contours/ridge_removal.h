// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_RIDGE_REMOVAL_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_RIDGE_REMOVAL_H__
#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include "contour_types.h"
#include <terrastream/common/tflow_types.h>

namespace terrastream{

	//Remove all ridges and duplets in the stream of signed_contour_segments. This assumes that the
	//sign indicates whether a segment lying on the edge of a triangle was generated as the lowest
	//or highest part of the triangle (the output from the intersect routine in intersect.h).
	//The segments are outputted with the sign set to false.
	void remove_ridges_and_duplets(stream<signed_contour_segment> &input,stream<labelling_signed_contour_segment> &output);

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_RIDGE_REMOVAL_H__*/
