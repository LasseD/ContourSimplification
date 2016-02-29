// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_CONTOUR_SIMPLIFICATION_H__
#define __TEST_CONTOUR_SIMPLIFICATION_CONTOUR_SIMPLIFICATION_H__
#include "io_contours/contour_types.h"
#include "decomposition.h"
#include <tpie/stream.h>
#include <set>

using namespace tpie::ami;
using namespace terrastream;
using namespace simplification;

typedef triangulated_ranked_labelled_signed_contour_segment ts;
typedef ranked_labelled_signed_contour_segment rlss;
typedef topology_edge topo;

/////////////////////////////////////////////////////////
///
///  Various methods for simplifying contours
///
/////////////////////////////////////////////////////////
namespace simplification {

	/////////////////////////////////////////////////////////
	///
	///  Constrained Douglas Peucker's algorithm for polygonal line simplification.
	///  The algorithm follows the simple original proposal with a running time between O(n) and O(n^2) (typically O(nlogn)).
	///  e_simplify is the allowed error margin
	///  This algorithm assumes the line segments of the contours are sorted and every contour forms a cycle.
	///
	/////////////////////////////////////////////////////////
	void constrained_dp(const float e_simplify,
						stream<contour_point> &input_segments,
						stream<topo> &topology,
						elev_t granularity, float e_granularity,
						stream<contour_point> &output);
}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_CONTOUR_SIMPLIFICATION_H__*/
