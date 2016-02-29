// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CW_ORDERING_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CW_ORDERING_H__
#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include "contour_types.h"

namespace terrastream{
  void print_stream_status(stream<labelled_signed_contour_segment> &in);

  /* Takes the input stream of labelled_signed_contour_segments and sorts it by their labels, and
   * secondarely sorts each contour into cw order. The first segment in the contour is the left-most
   * contour pointing the most upwards.
   * This step also removes the last type of degeneracy, namely contours intersecting in a single
   * point. Such contours are split into separate contours with their own labels. Furtermore, all
   * contours consisting completely of outer segments will be removed.
   */
  void cw_order(stream<labelled_signed_contour_segment> &ls,stream<ranked_labelled_signed_contour_segment> &out);

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CW_ORDERING_H__*/
