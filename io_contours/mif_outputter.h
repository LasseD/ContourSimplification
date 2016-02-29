// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_MIF_OUTPUTTER_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_MIF_OUTPUTTER_H__
#include <terrastream/common/common.h>
#include "contour_types.h"

namespace terrastream {
	////////////////////////////////////////////////////////
	/// Outputs a stream of ranked_labelled_signed_contour_segment
	/// into a .mif and .mid file
	///
	/// \param[in] stream segment stream to output
	/// \param[in] basename basename of .mif and .mid file to output
	////////////////////////////////////////////////////////
	void output_mif(stream<ranked_labelled_signed_contour_segment>& stream_in, stream<topology_edge> &topos, std::string basename);

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_MIF_OUTPUTTER_H__*/
