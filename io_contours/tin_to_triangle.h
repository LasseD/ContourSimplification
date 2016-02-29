// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_TIN_TO_TRIANGLE_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_TIN_TO_TRIANGLE_H__
#include <terrastream/common/common.h>
#include <terrastream/common/sort.h>
#include <terrastream/io_contours/contour_types.h>
#include <terrastream/common/tin_reader.h>

namespace terrastream{
	/////////////////////////////////////////////////////////
	/// \brief Converts a tin into a stream of type triangle
	///
	/// Takes as input a tin_read and converts the tin into
	/// a stream of type triangle, by sorting the tin_triangles
	/// three times and assigning the points
	///
	/// \param[in] in An opened tin reader, read the tin to convert
	/// \param[out] out Output stream for the triangles
	/////////////////////////////////////////////////////////
	void tin_to_triangle(tin_reader<elev_t> &in, stream<triangle> &out);
}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_TIN_TO_TRIANGLE_H__*/
