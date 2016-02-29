// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_CONTOUR_READER_H__
#define __TEST_CONTOUR_SIMPLIFICATION_CONTOUR_READER_H__
//#include "app_config.h"
#include "io_contours/contour_types.h"
#include <tpie/portability.h>
#include <cstdlib>
#include <tpie/stream.h>
#include <terrastream/common/grid_reader.h>

using namespace terrastream;
using namespace tpie::ami;

typedef elev_t height_type;

namespace contour_reader {
void read_grid(grid_reader<height_type> &reader,
			   float const contour_interval,
			   float const e_z,
			   stream<topology_edge>& os_topo,
			   stream<contour_point>& os_segs);
}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_CONTOUR_READER_H__*/
