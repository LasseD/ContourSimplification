// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_SIMPLIFY_H__
#define __TEST_CONTOUR_SIMPLIFICATION_SIMPLIFY_H__
#include <terrastream/common/common.h>
#include "io_contours/contour_types.h"

namespace terrastream {
namespace simplification {
void run(grid_reader<elev_t> &reader, float gran, float e_z, float e_dp, bool cdp);
}
}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_SIMPLIFICATION_H__*/
