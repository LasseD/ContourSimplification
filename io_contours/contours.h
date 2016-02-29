// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CONTOURS_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CONTOURS_H__

#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include "contour_types.h"
#include <terrastream/common/tflow_types.h>
#include <terrastream/common/wlabel.h>
#include <terrastream/common/labelling.h>
#include <set>

namespace terrastream{

template<typename T>
void print_segs_in_region(stream<T> &segs) {
	xycoord_t min_x = -1000000000;
	xycoord_t max_x = 10000000000;
	xycoord_t min_y = -1000000000;
	xycoord_t max_y = 10000000000;

	T *s;
	std::cerr << "Segments in stream within [" << min_x << ";" << max_x << "]X[" << min_y << ";" << max_y << "]: " << std::endl;
	while(segs.read_item(&s) == ami::NO_ERROR) {
		if(std::min(s->x1,s->x2) >= min_x &&
		   std::max(s->x1,s->x2) <= max_x &&
		   std::min(s->y1,s->y2) >= min_y &&
		   std::max(s->y1,s->y2) <= max_y){// && s->label<=1 && s->x1!=s->x2 && s->y1!=s->y2 ) {
			std::cerr << *s << std::endl;
		}
	}
	segs.seek(0);
}

template<typename A>
bool assert_eq(stream<A> &a, stream<labelling_signed_contour_segment> &b) {
	if(a.stream_len() == b.stream_len()) {
		std::cerr << "Equal stream lengths, a:" << a.stream_len() << ", b: " << b.stream_len() << std::endl;
		return true;
	}
	return false;
}

template<typename A, typename B>
bool assert_eq(stream<A> &a, stream<B> &b) {
	if(a.stream_len() == b.stream_len()) {
		std::cerr << "Equal stream lengths, a:" << a.stream_len() << ", b: " << b.stream_len() << std::endl;
		return true;
	}
	std::cerr << "Assertion failed, a:" << a.stream_len() << ", b: " << b.stream_len() << std::endl;
	std::cerr << "Not common:" << std::endl;
	a.seek(0);
	b.seek(0);
	
	std::set<contour_segment> a_set;
	int i = 0;

	A* ele;
	while(a.read_item(&ele) == NO_ERROR) {
		contour_segment cs = (contour_segment)(*ele);
		bool inserted = a_set.insert(cs).second;
		if(!inserted) {
			std::cerr << "Duplicate in a: " << *ele << " (" << ++i << ")" << std::endl;			
		}
	}
	B* eleb;
	while(b.read_item(&eleb) == NO_ERROR) {
		contour_segment csb((contour_segment)(*eleb));
		std::pair<std::set<contour_segment>::iterator, bool> p = a_set.insert(csb);
		if(p.second) {
			std::cerr << "In b but not a: " << *eleb << " (" << ++i << ")" << std::endl;			
		}
		a_set.erase(p.first);
	}
	for(std::set<contour_segment>::iterator it = a_set.begin(); it != a_set.end(); ++it) {
		std::cerr << "In a but not b: " << *it << " (" << ++i << ")" << std::endl;					
	}	
	return false;
}

void compute_contours(stream<triangle> &triangulation,
					  elev_t granularity,
					  stream<contour_point> &out_segs,
					  stream<topology_edge> &out_topo);

// Computes contours with heights of the given granularity, and for every height t, contours for heights
//t-z_diff and t+z_diff are added as well.
void compute_contours(stream<triangle> &triangulation,
					  elev_t granularity,
					  float z_diff,
					  stream<contour_point> &out_segs,
					  stream<topology_edge> &out_topo);

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CONTOURS_H__*/
