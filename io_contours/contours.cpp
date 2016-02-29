// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "contours.h"
#include "intersect.h"
#include "ridge_removal.h"
#include "outer_curves.h"
#include <terrastream/common/labelling.h>
#include "cw_ordering.h"
#include "topology.h"

//#define DEBUG_CONTOURS

using namespace std;
using namespace terrastream;
using namespace tpie;

typedef signed_contour_segment ss;
typedef labelled_signed_contour_segment lss;
typedef ranked_labelled_signed_contour_segment rlss;
typedef topology_edge topo;
typedef contour_point cp;

void terrastream::compute_contours(stream<triangle> &tris,
								   elev_t gran,
								   stream<cp> &out_segs,
								   stream<topo> &out_topo) {
	compute_contours(tris,gran,0.0f,out_segs,out_topo);
}

// sets x,y to the common point of own and other if set_common, 
// else other of own.
inline void common_point(rlss &own, rlss &other, xycoord_t &x, xycoord_t &y, bool set_common) {
	if(((own.x1 == other.x1 && own.y1 == other.y1) || 
		(own.x1 == other.x2 && own.y1 == other.y2)) == set_common) {
		x = own.x1;
		y = own.y1;
	}
	else {
		x = own.x2;
		y = own.y2;
	}
}

void build_small_output(stream<rlss> &large,stream<cp> &small) {
	rlss prev, prevprev; // labels are -1 when initialized here.
	xycoord_t x,y;

	large.seek(0);
	rlss *l;
	// first segment:
	if(large.read_item(&l) != ami::NO_ERROR)
		return;
	prev = *l;
	prevprev.label = prev.label-10;

	int rank = 0;
	contour_point first;
	while(large.read_item(&l) == ami::NO_ERROR) {
		if(l->x1 == l->x2 && l->y1 == l->y2) {
			cerr << "ERROR! POINT SEG" << endl;
			continue;
		}
//		cerr << "Adding " << *l << endl;

		if(prevprev.label != prev.label) { // new contour:
			// add extra starting node:
			rank = 0;
			common_point(prev, *l, x, y, false); // not common point.
			small.write_item(first = cp(x,y,rank++,prev.label));
		}
		if(l->label == prev.label) { // normal add.
			common_point(prev, *l, x, y, true);
			small.write_item(cp(x,y,rank++,prev.label));
		}
		else { // last point add.
			common_point(prev, prevprev, x, y, false);
			contour_point last = cp(x,y,rank++,prev.label);
//			if(!(first == last))
//				cerr << first << " @ " << last << std::endl;
			assert(first == last);
			small.write_item(last);
		}

		prevprev = prev;
		prev = *l;
	}

	// handle very last:
	common_point(prev, prevprev, x, y, false);
	small.write_item(cp(x,y,prev.rank+1,prev.label));
}

void print_labelling_segs_in_region(stream<labelling_signed_contour_segment> &segs) {
	xycoord_t min_x = -1000000000;
	xycoord_t max_x = 10000000000;
	xycoord_t min_y = -1000000000;
	xycoord_t max_y = 10000000000;

	labelling_signed_contour_segment *s;
	std::cerr << "Labelling segments in stream within [" << min_x << ";" << max_x << "]X[" << min_y << ";" << max_y << "]: " << std::endl;
	while(segs.read_item(&s) == ami::NO_ERROR) {
		if(std::min(s->p[0].x,s->p[1].x) >= min_x &&
		   std::max(s->p[0].x,s->p[1].x) <= max_x &&
		   std::min(s->p[0].y,s->p[1].y) >= min_y &&
		   std::max(s->p[0].y,s->p[1].y) <= max_y) {
			std::cerr << *s << std::endl;
		}
	}
	segs.seek(0);
}

void terrastream::compute_contours(stream<triangle> &tris,
								   elev_t gran,float z_diff,
								   stream<cp> &o_segs,
								   stream<topo> &out_topo) {
  //Prepare
  tris.seek(0);
  o_segs.truncate(0);
  out_topo.truncate(0);

  //Intersect triangulation
  
  stream<ss> segs;

  map_info inf = intersect(tris,gran,z_diff,segs);
  segs.seek(0);
//  cerr << "After intersect:" << endl;
//  print_segs_in_region(segs);

  //Remove duplicates
  stream<labelling_signed_contour_segment> no_duplets;
  remove_ridges_and_duplets(segs,no_duplets);
  segs.truncate(0);

  cerr << "#seg After remove ridges and duplets: " << no_duplets.stream_len() << endl;
//  print_labelling_segs_in_region(no_duplets);

  //Add the outer contours
  add_outer_curves(tris,gran,inf,no_duplets);
  no_duplets.seek(0);

  cerr << "#segs after add outer curves:" << no_duplets.stream_len() << endl;
  //print_labelling_segs_in_region(no_duplets);

  //Label the segments
  stream<lss> lbl;
  Labelling<labelling_signed_contour_segment, xycoord_t> l;
  
  tflow_progress dummy("DUMMY", "DUMMY:",0,1,1);

  l.label<lss, DummyLabel<segment_point>, DummyCH>(no_duplets, &lbl, NULL, NULL, dummy);

  lbl.seek(0);
//  cerr << "After labelling:" << endl;
//  print_segs_in_region(lbl);

  assert(assert_eq(lbl, no_duplets));
  no_duplets.truncate(0);

#ifdef DEBUG_CONTOURS
  cerr << "cc- cw ordering" << endl;
#endif
  //Order the contours
  stream<rlss> out_segs;
  cw_order(lbl,out_segs);
  //assert(assert_eq(lbl, out_segs)); // FAIL: Will remove all contours wo. real segs (sign 0)
  lbl.truncate(0);
  out_segs.seek(0);

#ifdef DEBUG_CONTOURS
  cerr << "After cw_ordering:" << endl;
  print_segs_in_region(out_segs);

  cerr << "cc- build topo" << endl;
#endif
  //Compute topology
  stream<rlss> out_segs2;
  relabel(out_segs,out_segs2);
  out_segs2.seek(0);
  out_segs.truncate(0);
  build_topology(out_segs2,out_topo); // This assumes that input is sorted by label,rank!

  stream<cp> o_segs2;
  build_small_output(out_segs2,o_segs2);  // Minimize output: (TODO FIXME: Merge this upwards)
  out_segs2.truncate(0);
  out_topo.seek(0);
  o_segs2.seek(0);

  // make real output:
  order_for_simplification(out_topo, o_segs2, o_segs);
  o_segs2.truncate(0);
}
