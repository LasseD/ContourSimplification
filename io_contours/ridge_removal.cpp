// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "ridge_removal.h"
#include <terrastream/common/sort.h>

using namespace std;
using namespace tpie;
using namespace terrastream;

struct signed_segment_sorter{
	int compare(const signed_contour_segment &s1,const signed_contour_segment &s2) const{
		if (s1.x1!=s2.x1) return s1.x1<s2.x1 ? -1 : 1;
		if (s1.y1!=s2.y1) return s1.y1<s2.y1 ? -1 : 1;
		if (s1.x2!=s2.x2) return s1.x2<s2.x2 ? -1 : 1;
		if (s1.y2!=s2.y2) return s1.y2<s2.y2 ? -1 : 1;
		return 0;
	}
};

void terrastream::remove_ridges_and_duplets(stream<signed_contour_segment> &in,stream<labelling_signed_contour_segment> &out) {
	//Sort the input stream such that identical segments occur consecutively
	if (in.stream_len()==0) return;
	in.seek(0);
	tflow_progress sort_progress("Sorting segments", "Sorting segments", 0, in.stream_len(), 1);
	signed_segment_sorter cmp;
	ts_sort(&in,&cmp,&sort_progress);
	in.seek(0);
	//Start the removal of duplicates and ridges
	signed_contour_segment s, *t;
	in.read_item(&t);
	s = *t;
	tflow_progress progress("Removing ridges", "Removing ridges", 0, in.stream_len(), 1);
	while (true) {
		bool done=true;
		bool write_it=true;
		bool alone=true;
		while (in.read_item(&t)==ami::NO_ERROR) {
			progress.step();
			//Check if s and t differ
			if (t->x1!=s.x1 || t->x2!=s.x2 || t->y1!=s.y1 || t->y2!=s.y2) {
				done=false;
				break;
			}
			alone=false;
			//The segment t is a duplicate of s
			if (t->sign==s.sign) {
				//This is a ridge
				write_it=false;
			}
		}
		//A lonely segment lies on the boundary if it is signed, and will have same sign if generated from outside.
		if (alone && s.sign) write_it=false;
		//Check if the segment s should be written
		if (write_it) {
			//Write the segment s
			labelling_signed_contour_segment o(s.x1,s.y1,s.x2,s.y2,s.z);
			out.write_item(o);
		}
		//Check if we are done with the stream
		if (done) break;
		//Prepare for next segment
		s=*t;
	}
	progress.done();
}
