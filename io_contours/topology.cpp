// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "topology.h"
#include <terrastream/common/sort.h>
#include <tpie/priority_queue.h>
#include <tpie/queue.h>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <iostream>

//#define DEBUG_SWEEP
//#define DEBUG_OFS
//#define DEBUG_RL

using namespace terrastream;
using namespace tpie;

typedef std::pair<xycoord_t,xycoord_t> pt;
typedef ranked_labelled_signed_contour_segment rlss;
typedef topology_edge topo;
typedef std::pair<int,int> int_int;

struct label_point {
	xycoord_t x, y;
	int label;

	label_point(xycoord_t _x,xycoord_t _y,int l) : x(_x),y(_y),label(l) {}
	label_point(const label_point &lp) : x(lp.x),y(lp.y),label(lp.label) {}
	label_point() {}

	const bool operator==(const label_point &l) const {
		return x==l.x&&y==l.y;
	}

	bool operator<(const label_point &l) const {
		if(x != l.x) {
			return x < l.x;
		}
		if(y != l.y) {
			return y < l.y;
		}		
		return label < l.label;
	}

	friend std::ostream& operator << (std::ostream &ostr, const label_point &l) {
		return ostr << "(" << l.x << "," << l.y << "," << l.label << ")";
	}
};

struct relabel_cmp {
	int compare(const label_point &a,const label_point &b) const{
		if(a < b) 
			return -1;
		if(b < a) 
			return 1;
		return 0;
	}
};
struct first_cmp {
	bool operator() (const int_int& a, const int_int& b) const {
		return compare(a,b) < 0;
	}
	int compare(const int_int &a,const int_int &b) const {
		if(a.first != b.first)
			return a.first - b.first;
		return a.second - b.second;
	}
};
struct label_rank_cmp {
	int compare(const rlss &a,const rlss &b) const{
		if(a.label < b.label || (a.label == b.label && a.rank < b.rank)) 
			return -1;
		if(a.label > b.label || (a.label == b.label && a.rank > b.rank)) 
			return 1;
		return 0;
	}
};

void terrastream::relabel(stream<rlss> &ls, stream<rlss> &out) {
	ls.seek(0);

#ifdef DEBUG_RL
	std::cerr << "Relabel start from: " << std::endl;
	print(ls);   
#endif

	// Find min points of all contours:
	stream<label_point> label_points; // to sort for relabeling.
	int prev_label = -1000; // no such label.
	rlss *s;
	while(ls.read_item(&s) == NO_ERROR) {
		if(s->label != prev_label) {
			assert(s->label > prev_label); // increasing from cw_ordering.
			// Read contour and find min:
			label_point min_point(s->x1,s->y1,s->label);
			err ae;
			while((ae = ls.read_item(&s)) == NO_ERROR && s->label == prev_label) {
				if(s->x1 < min_point.x) {
					min_point.x = s->x1;
					min_point.y = s->y1;
				}
			} 
			label_points.write_item(min_point);
			if(ae != NO_ERROR)
				break;
		}
		prev_label = s->label;
	}
	label_points.seek(0);
	ls.seek(0);

#ifdef DEBUG_RL
	std::cerr << "Relabel points: " << std::endl;
	print(label_points);   
#endif

	// sort:
	relabel_cmp cmp;
	ts_sort(&label_points,&cmp);
	label_points.seek(0);

#ifdef DEBUG_RL
	std::cerr << "Relabel points sorted: " << std::endl;
	print(label_points);   
#endif

	// Set new labels:
	stream<int_int> new_labels;
	label_point *lp;
	int i=0;
	while(label_points.read_item(&lp) == NO_ERROR) {
		new_labels.write_item(int_int(lp->label,i++));
	}	
	new_labels.seek(0);
	label_points.truncate(0);
	
	// sort on original labels:
	first_cmp cmp2;
	ts_sort(&new_labels,&cmp2);
	new_labels.seek(0);

	// Augment new labels:
	err e_segs = ls.read_item(&s);
	int_int *nid;
	while(new_labels.read_item(&nid) == NO_ERROR) {
		int contour = nid->first;
		while(e_segs == NO_ERROR && s->label <= contour) {
			if(s->label == contour) {
				rlss p(*s);
				p.label = nid->second;
				out.write_item(p);
			}
			else {
				std::cerr << "Warning: AX3  Ignoring seg: " << *s << std::endl;
			}
			e_segs = ls.read_item(&s);
		}
	}
	ls.seek(0);	
	out.seek(0);
	new_labels.truncate(0);

	// sort on new labels:
	label_rank_cmp cmp3;
	ts_sort(&out,&cmp3);
	out.seek(0);

#ifdef DEBUG_RL
	std::cerr << "Relabel output: " << std::endl;
	print(out);   
#endif
}

struct faced_ranked_labelled_signed_contour_segment : rlss{
	bool outside_up;
	faced_ranked_labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,
												 bool sign,int label,int rank,bool fce) :
		rlss(x1,y1,x2,y2,z,sign,label,rank),outside_up(fce) {}

	faced_ranked_labelled_signed_contour_segment(rlss &r,bool fce) :
		rlss(r.x1,r.y1,r.x2,r.y2,r.z,r.sign,r.label,r.rank) , outside_up(fce) {}

	faced_ranked_labelled_signed_contour_segment() {}

	const bool operator==(const faced_ranked_labelled_signed_contour_segment &l) const { // for Sweep line
		return (p1() == l.p1() && p2() == l.p2()) ||
			(p1() == l.p2() && p2() == l.p1());
	}

	bool operator<(const faced_ranked_labelled_signed_contour_segment &l) const {
		if (x1==l.x1) {
			if (y1==l.y1) {
				if (x2==l.x2) {
					return y2 < l.y2;
				} else
					return x2 < l.x2;
			} else
				return y1 < l.y1;
		} else
			return x1 < l.x1;
	}

	friend std::ostream& operator << (std::ostream &ostr, const faced_ranked_labelled_signed_contour_segment &l) {
		return ostr << l.rank << "@" << l.label << " (" << l.x1 << "," << l.y1 << ") - (" << l.x2 << "," << l.y2 << ") z=" << l.z <<
			" sign=" << l.sign << (l.outside_up ? " outside up" : " inside up");
	}
};

typedef faced_ranked_labelled_signed_contour_segment frlss;

inline bool is_vertical(rlss &s) {
	return s.x1 == s.x2;
}

void augment_contour(std::vector<rlss> &con,stream<frlss> &out) {
	size_t i=0;
	while (is_vertical(con[i])) //Ignore first vertical segments
		i++;
	//Because of ordering, the first segment will always have outside as upwards face
	bool outside_up=true;
	out.write_item(frlss(con[i++],true));
	int prev_j;
	for (;i<con.size();i++) {
		if (is_vertical(con[i]) && is_vertical(con[i-1]))
			continue;
		//Find common point
		int j=0;
		int k;
		for (;j<2;j++) {
			pt p1 = j==0 ? pt(con[i-1].x1,con[i-1].y1) : pt(con[i-1].x2,con[i-1].y2);
			bool found=false;
			for (k=0;k<2;k++) {
				pt p2 = k==0 ? pt(con[i].x1,con[i].y1) :
					pt(con[i].x2,con[i].y2);
				if (p1==p2) {
					found=true;
					break;
				}
			}
			if (found) break;
		}
		if (is_vertical(con[i])) {
			prev_j = j;
		} else {
			if (is_vertical(con[i-1])) //if previous segment was vertical, compare to orientation of last seen non-vertical segment instead.
				j=prev_j;
			//If same side, upwards face changes
			if (j==k) outside_up = !outside_up;
			//Write the segment
			out.write_item(frlss(con[i],outside_up));
		}
	}
}

void augment_with_faces(stream<rlss> &in,stream<frlss> &out) {
	std::vector<rlss> contour;
	//Read in the single contours
	rlss* r;
	tflow_progress progress("Augmenting contours", "Augmenting contours", 0, in.stream_len(), 1);
	err ae=in.read_item(&r);
	while (true) {
		contour.clear();
		contour.push_back(*r);
		while ((ae=in.read_item(&r))==NO_ERROR && r->label==contour[0].label) {
			progress.step();
			contour.push_back(*r);
		}
		//We have extracted the single contour
		augment_contour(contour,out);
		//Check if we should continue
		if (ae!=NO_ERROR) break;
	}
	progress.done();
}

xycoord_t interpolate(const frlss &l,xycoord_t x) {
	xycoord_t a = (l.y2-l.y1)/(l.x2-l.x1);
	xycoord_t b = l.y1-a*l.x1;
	return a*x+b;
}

struct sweeporder_cmp{
	int compare(const frlss &l1,const frlss &l2) const{
#ifdef DEBUG_SWEEP
		if(l1 == l2) {
			std::cerr << "Comparing duplicates: " << std::endl << l1 << std::endl << l2 << "\n";
			assert(false);
		}
#endif

		if (l1.x1!=l2.x1) return l1.x1<l2.x1 ? -1 : 1;
		if (l1.y1!=l2.y1) return l2.y1<l1.y1 ? -1 : 1; //Topmost first
		return (interpolate(l2, l1.x2+1.0) < interpolate(l1,l1.x2+1.0)) ? -1 : 1; // TODO: Use determinants!
        //return l1.rank-l2.rank; //Make sure we get the topmost segment first when shooting*/
	}
};

struct sweepline_cmp{
	xycoord_t &x;

	sweepline_cmp(xycoord_t &x1) : x(x1) {}

	inline xycoord_t assign_y(const frlss &l) const{
		// TODO: DON'T MAKE NEW POINTS!
		return l.x1 == x ? l.y1 : (l.x2 == x ? l.y2 : interpolate(l,x));
	}

	bool operator()(const frlss &l1,const frlss &l2) const{
		if (l1==l2) return false;

/*		// this: p1,p2, l: p3,p4
		segment_point p1 = l1.p1();
		segment_point p2 = l1.p2();
		segment_point p3 = l2.p1();
		segment_point p4 = l2.p2();

		xycoord_t minx12 = std::min(p1.x,p2.x);
		xycoord_t maxx12 = std::max(p1.x,p2.x);
		xycoord_t minx34 = std::min(p3.x,p4.x);
		xycoord_t maxx34 = std::max(p3.x,p4.x);
		std::cerr << "Comparing " << l1 << " and " << l2 << std::endl;
		if(maxx12 <= minx34) {
			std::cerr << "a" << std::endl;
			return true;
		}
		if(maxx34 <= minx12) {
			std::cerr << "b" << std::endl;
			return false;
		}
		xycoord_t miny12 = std::min(p1.y,p2.y);
		xycoord_t maxy12 = std::max(p1.y,p2.y);
		xycoord_t miny34 = std::min(p3.y,p4.y);
		xycoord_t maxy34 = std::max(p3.y,p4.y);
		if(maxy12 <= miny34) {
			std::cerr << "c" << std::endl;
			return true;
		}
		if(maxy34 <= miny12) {
			std::cerr << "d" << std::endl;
			return false;
		}

		if(p1 == p3 || p1 == p4) {
			std::cerr << "e" << std::endl;
			return !p1_above_l(p2, p3, p4);
		}
		if(p2 == p3 || p2 == p4) {
			std::cerr << "f" << std::endl;
			return !p1_above_l(p1, p3, p4);
		}

		if(p1_above_l(p1, p3, p4) && p1_above_l(p2, p3, p4)) {
			std::cerr << "g" << std::endl;
			return false;
		}
		std::cerr << "h" << std::endl;
		return p1_above_l(p3, p1, p2) && p1_above_l(p4, p1, p2);
//*/


		xycoord_t y1 = this->assign_y(l1);
		xycoord_t y2 = this->assign_y(l2);

		//Assume segments are non-vertical.. can't assume that any more.
		assert(l1.x1 != l1.x2);
		assert(l2.x1 != l2.x2);
/*		if(l1.x1 == l1.x2) {
			assert(l1.x1 == x);
			if (l1.y1!=y2) 
				return l1.y1<y2;			
			assert(l1.y2!=y2);
			return l1.y2<y2;			
		}
		if(l2.x1 == l2.x2) {
			assert(l2.x1 == x);
			if (l2.y1!=y1) 
				return y1<l2.y1;			
			assert(l2.y2!=y1);
			return y1<l2.y2;			
			}*/

		if (y1!=y2) return y1<y2;
		if (x == l1.x1 && x == l2.x1) //Segments have same left endpoint
			return interpolate(l1, x+1.0)<interpolate(l2, x+1.0);
// new:		
		else if(x == l1.x2 && l1.x2 == l2.x2) // same right end point:
			return interpolate(l1, x-1.0)<interpolate(l2, x-1.0);
		else if(l1.x2 == l2.x1)
			return true;
		else {
			assert(l1.x1 == l2.x2);
			return false;
		}

//		assert(l1.x2 == l2.x2);//Segments have same right endpoint
//		return interpolate(l1, l1.x1-1.0)<interpolate(l2, l1.x1-1.0);
	}
};

struct extract_cmp{
	//Reverse the comparisons for max pq
	bool operator()(const frlss &l1,const frlss &l2) const{
		return l2.x2<l1.x2;
	}
};

void sweep(stream<frlss> &in,stream<topo> &out) {
	//Sort into sweep order
	sweeporder_cmp sweep_order;
	ts_sort(&in,&sweep_order); // left to right.
	in.seek(0);

	//Setup the sweepline structure
	xycoord_t x=0;
	sweepline_cmp cmp(x);
	std::set<frlss,sweepline_cmp> sweepline(cmp);
	std::priority_queue<frlss,std::vector<frlss>,extract_cmp> extract_pq;
	std::map<int,int> active_labels;
	std::map<int,int> parents_assigned;
	std::vector<frlss> add_to_sweepline; // buffer to be added for each x.

	//Do the sweep
	frlss *f;

#ifdef DEBUG_SWEEP
	std::cerr << "Checking ordered in stream:" << std::endl;
	xycoord_t prev_x=-10000;
	while (in.read_item(&f)==NO_ERROR) {
		assert(prev_x <= f->x1);
		prev_x = f->x1;
	}
	in.seek(0);
	std::cerr << "OK" << std::endl;
#endif

	err ae = in.read_item(&f);
	tflow_progress sweep_progress("Sweeping", "Sweeping", 0, in.stream_len(), 1);
	while (true) {
		//Read all segments beginning now
		add_to_sweepline.clear();
		xycoord_t cur_x = f->x1;
		
		if(!extract_pq.empty() && (extract_pq.top().x2<cur_x || ae != NO_ERROR)) {
			cur_x = extract_pq.top().x2;
#ifdef DEBUG_SWEEP
			std::cerr << "Handling pq x=" << cur_x << std::endl;
#endif
		}
		else {
#ifdef DEBUG_SWEEP
			std::cerr << "Handling normal x=" << cur_x << std::endl;
#endif
			add_to_sweepline.push_back(*f);
			while ((ae = in.read_item(&f))==NO_ERROR && f->x1==cur_x) {
				add_to_sweepline.push_back(*f);
				sweep_progress.step();
			}
		}

		//Update active_labels (count on sweepline)
		for (size_t i=0;i<add_to_sweepline.size();i++) {
			int lbl = add_to_sweepline[i].label;
			if (active_labels.count(lbl)==0)
				active_labels[lbl]=1;
			else
				active_labels[lbl]++;
		}

		//Insert segments beginning now
		x=cur_x;
		for (size_t i=0;i<add_to_sweepline.size();i++) {
			frlss &f2 = add_to_sweepline[i];
//			if (parents_assigned.count(f.label)==0) {
#ifdef DEBUG_SWEEP
			std::cerr << "Shooting up from " << f2 << "\n";
#endif
			std::set<frlss,sweepline_cmp>::iterator above = sweepline.upper_bound(f2);
			
//			xycoord_t fy = x == f.x1 ? f.y1 : f.y2;
			while(above != sweepline.end() && 
				  f2.outside_up && 
				  above->label != f2.label &&
				  std::max(above->x1, above->x2) == x) {
#ifdef DEBUG_SWEEP
				std::cerr << " Walking up from fail point " << *above << "\n";
#endif
				++above;
			}//*/
			
		    int lbl;
			if (above==sweepline.end()) { //Outside is my upwards face
#ifdef DEBUG_SWEEP
				std::cerr << "Hit external face" << std::endl;
#endif
				lbl = -1;
			} else {
				//I hit another segment
				frlss a = *above;
				if(!f2.outside_up || a.label == f2.label) {
					//Insert into sweepline
#ifdef DEBUG_SWEEP
					std::cerr << " Inserting to sweep line" << std::endl;
#endif
					sweepline.insert(f2);
					//Add to extraction queue
					extract_pq.push(f2);
					continue; // WTF!?
				}
				
				if (a.outside_up) {
					//I hit the inside of a contour
					lbl = a.label;
#ifdef DEBUG_SWEEP
					std::cerr << "Hit inside of parent " << a << "\n";
#endif
				} else {
#ifdef DEBUG_SWEEP
					std::cerr << "Hit sibling " << a << "\n";
#endif
					//I hit the outside of a contour, we have the same parent
					if(parents_assigned.count(a.label)==0) {
						std::cerr << "WARNING: orphan a=" << a << ", f=" << f2 << std::endl;
						sweepline.insert(f2);
						extract_pq.push(f2);
						continue;							
					}
					lbl = parents_assigned[a.label];
				}
			}

			if(parents_assigned.count(f2.label)>0) {
				if(parents_assigned[f2.label] != lbl) {
					// Output debug info:
					std::cerr << "Error: " << parents_assigned[f2.label] << "!= " << lbl << std::endl;
					std::cerr << "f:" << f2 << std::endl;
					std::cerr << "i:" << i << ",x=" << x << ", above:" << *above << std::endl;
					std::cerr << "sweep line:" << std::endl;
					for(std::set<frlss,sweepline_cmp>::iterator it = sweepline.begin(); it != sweepline.end(); ++it) {
						std::cerr << *it << std::endl;
					}
				}//*/

				assert(parents_assigned[f2.label]==lbl);
			}
			else {
				parents_assigned[f2.label]=lbl; // maintain only these 2 lines for release! (and uncomment the if...)
				assert(f2.label > lbl);
				out.write_item(topo(f2.label,lbl,f2.z));
			}
			
//			}
			//Insert into sweepline
#ifdef DEBUG_SWEEP
				std::cerr << " Inserting b " << f2 << std::endl;
#endif
			bool inserted = sweepline.insert(f2).second;
			assert(inserted);
			//Add to extraction queue
			extract_pq.push(f2);
		}

		//Remove segments ending now
		while (!extract_pq.empty() && extract_pq.top().x2<=cur_x) {
			frlss t = extract_pq.top();
			extract_pq.pop();
			bool erased = sweepline.erase(t);
#ifdef DEBUG_SWEEP
			if(!erased) {
				std::cerr << "Extracting " << t << "\n";
				std::cerr << "Sweepline is:\n";
				for(std::set<frlss,sweepline_cmp>::iterator it = sweepline.begin(); it!=sweepline.end(); ++it) { 
					std::cerr << " " << *it << "\n";
				}
			}
#endif

			assert(erased);
			int remain = --active_labels[t.label];
			if (remain==0) {
				//Clear memory for contours no longer intersecting sweepline
				//cout << "Label " << t.label << " is no longer on sweepline\n";
				active_labels.erase(t.label);
				parents_assigned.erase(t.label);				
			}
		}

		//Advance to next x location
		if(ae!=NO_ERROR && extract_pq.empty()) 
			break;
	}
	sweep_progress.done();
}

void terrastream::build_topology(stream<rlss> &in,stream<topo> &out) {
	in.seek(0);
	stream<frlss> faced;
	augment_with_faces(in,faced);
	faced.seek(0);

	sweep(faced,out);
}


pq_entry::pq_entry(int _p, int _lv, int _c) : p(_p), lv(_lv), c(_c) {}
pq_entry::pq_entry() {}
pq_entry::pq_entry(const pq_entry &p) : p(p.p), lv(p.lv), c(p.c) {}

struct level_orderer {
	int compare(const pq_entry &a,const pq_entry &b) const {
		if(a.lv != b.lv)
			return a.lv - b.lv;
		if(a.p != b.p)
			return a.p - b.p;
		return a.c - b.c;
	}
	bool operator() (const pq_entry& a, const pq_entry& b) const {
		return compare(a,b) < 0;
	}
};
struct basic_orderer {
	bool operator() (const pq_entry& a, const pq_entry& b) const {
		return a < b;
	}
	int compare(const pq_entry &a,const pq_entry &b) const{
		if(a < b)
			return -1;
		if(b < a)
			return 1;
		return 0;
	}
};
struct topo_parent_order {
	int compare(const topo &a,const topo &b) const{
		if(a.p != b.p)
			return a.p - b.p;
		return a.c - b.c;
	}
};
struct topo_child_order {
	int compare(const topo &a,const topo &b) const{
		if(a.c != b.c)
			return a.c - b.c;
		return a.p - b.p;
	}
};
struct seg_order {
	int compare(const contour_point &a,const contour_point &b) const{
		if(a.label != b.label)
			return a.label - b.label;
		return a.rank - b.rank;
	}
};

void print_stream(stream<contour_point> &s) {
	std::cerr << "|segs| = " << s.stream_len() << std::endl;	
	contour_point *t, prev;	
	while (s.read_item(&t) == NO_ERROR) {
		if(t->label != prev.label)
			std::cerr << *t << std::endl;
		prev = *t;
	}
	s.seek(0);
}

void terrastream::order_for_simplification(stream<topo> &topology,
										   stream<contour_point> &segs, 
										   stream<contour_point> &segments2) {
	// Setup:
	topology.seek(0);
	segs.seek(0);
	topo_child_order topo_child_orderer;

#ifdef DEBUG_OFS
	std::cerr << "TOPOLOGY from start:" << std::endl;
	topo *td;
	while (topology.read_item(&td) == NO_ERROR) {
		assert(td->p < td->c);
		std::cerr << *td << std::endl;
	}
	topology.seek(0);

	std::cerr << "SEGS from start:" << std::endl;
	print_stream(segs);
#endif

	// "Time forward": See paper.
	std::cerr << "Time forwarding topo tree for BFS labels" << std::endl;
	std::cout << "Time forwarding topo tree for BFS labels" << std::endl;

	stream<pq_entry> levels;
	ami::priority_queue<pq_entry> pq(0.9); 
	// insert all:
	const int INF = 1999999999; // Level 1 999 999 999 is inf.
	topo *t;
	while (topology.read_item(&t) == NO_ERROR) {
		assert(t->p < t->c);
		pq.push(pq_entry(t->p,INF,t->c));
	}
	topology.seek(0);	

    // Do time forward:
	int level = 0;
	while(!pq.empty()) {
		pq_entry e = pq.top();
		pq.pop();
		if(e.lv != INF) { // Set level.
			level = e.lv;
		}
		else {
			levels.write_item(pq_entry(e.p, level, e.c)); // output self
			pq.push(pq_entry(e.c, level+1, -2)); // Set level marker for child. -2 => Before other entries of c.
		}
	}
	levels.seek(0);

	std::cout << " Levels created." << std::endl;
#ifdef DEBUG_OFS
	std::cerr << "LEVELS created:" << std::endl;
	pq_entry *td3;
	while (levels.read_item(&td3) == NO_ERROR)
		std::cerr << *td3 << std::endl;
	levels.seek(0);
#endif

	// sort levels by level => BFS (new ids).
	std::cout << "Time forwarding level sort" << std::endl;
	std::cerr << "Time forwarding level sort" << std::endl;
	level_orderer lv_orderer;
	ts_sort(&levels,&lv_orderer);
	levels.seek(0);

#ifdef DEBUG_OFS
	std::cerr << "LEVELS sorted:" << std::endl;
	while (levels.read_item(&td3) == NO_ERROR)
		std::cerr << *td3 << std::endl;
	levels.seek(0);
#endif

	int nid_index = 0; // 10000
	stream<int_int> nids;
	nids.write_item(int_int(-1, -1)); // Special marker contour.

	// Sort every level:
	ami::queue<int_int> q; // contour,li,  enqueue(T), dequeue(**T), peek(**T)
	q.enqueue(int_int(-1,-1));
	level = 0;
	pq_entry *td4;	
	err err_lv = levels.read_item(&td4);
	const int_int *from_q = NULL;
	std::cout << " Handling levels" << std::endl;
	while(err_lv == NO_ERROR) { // handle every level:
#ifdef DEBUG_OFS
		std::cerr << "Handling level " << level << std::endl;
#endif

		int li = 0;
		expand_set<pq_entry> v(50000);
//		std::set<pq_entry> v; // p=p_li, c=self, lv=li
		// get all:
		bool first = true;
		while((first || (err_lv = levels.read_item(&td4)) == NO_ERROR) && td4->lv == level) {
#ifdef DEBUG_OFS
			std::cerr << " From level stream: " << *td4 << std::endl;
#endif
			while(from_q == NULL || from_q->first != td4->p) {
				assert(!q.is_empty());
#ifdef DEBUG_OFS
				if(from_q != NULL)
					std::cerr << "  Skipping from " << from_q->first << ", li:" << from_q->second << std::endl;
#endif
				q.dequeue(&from_q);
			}
#ifdef DEBUG_OFS
			std::cerr << " Parent: " << from_q->first << ", li: " << from_q->second << std::endl;
#endif
			assert(td4->p == from_q->first);
			
			v.insert(pq_entry(from_q->second, li++, td4->c));
			first = false;
		}
		
        // sort
		basic_orderer b_orderer;
		v.sort(b_orderer);
		//std::set<pq_entry,level_orderer> v2; // p=p_li, c=self, lv=li, order by li
		expand_set<pq_entry> v2(4000);
		pq_entry item;
		while(v.next(item)) {
#ifdef DEBUG_OFS
			std::cerr << " Out: 'pq entry':" << item << std::endl;
#endif
			nids.write_item(int_int(item.c, nid_index));
			v2.insert(pq_entry(item.c,item.c,nid_index));
			nid_index++;	
		}

		v2.sort(lv_orderer);

		while(v2.next(item)) {
#ifdef DEBUG_OFS
			std::cerr << " To Q: 'pq entry p,c'" << item << std::endl;
#endif
			q.enqueue(int_int(item.p,item.c));
		}
		level++;
	}
	
	// make id change stream:

	levels.truncate(0);
	nids.seek(0);
	std::cout << "Time forwarding level sort 2" << std::endl;
	std::cerr << "Time forwarding level sort 2" << std::endl;
	first_cmp fc;
	ts_sort(&nids,&fc);
	nids.seek(0);
	std::cout << "DONE: Time forwarding topo tree for BFS labels" << std::endl;
	std::cerr << "DONE: Time forwarding topo tree for BFS labels" << std::endl;
#ifdef DEBUG_OFS
	std::cerr << "new ids:" << std::endl;	
	int_int *td2;
	while (nids.read_item(&td2) == NO_ERROR) {
		std::cerr << td2->first << " -> " << td2->second << std::endl;
	}
	nids.seek(0);
#endif

	// Set new ids:
	// sort/scan 1:
	// - sort parent in topo-stream:
	std::cout << "Sorting topology on old parents" << std::endl;
	std::cerr << "Sorting topology on old parents" << std::endl;
	topo_parent_order topo_parent_orderer;
	ts_sort(&topology,&topo_parent_orderer);
	topology.seek(0);	

#ifdef DEBUG_OFS
	std::cerr << std::endl << "TOPOLOGY after sort on old parent:" << std::endl;
	while (topology.read_item(&td) == NO_ERROR) {
		std::cerr << *td << std::endl;
	}
	topology.seek(0);
#endif

	// - Step and replace (scan):
	contour_point *s;
	err e_topo = topology.read_item(&t);
	err e_segs = segs.read_item(&s);

	stream<topo> topo2;

 	std::cout << "Updating labels" << std::endl;
 	std::cerr << "Updating labels" << std::endl;
	int_int *nid;
	while(nids.read_item(&nid) == NO_ERROR) {
		int contour = nid->first;
		while(e_topo == NO_ERROR && t->p <= contour) {
			if(t->p == contour)
				topo2.write_item(topo(t->c, nid->second, t->c_z));
			e_topo = topology.read_item(&t);
		}
		while(e_segs == NO_ERROR && s->label <= contour) {
			if(s->label == contour) {
				contour_point p(*s);
				p.label = nid->second;
				segments2.write_item(p);
			}
#ifdef DEBUG_OFS
			else
				std::cerr << "Warning: AX2  Ignoring seg: " << *s << std::endl;
#endif
			e_segs = segs.read_item(&s);
		}
	}
	topo2.seek(0);	
	segments2.seek(0);		
	
#ifdef DEBUG_OFS
	std::cerr << "SEGS 2:" << std::endl;
	print_stream(segments2);
#endif

	// prepare out streams:
	segs.truncate(0);
	topology.truncate(0);
	// - sort segments:
	seg_order seg_orderer;
	std::cerr << "Sorting segments on new labels" << std::endl;
	ts_sort(&segments2,&seg_orderer);
	segments2.seek(0);		

	// sort/scan 2: (topo stream)
	nids.seek(0);
//	topo_child_order topo_child_orderer;
	std::cerr << "Sorting topology on old child labels" << std::endl;
	ts_sort(&topo2,&topo_child_orderer);
	topo2.seek(0);	
	e_topo = topo2.read_item(&t);
	std::cerr << "Updating topology child labels" << std::endl;
	while(nids.read_item(&nid) == NO_ERROR) {
		int contour = nid->first;
		while(e_topo == NO_ERROR && t->c <= contour) {
			if(t->c == contour)
				topology.write_item(topo(nid->second, t->p, t->c_z)); 
			e_topo = topo2.read_item(&t);
		}
	} // topology done.
	topology.seek(0);	
	std::cerr << "Sorting topology for new labels" << std::endl;
	ts_sort(&topology,&topo_child_orderer);
	topology.seek(0);	
	topo2.truncate(0);
	nids.truncate(0);

#ifdef DEBUG_OFS
	std::cerr << std::endl << "TOPOLOGY end:" << std::endl;
	while (topology.read_item(&td) == NO_ERROR) {
		std::cerr << *td << std::endl;
	}
	topology.seek(0);

	std::cerr << "SEGS2 end:" << std::endl;
	print_stream(segments2);
#endif

	topo *td10, prev(-10,-10,-10);
	while(topology.read_item(&td10) == NO_ERROR) {
		assert(prev.p <= td10->p);
		assert(prev.c < td10->c);
		prev = *td10;
	}//*/
	topology.seek(0);

	std::cerr << "DONE" << std::endl;	
}
