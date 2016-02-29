// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

//#define DEBUG_SIMPLIFICATION
//#define DEBUG_CI2
//#define DEBUG_RP
 
#include "io_contours/contour_types.h"
#include "contour_simplification.h"
#include "decomposition.h"
#include "util.h"
#include <set>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <tpie/queue.h>

using namespace std;
using namespace terrastream;
using namespace tpie;
using namespace simplification;
using namespace boost::posix_time;

typedef labelled_signed_contour_segment lss;
typedef ranked_labelled_signed_contour_segment rlss;
typedef triangulated_ranked_labelled_signed_contour_segment ts;
typedef topology_edge topo;
typedef vector<contour_point> contour;

typedef pair<xycoord_t,xycoord_t> pt;
typedef pair<contour_point,contour_point> pt2;
typedef unsigned int CW_SIZE_T; // The size_t type used throughout the code -- kept as 32-bit to avoid doubling memory usage

// distance between the line segment p1-p2 and the point a.
double distanceLineToPoint(vector<contour_point> &points, int pp1, int pp2, int pa) {
	assert(pp1 >= 0);
	assert(pp2 >= 0);
	assert(pa >= 0);
	assert(pp1 < points.size());
	assert(pp2 < points.size());
	assert(pa < points.size());

	contour_point a = points[pa];
	contour_point p1 = points[pp1];
	contour_point p2 = points[pp2];
	if (p1.x == p2.x && p1.y == p2.y) {
		return (a-p1).length();
	}
	contour_point dir = p2 - p1;
	contour_point v = a - p1;
	assert(dot(dir, dir) != 0);
	double t = (dot(v, dir) / dot(dir, dir));
	t = min(max(t, 0.0), 1.0);
	contour_point p = p1 + dir*t;
	return (a-p).length();
}

bool near(contour_point &p, contour_point &i1, contour_point &i2) {
	return i1.rank <= p.rank && p.rank <= i2.rank;
}

int print_intersection(contour &c, 
					   contour_point &i1, contour_point &i2, 
					   contour_point &j1, contour_point &j2) {
	cerr << "Contour segs" << endl;
	contour_point prev = *(c.begin());
	int cnt = 0;
	for(contour::iterator it = ++(c.begin()); it != c.end(); ++it) {
		if(near(prev, i1, j2) && near(*it, i1, j2)) {
			cerr << prev << " >> " << *it << endl;
			cnt++;
		}
		prev = *it;
	}
	return cnt;
}
void print_intersection(contour &c1, 
						contour &c2, contour_point &i1, contour_point &i2, 
						contour_point &j1, contour_point &j2) {
	if(i1.rank > i2.rank)
		swap(i1,i2);
	if(j1.rank > j2.rank)
		swap(j1,j2);
	if(i1.rank > j2.rank) {
		swap(j1,i1);
		swap(j2,i2);
	}
	cerr << "Intersection discovered. Printing involved segments" << endl;
	cerr << "Crossing lines:" << endl;
	cerr << "Color:  0x0000ff" << endl;
	cerr << i1 << " >> " << i2 << endl;
	cerr << j1 << " >> " << j2 << endl;
	cerr << "From part of simplified contour:" << endl;
	int s_cnt = print_intersection(c1, i1, i2, j1, j2);
	cerr << "Color:  0x000000" << endl;
	int o_cnt = print_intersection(c2, i1, i2, j1, j2);
	cerr << "#original: " << o_cnt << ", #simplified: " << s_cnt << endl;
//	int *a = NULL;
//	cerr << *a;
}

// Credits: http://www.math.niu.edu/~rusin/known-math/95/line_segs
bool lines_cross(contour_point const &p1,contour_point const &p2, 
				 contour_point const &p3,contour_point const &p4) {
#ifdef DEBUG_CI2
	if(start_debug(p1)) {
		cerr << " Testing cross of " << p1 << "-" << p2 << " vs " << p3 << "-" << p4 << endl;
	}
#endif
	if(abs(p1.rank-p2.rank) == 1 && abs(p3.rank-p4.rank) == 1) {
		return false; // original lines always OK.
	}

	// Overlapping points:
	assert(p1 != p2);
	assert(p3 != p4);
	if(p1 == p3 && p2 == p4) {
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << "  Overlapping a!" << endl;
		}
#endif
		return !(p1.rank == p3.rank && p2.rank == p4.rank); // overlapping segments are intersecting.
	}
	else if(p1 == p4 && p2 == p3) {
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << "  Overlapping b!" << endl;
		}
#endif
		return !(p1.rank == p4.rank && p2.rank == p3.rank); // -||-
	}

	// Bounding box exclusion:
	xycoord_t minx12 = min(p1.x,p2.x);
	xycoord_t maxx12 = max(p1.x,p2.x);
	xycoord_t minx34 = min(p3.x,p4.x);
	xycoord_t maxx34 = max(p3.x,p4.x);
	if((maxx12 <= minx34 && minx12 != maxx34) || 
	   (maxx34 <= minx12 && minx34 != maxx12)) {
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << " Cross bail x" << endl;
		}
#endif
		return false;
	}
	xycoord_t miny12 = min(p1.y,p2.y);
	xycoord_t maxy12 = max(p1.y,p2.y);
	xycoord_t miny34 = min(p3.y,p4.y);
	xycoord_t maxy34 = max(p3.y,p4.y);
	if((maxy12 < miny34 && miny12 != maxy34) || 
	   (maxy34 < miny12 && miny34 != maxy12)) {
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << " Cross bail y" << endl;
		}
#endif
		return false;
	}

	// Determinants:
	xycoord_t x1 = p1.x;
	xycoord_t y1 = p1.y;
	xycoord_t x2 = p2.x;
	xycoord_t y2 = p2.y;
	xycoord_t x3 = p3.x;
	xycoord_t y3 = p3.y;
	xycoord_t x4 = p4.x;
	xycoord_t y4 = p4.y;

	float det123 = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
	float det124 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1);
	float det341 = (x3 - x1)*(y4 - y1) - (x4 - x1)*(y3 - y1);
	float det342 = det123 - det124 + det341;
	int s123 = simplification::int_sign(det123, SQRT_EPSILON);
	int s124 = simplification::int_sign(det124, SQRT_EPSILON);
	int s341 = simplification::int_sign(det341, SQRT_EPSILON);
	int s342 = simplification::int_sign(det342, SQRT_EPSILON);

	// Do not allow lines to overlap, except at end points:
	if((s123 == 0 && !(p1 == p3 || p2 == p3)) ||
	   (s124 == 0 && !(p1 == p4 || p2 == p4)) ||
	   (s341 == 0 && !(p3 == p1 || p4 == p1)) ||
	   (s342 == 0 && !(p3 == p2 || p4 == p2))) {
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << " Spurious overlap!" << endl;
		}
#endif
		return true;
	}

	// Vertical lines:
	if(x1 == x2) { // First line vertical!
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << " sVertical 12" << endl;
		}
#endif
		if(x3 == x4) { // Other line vertical!
#ifdef DEBUG_CI2
			if(start_debug(p1)) {
				cerr << " Vertical both!" << endl;
			}
#endif
			return !(max(y3,y4) <= min(y1,y2) || max(y1,y2) <= min(y3,y4));
		}
		return s341 != s342;
	}
	else if(x3 == x4) { // Other line vertical!
#ifdef DEBUG_SLC
		if(start_debug(p1)) {
			cerr << " Vertical 34" << endl;
		}
#endif
		return s123 != s124;
	}

	// Overlapping lines:
	if(s123 == 0 && s124 == 0 || s341 == 0 && s342 == 0) {
#ifdef DEBUG_CI2
		if(start_debug(p1)) {
			cerr << "  Overlapping 12 / 34" << endl;
		}
#endif
		return true;
	}

#ifdef DEBUG_CI2
	if(start_debug(p1)) {
		cerr << " Testing cross of " << p1 << "-" << p2 << " vs " << p3 << "-" << p4 << endl;
		cerr << " det123 " << det123 << ", det124" << det124 << ", det341" << det341 << ", det342" << det342 << endl;
	}
#endif
	return det123*det124 < 0 && det341*det342 < 0;
}

bool lines_cross(lp_seg const &s1, lp_seg const &s2) {
	contour_point p1 = s1.p1();
	contour_point p2 = s1.p2();
	contour_point p3 = s2.p1();
	contour_point p4 = s2.p2();
	return lines_cross(p1, p2, p3, p4);
}
bool lines_cross(pt2 const &s1, pt2 const &s2) {
	contour_point p1 = s1.first;
	contour_point p2 = s1.second;
	contour_point p3 = s2.first;
	contour_point p4 = s2.second;
	return lines_cross(p1, p2, p3, p4);
}

int cnt_intersections = 0;
/*bool contains_intersections1(contour &c, contour &original) {
	contour_point prev1;
	bool first1 = true;
	contour_point i1, i2, j1, j2;
	bool does = false;
	for(contour::iterator it1 = c.begin(); it1 != c.end(); ++it1) {
		contour_point p1 = *it1;

		contour_point prev2;
		if(!first1) {
			bool first2 = true;
			for(contour::iterator it2 = it1; it2 != c.end(); ++it2) {
				contour_point p2 = *it2;
				if(!first2) {
					if(lines_cross(prev1,p1, prev2,p2)) {
						cerr << "Intersection discovered! Total: " << ++cnt_intersections << endl;
						does = true;
						print_intersection(c, original, prev1, p1, prev2, p2);
					}
				}
				first2 = false;
				prev2 = p2;
			}			
		}

		first1 = false;
		prev1 = p1;
	}
	return does;
}//*/

typedef set<lp_seg,sl_cmp> lpsegpset;

pt2 to_pt2(lp_seg s) {
	return pt2(s.p1(),s.p2());
}

bool update_sweep_line(lpsegpset &sweep_line, lp_seg &seg, bool is_ending, pt2 &cross1, pt2 &cross2) {
	if(!is_ending) { // starting:
#ifdef DEBUG_CI2
		if(start_debug(seg.p1())) {
			cerr << " Handling starting " << seg << endl;
		}
#endif		
		lpsegpset::iterator already = sweep_line.find(seg);
		if(already != sweep_line.end()) {
#ifdef DEBUG_CI2
			if(start_debug(seg.p1())) {
				cerr << " Crossing 0 found vs " << *already << endl;
			}
#endif		
			cross1 = to_pt2(*already);
			cross2 = to_pt2(seg);
			return true;
		}

		pair<lpsegpset::iterator,bool> pair = sweep_line.insert(seg);

		if(!pair.second) {
#ifdef DEBUG_CI2
			cerr << " Warning: not inserted, hope comparator raises! " << endl;
#endif		
			return false; // WARNING! TODO: SAFE?
		}
		assert(pair.second);

		lpsegpset::iterator it = ++(pair.first);
		if(it != sweep_line.end()) {
			lp_seg above = *it;
			if(lines_cross(above,seg)) {
#ifdef DEBUG_CI2
				if(start_debug(seg.p1())) {
					cerr << " Crossing 2 found vs " << above << endl;
				}
#endif		
				cross1 = to_pt2(above);
				cross2 = to_pt2(seg);
				return true;
			}
		}
#ifdef DEBUG_CI2
		else
			if(start_debug(seg.p1()))
				cerr << " At the end! " << seg << endl;
#endif		

		--it;
		if(it == sweep_line.begin()) {
#ifdef DEBUG_CI2
			if(start_debug(seg.p1())) {
				cerr << " At the beginning! " << seg << endl;
			}
#endif		
			return false;
		}
		lp_seg below = *(--it);
		if(lines_cross(below,seg)) {
#ifdef DEBUG_CI2
			if(start_debug(seg.p1())) {
				cerr << " Crossing 3 found vs " << below << endl;
			}
#endif		
			cross1 = to_pt2(below);
			cross2 = to_pt2(seg);
			return true;
		}
	}
	else { // is ending:
#ifdef DEBUG_CI2
		if(start_debug(seg.p1())) {
			cerr << " Handling ending " << seg << endl;
		}
#endif		
		bool erased = sweep_line.erase(seg);
		assert(erased);
		lpsegpset::iterator it = sweep_line.upper_bound(seg);
		if(it == sweep_line.end())
			return false;
		lp_seg above = *it;
		if(it == sweep_line.begin())
			return false;
		lp_seg below = *(--it);
		if(lines_cross(above,below)) {
			cross1 = to_pt2(below);
			cross2 = to_pt2(above);
			return true;
		}
	}
	return false;
}

// Also works for vertical lines: below is "left"
inline bool is_left_of_point(lp_seg &seg, link_point *point) {
	assert(point != NULL);
	link_point *other = seg.lp2();
	if(other->p == point->p) {
		other = seg.lp1();
	}
	else {
		assert(point->p == seg.p1());
	}
	return other->p.x < point->p.x || (other->p.x == point->p.x && other->p.y < point->p.y);
}

bool ptr_cmp2(link_point* a, link_point* b) {
	return *a < *b;
}

struct ep{
	pt p;
	vector<pair<CW_SIZE_T,CW_SIZE_T> > ns; // adjacent point, seg index in contour
	ep(pt _p) : p(_p) {}
};

const xycoord_t PI = acos(-1.0);

bool build_eps(vector<ep> &eps, contour &c, int &c1, int &c2) {
	map<pt,CW_SIZE_T> ep_map; // pt: pair<xycoord_t,xycoord_t>, ep_map: point->index
	for (CW_SIZE_T i=1;i<c.size();i++) {
		contour_point l = c[i];
		contour_point prev = c[i-1];
		if(l == prev) {
			c1 = i-1;
			c2 = i;
			return false;
		}
		pt p1 = pt(prev.x,prev.y);
		pt p2 = pt(l.x,l.y);
		for (int j=0;j<2;j++) {
			pt cur = j==0 ? p1 : p2;
			if (ep_map.count(cur)==0) {
				eps.push_back(ep(cur)); // ns not set.
				ep_map[cur]=eps.size()-1; // set index
			}
		}
		CW_SIZE_T ip1 = ep_map[p1];
		CW_SIZE_T ip2 = ep_map[p2];
		eps[ip1].ns.push_back(make_pair(ip2,i));
		eps[ip2].ns.push_back(make_pair(ip1,i));
	}
	return true;
}

bool is_original(contour_point &p1, contour_point &prev) {
	return p1.rank == prev.rank+1;
}

//Sort neighbors of nodes of degree > 2 (in clockwise order)
bool order_points(vector<ep> &eps, contour &c, int &c1, int &c2) {
#ifdef DEBUG_CI2
	if(start_debug()) {
		cerr << "Starting Check order around points. ||=" << eps.size() << endl;
	}
#endif
	for (CW_SIZE_T i=0;i<eps.size();i++) {
		if (eps[i].ns.size()>2) {
#ifdef DEBUG_CI2
			if(start_debug()) {
				cerr << "Checking " << i << ": " << eps[i].p.first << "," << eps[i].p.second << " of size " << eps[i].ns.size() << endl;
			}
#endif
			vector<xycoord_t> tmp; // cos, index.
			xycoord_t min_arg = 100000, max_arg = -10000000;
			for (CW_SIZE_T j=0;j<eps[i].ns.size();j++) {
				pair<CW_SIZE_T,CW_SIZE_T> n = eps[i].ns[j];
#ifdef DEBUG_CI2
				if(start_debug())
					cerr << "vs " << n.first << ": " << eps[n.first].p.first << "," << eps[n.first].p.second << endl;
#endif
				xycoord_t dx = eps[n.first].p.first-eps[i].p.first;
				xycoord_t dy = eps[n.first].p.second-eps[i].p.second;
				xycoord_t l = sqrt(dx*dx+dy*dy);
				dx/=l;
				xycoord_t acos_r = acos(dx);
				if (acos_r<0) 
					acos_r += 2*PI;
				if ((dy>0 && acos_r>PI) || (dy<0 && acos_r<PI))
					acos_r = 2*PI-acos_r;
				min_arg = min(min_arg,acos_r);
				max_arg = max(max_arg,acos_r);
				tmp.push_back(acos_r);
#ifdef DEBUG_CI2
				if(start_debug())
					cerr << " Added segment " << j << " for angle " << acos_r << endl;
#endif
			}
#ifdef DEBUG_CI2
			if(start_debug())
				cerr << "Added segments ||=" << tmp.size() << " min:" << min_arg << " vs max " << max_arg << endl;
#endif
			
			for (int j=0;j<tmp.size();j++) {
				int jj = (j-1+tmp.size())%tmp.size();
#ifdef DEBUG_CI2
				if(start_debug()) {
					cerr << "Testing from " << jj << " to " << j << ":" << tmp[jj] << " / " << tmp[j] << endl;
					cerr << " from jj max:" << (tmp[jj] == max_arg) << " to j min: " << (tmp[j] == min_arg) << endl;
				}
#endif
				if(tmp[j] > tmp[jj]) {
					if(tmp[j] == min_arg && tmp[jj] == max_arg) {
#ifdef DEBUG_CI2
						if(start_debug())
							cerr << "Max/Min bail." << endl;
#endif
						continue;
					}
				
					c1 = eps[i].ns[jj].second; 
					c2 = eps[i].ns[j].second; 
					if(is_original(c[c1],c[c1-1]) && is_original(c[c2],c[c2-1])) {
#ifdef DEBUG_CI2
						if(start_debug())
							cerr << "original bail." << endl;
#endif
						continue;
					}
					return false;
				}
			}
		}
	}
	return true;
}

bool contains_intersections2(contour &c, pt2 &cross1, pt2 &cross2) {
#ifdef DEBUG_CI2
	if(start_debug())
		cerr << "Starting Contains Intersections 2" << endl;
#endif
	if(c.size() < 4) {
#ifdef DEBUG_CI2
		if(start_debug())
			cerr << "Contour too small to contain intersections" << endl;
#endif
		return false;
	}

	// Check for cw around points:
	//Build the graph

#ifdef DEBUG_CI2
	if(start_debug())
		cerr << " cw_ordering tie in starting" << endl;
#endif
	vector<ep> eps; // ep: p=pt point, ns=vector(unref neighbour index,seg index in contour)
	int c1, c2;
	bool ok = build_eps(eps, c, c1, c2);
#ifdef DEBUG_CI2
	if(start_debug())
		cerr << " Built eps" << endl;
#endif
	if(!ok || !order_points(eps, c, c1, c2)) {
		assert(c1 > 0);
		assert(c2 > 0);
		contour_point a1 = c[c1], b1 = c[c1-1];
		contour_point a2 = c[c2], b2 = c[c2-1];
		
		cross1 = pt2(a1, b1);
		cross2 = pt2(a2, b2);

#ifdef DEBUG_CI2
		if(start_debug()) {
			cerr << "Not cw ordered: " << endl;
			cerr << "seg 1: " << c1 << ":" << a1 << "-" << b1 << endl;
			cerr << "seg 2: " << c2 << ":" << a2 << "-" << b2 << endl;
			for (CW_SIZE_T i=0;i<c.size();i++) {
				cerr << i << ":" << c[i] << endl;
			}
		}
#endif
		return true;
	}
	eps.clear();//*/
#ifdef DEBUG_CI2
	if(start_debug())
		cerr << " cw_ordering tie in clear!" << endl;
#endif

    // Add segs to event queue:
	set<link_point*,bool(*)(link_point*,link_point*)> pq(ptr_cmp2);
	link_point *first_point = decomposition::addContourToPQ(pq, c, false);
	assert(first_point != NULL);

#ifdef DEBUG_CI2
	if(start_debug()) {
		cerr << "Starting sweep on |PQ| = " << pq.size() << endl;
		cerr << "First point: " << *first_point << endl;
		print_set(pq);
		cerr << "---------------------------" << endl;
	}
#endif
	
	xycoord_t x;
	sl_cmp sweep_line_cmp(0, &x);
	lpsegpset sweep_line(sweep_line_cmp);

	for(set<link_point*,bool(*)(link_point*,link_point*)>::iterator it_pq = pq.begin(); it_pq != pq.end();) {
		link_point *p = *it_pq;
		x = p->p.x;
#ifdef DEBUG_CI2
		if(start_debug())
			cerr << " Handling point " << *p << " at " << x << endl;
#endif		

		lp_seg seg(p);
		while(true) { // Insert into left/right:
#ifdef DEBUG_CI2
			if(start_debug())
				cerr << " Updating l/r of " << *p << endl;
#endif		
			// Update l/r:
			bool is_left = is_left_of_point(seg, p);
			bool fail = update_sweep_line(sweep_line, seg, is_left, cross1, cross2);
			if(fail) {
				return true;
			}
			contour_point c1,c2,c3,c4;
			if(sweep_line_cmp.has_assertion_raised(c1, c2, c3, c4)) {
				cross1 = pair<contour_point,contour_point>(c1,c2);
				cross2 = pair<contour_point,contour_point>(c3,c4);
				return true;				
			}

			seg = lp_seg(p->next);
			is_left = is_left_of_point(seg, p);
			fail = update_sweep_line(sweep_line, seg, is_left, cross1, cross2);
			if(fail) {
				return true;
			}
			if(sweep_line_cmp.has_assertion_raised(c1, c2, c3, c4)) {
				cross1 = pair<contour_point,contour_point>(c1,c2);
				cross2 = pair<contour_point,contour_point>(c3,c4);
				return true;				
			}

			++it_pq;
			if(it_pq == pq.end() || !(*p == **it_pq)) {
				break;
			}
			p = *it_pq;
			seg = lp_seg(p);
		}

#ifdef DEBUG_CI2
		if(start_debug()) {
			cerr << " Sweep line:" << endl;
			print_set(sweep_line);
			if(!sweep_line.empty()) {
				lp_seg prev = *(sweep_line.begin());
				xycoord_t prev_y = prev.eval(x);
				for(lpsegpset::iterator it_sl = ++sweep_line.begin(); it_sl != sweep_line.end(); ++it_sl) {
					lp_seg sl_seg = *it_sl;
					assert(!lines_cross(prev, sl_seg));
					if(prev_y > sl_seg.eval(x)) {
						cerr << " SWEEP LINE INCONSISTENCY x ON: " << sl_seg << " vs prev at y " << prev_y << endl;
						assert(false);
					}
					if(!sweep_line_cmp(prev,sl_seg)) {
						cerr << " SWEEP LINE INCONSISTENCY x2: " << endl;
						cerr << "  compare prev,sl_seg: " << sweep_line_cmp(prev,sl_seg) << endl;
						cerr << "  compare sl_seg,prev: " << sweep_line_cmp(sl_seg,prev) << endl;
						cerr << prev << endl;
						cerr << sl_seg << " at y " << prev_y << endl;
						assert(false);
					}
					if(sweep_line_cmp(sl_seg,prev)) {
						cerr << " SWEEP LINE INCONSISTENCY x3 (too late): " << endl;
						cerr << "  compare prev,sl_seg: " << sweep_line_cmp(prev,sl_seg) << endl;
						cerr << "  compare sl_seg,prev: " << sweep_line_cmp(sl_seg,prev) << endl;
						cerr << prev << endl;
						cerr << sl_seg << " at y " << prev_y << endl;
						assert(false);
					}
					prev_y = sl_seg.eval(x);
					prev = sl_seg;
				}
			}
		}
#endif		
	}
#ifdef DEBUG_CI2
	if(start_debug())
		cerr << " Did not find any intersections " << endl;
#endif		
	return false;
}

// Douglas Peucker helper method. This method does the actual work of finding the point with the largest distance and recurses.
int cnt_bail_e = 0, cnt_bail_d = 0;
int dp(decomposition *d, float e, contour &points, int start, int end, contour *c) {
	if(start == end-1) {
		return 0;
	}
	assert(points.size() >= end);
	assert(start >= 0);
	assert(start < points.size());
	assert(end > 0);
	assert(end < points.size()+1);
	assert(end > start);
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug())
		cerr << "  DP from " << start << " (" << points[start] << ") to " << end-1 << " (" << points[end-1] << ") and e=" << e << endl;
#endif
	// Find largest/max dist:
	int max = start;
	double max_dist = -1;
	for(int i = start+1; i < end; ++i) {			
		double dist = distanceLineToPoint(points, start, end-1, i);
		if (dist > max_dist) {
			max_dist = dist;
			max = i;
		}
	}
	assert(max < end);
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug()) {
		cerr << "  longest distance to a point on " << points[start] << "->" << points[end-1] << ": " << max_dist << " on " << max << ": " << points[max] << endl;
		if(d == NULL)
			cerr << "(d NULL)" << endl;
	}
#endif

	bool ok_contains = max_dist <= e;
	if(d != NULL && ok_contains) {
		if (points[start] == points[end-1]) {
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug(points[start])) {
				cerr << " Split investigation, start " << points[start] << ", max " << points[max] << ", end-1: " << points[end-1] << endl;
			}
#endif
			ok_contains = d->contains_line(points, start, max);
			if(ok_contains) {
#ifdef DEBUG_SIMPLIFICATION
				if(start_debug(points[start])) {
					cerr << " Second part " << endl;
				}
#endif
				ok_contains = d->contains_line(points, max, end-1);
			}
			if(!ok_contains)
				cnt_bail_d++;
		} 
		else if(max < end-1) {
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug(points[start])) {
				cerr << " Full investigation, start " << points[start] << ", max " << points[max] << ", end-1: " << points[end-1] << endl;
			}
#endif
			ok_contains = d->contains_line(points, start, end-1);
			if(!ok_contains)
				cnt_bail_d++;
		}
		else {
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug(points[start])) {
				cerr << " No investigation, start " << points[start] << ", max " << points[max] << ", end-1: " << points[end-1] << endl;
			}
#endif
			bool actually_ok = d->contains_line(points, start, end-1);
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug(points[start])) {
				cerr << " Actually ok: " << actually_ok << endl;
			}
#endif
		}
	}
	else {
		cnt_bail_e++;
	}
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug(points[start]) && ok_contains) {
		cerr << " OK" << endl;
	}
#endif
	
	if (!ok_contains) {
		int rank_add = dp(d, e, points, start, max+1, c);
		return rank_add + dp(d, e, points, max, end, c);
	} else {
		contour_point end_point = points[end-1];
		if(c != NULL) {
			c->push_back(end_point);
		}
		return 1;
	}
}

void remove_pins(contour *c) {
	assert(c != NULL);
	assert(c->front() == c->back());
	if(c->size() <= 3)
		return;
#ifdef DEBUG_RP
	if(start_debug())
		cerr << " Removing pins from contour of size " << c->size() << endl;
#endif

	contour points;
	contour::iterator end = --(c->end());
	for(contour::iterator it = c->begin(); it != end; ++it) {
		points.push_back(*it);
	}
	c->clear();

	contour_point *prev = &(points.front()), *prevprev = &(points.back());
	contour::iterator it = points.begin();
	++it;
	assert(it != points.end());
#ifdef DEBUG_RP
	if(start_debug()) {
		cerr << " front: " << *prev << endl;
		cerr << " back: " << *prevprev << endl;
	}
#endif

	for(; it != points.end(); ++it) {
#ifdef DEBUG_RP
		if(start_debug())
			cerr << "  Handling: " << *it << endl;
#endif
		assert(prevprev != NULL && prev != NULL);
		if(*prevprev != *it && *prev != *it && *prevprev != *prev) {
			c->push_back(*prev);	
			prevprev = prev;
		}
#ifdef DEBUG_RP
		else if(start_debug())
			cerr << "  Removing: " << *prev << endl;
#endif
		prev = &(*it);		
	}
	assert(prevprev != NULL && prev != NULL);
	if(*prevprev != points.front() && *prev != points.front()) {
		c->push_back(points.back()); // last point.
	}
	
	c->push_back(c->front()); // Cyclic contour.
}

int fix_crossing2(contour *c, 
				  contour::iterator &it_o, contour::iterator const &o_end, 
				  contour::iterator &it_s, contour::iterator const &s_end, pt2 &cross) {
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << "Fixing crossing segment " << cross.first << "," << cross.second << endl;
#endif
//	assert(cross.lp1()->p.rank > cross.lp2()->p.rank); No can do!
	int added = 0;
	while(it_s != s_end, it_s->rank != cross.second.rank) {
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug(*it_s))
			cerr << "Adding simplified " << *it_s << endl;
#endif
		c->push_back(*it_s);
		added = it_s->rank;
		++it_s;
	}
	// Spool it_o
	while(it_o != o_end && it_o->rank < it_s->rank) {
		++it_o;
	}
	// add it_o:
	while(it_o != o_end && it_o->rank != cross.first.rank) {
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << "Adding original " << *it_o << endl;
#endif
		c->push_back(*it_o);
		added = it_o->rank;
		++it_o;
	}
	++it_s;
	return added;
}

void fix_crossing(contour *c, contour &orig, pt2 cross1, pt2 cross2) {
	int start_rank1 = cross1.second.rank;
	int start_rank2 = cross2.second.rank;
	assert(start_rank1 != start_rank2);
	if(start_rank1 > start_rank2) {
		swap(cross1,cross2);
		swap(start_rank1, start_rank2);
	}
	// Copy c:
	contour copy;
	for(contour::iterator it = c->begin(); it != c->end(); ++it) {
		copy.push_back(*it);
	}		

	// Move both to c:
	c->clear();
	contour::iterator it_o = orig.begin();
	contour::iterator it_s = copy.begin();
	// Fix cross 1:
	                 fix_crossing2(c, it_o, orig.end(), it_s, copy.end(), cross1);
	int last_added = fix_crossing2(c, it_o, orig.end(), it_s, copy.end(), cross2);

	while(it_s != copy.end()) {
		if(it_s->rank > last_added) {
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug())
				cerr << "Adding rest from simplified " << *it_s << endl;
#endif
			c->push_back(*it_s);
		}
		++it_s;
	}	
}

int cnt_rd = 0;
/*
  Douglas Peucker for a single contour.
 */
// current_contour, e_simplify, &d
void cdp(contour *c, const float e, decomposition *d,
		 stream<contour_point> &output_stream) {
	assert(d != NULL);
	assert(c != NULL);
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug(*(c->begin()))) {
		cerr << "Initiating CDP on |points|:" << c->size() << endl;
		for(contour::iterator it = c->begin(); it != c->end(); ++it) {
			cerr << " " << *it << endl;
		}
	}
#endif
	if(c->empty())
		return;

	// first point:
	contour_point p = c->front();

	contour points;
	for(contour::iterator it = c->begin(); it != c->end(); ++it) {
		points.push_back(*it);
	}
	c->clear();
	c->push_back(p);
	
	int size = 1+dp(d, e, points, 0, points.size(), c);
	
	pt2 cross1, cross2;
	if(size <= 3) {// || contains_intersections1(*c, points)) { // Load old points: // , points
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << "WARNING: Contour too small " << size << "). Reverting. " << endl;
#endif
		c->clear();
		for(contour::iterator it = points.begin(); it != points.end(); ++it) {
			c->push_back(*it);
			output_stream.write_item(*it);		
		}
	}
	else if(contains_intersections2(*c, cross1, cross2)) {
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << "WARNING: Contour contains self intersections (size " << size << "). Reverting. " << endl;
#endif
		int rd = 0;
		do {
			cnt_intersections++;
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug(cross1.first)) {
				cerr << "Fixing crossing of" << endl;
				cerr << cross1.first <<","<< cross1.second << endl;
				cerr << cross2.first <<","<< cross2.second << endl;
				for(contour::iterator it = c->begin(); it != c->end(); ++it) {
					cerr << *it << endl;
				}
				cerr << ":" << endl;
			}
#endif
			int before_size = c->size();
			fix_crossing(c, points, cross1, cross2);
			assert(c->size() > before_size);
			if(rd > cnt_rd)
				cnt_rd = rd;
			++rd;
		}
		while(contains_intersections2(*c, cross1, cross2));

		for(contour::iterator it = c->begin(); it != c->end(); ++it) {
			output_stream.write_item(*it);		
		}		
	}//*/
	else {
		for(contour::iterator it = c->begin(); it != c->end(); ++it) {
			output_stream.write_item(*it);	
		}
	}
	points.clear();
#ifdef DEBUG_SIMPLIFICATION
	assert(c != NULL);
	if(start_debug())
		cerr << "Done CDP to |c|:" << c->size() << endl;
#endif
}

void loadParentAndSiblings(int &parent, 
						   ami::queue<topo> &q_topo, ami::queue<contour_point> &q_segs, 
						   map<int,topo> &sibling_topos, map<int,contour*> &contours) {
	assert(!q_topo.is_empty());
	const topo *t;
	q_topo.dequeue(&t);
	parent = t->c;
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug())
		cerr << "Loading Q->M parent and children of parent, parent: " << parent << endl;
#endif		
	
	// Read all sibling topos from queue:
	while(!q_topo.is_empty()) { 
		q_topo.peek(&t);
		if(t->p != parent)
			break;

		q_topo.dequeue(&t);
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " Q->M: " << *t << endl;
#endif		
		sibling_topos.insert(pair<int,topo>(t->c,*t));
	}

	// load contours (parent, self, siblings):
	const contour_point *p;
	contour_point prev;
	prev.label = parent;
	int contour_to_write = parent;
	vector<contour_point> *v = new vector<contour_point>;
	while(!q_segs.is_empty()) { 
		q_segs.peek(&p);
		if(p->label != prev.label && sibling_topos.find(p->label) == sibling_topos.end())
			break;

		q_segs.dequeue(&p);
		if(p->label != contour_to_write) {
			contours.insert(pair<int,vector<contour_point>* >(contour_to_write, v));
			v = new vector<contour_point>;
		}
		v->push_back(*p);				
		contour_to_write = p->label;
		prev = *p;
	}
	contours.insert(pair<int,vector<contour_point>* >(contour_to_write, v));
}

void loadChildrenFromStream(topo* &top, contour_point* &cp, int current, stream<topo> &topology,
							map<int,topo> &children_topos,
							stream<contour_point> &input_segments, map<int,contour*> &contours) {
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug()) {
		cerr << "Loading stream->M, children of " << current << endl;
		if(top != NULL)
			cerr << " Top: " << *top << endl;
	}
#endif		
	bool halt = (top != NULL);
	while((halt || topology.read_item(&top) == NO_ERROR) && top->p == current) {
		halt = false;
		children_topos.insert(pair<int,topo>(top->c,*top));
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " Stream->M " << top->c << endl;
#endif		
	} 
#ifdef DEBUG_SIMPLIFICATIONx
	if(start_debug()) {
		cerr << "Children topos: " << endl;
		for(map<int,topo>::iterator it = children_topos.begin(); it != children_topos.end(); ++it) {
			cerr << it->first << ": " << it->second << endl;
		}
		cerr << "Children topos end " << endl;
	}
#endif		
	
	contour* c = NULL;
	contour_point last;
	// not last point && (same label || found label)
	halt = (cp != NULL);
	while((halt || input_segments.read_item(&cp) == NO_ERROR) && 
		  (cp->label == last.label || children_topos.find(cp->label) != children_topos.end())) {
		halt = false;
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << "Iterating on " << *cp << endl;
#endif		
		if(c == NULL) {
			c = new contour;
		}
		else if(last.label != cp->label) {
			contours.insert(pair<int,contour*>(last.label, c));
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug())
				cerr << " contour treamS->M: " << last.label << " , ||=" << c->size() << endl;
#endif		
			c = new contour;
		}
		c->push_back(*cp);
		assert(last.label == -1 || cp->rank != last.rank);
		last = *cp;
	}
	if(c != NULL) {
		contours.insert(pair<int,contour*>(last.label, c));
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " last contour Stream->M: " << last.label << " , ||=" << c->size() << endl;
#endif		
	}
#ifdef DEBUG_SIMPLIFICATION
	else if(start_debug()) {
		cerr << "Warning: Nothing loaded from segs!" << endl;
		if(cp != NULL)
			cerr << " cp: " << *cp << " where label is " << cp->label << endl;
	}
#endif	
}

void toQueues(contour *current_contour, topo current_topo, 
			  ami::queue<contour_point> &q_segs, 
			  ami::queue<topo> &q_topo,
			  map<int,topo> &children_topos, map<int,contour*> &contours,
			  stream<contour_point> &output) {
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug())
		cerr << "From memory To Q_segs: " << endl;
#endif		
	if(current_contour != NULL) {
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " M->Q_segs,out " << current_topo << ", ||=" << current_contour->size();
#endif		
		q_topo.enqueue(current_topo);
		int s = 0;
		for(contour::iterator it2 = current_contour->begin(); it2 != current_contour->end(); ++it2) {
			q_segs.enqueue(*it2);
			s++;
//			output.write_item(*it2); // TODO: OK?
		}
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " ||=" << s << endl;
#endif		
	}
	for(map<int,topo>::iterator it2 = children_topos.begin(); it2 != children_topos.end(); ++it2) {
		int child = it2->first;
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " M->Q_segs " << it2->second;
#endif
		q_topo.enqueue(it2->second);
		map<int,contour*>::iterator it = contours.find(child);
		assert(it != contours.end());
		current_contour = it->second;
		int s = 0;
		for(contour::iterator it3 = current_contour->begin(); it3 != current_contour->end(); ++it3) {
			q_segs.enqueue(*it3);
//			output.write_item(*it3); // TODO: OK?
			s++;
		}
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << " ||=" << s << endl;
#endif		
		// Unload own children from memory:
		contours.erase(child);
		delete current_contour;
	}
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug())
		cerr << endl;
#endif		
}

void simplification::constrained_dp(const float e_simplify,
									stream<contour_point> &input_segments,
									stream<topo> &topology,
									elev_t granularity, float e_granularity,
									stream<contour_point> &output) {
#ifdef DEBUG_SIMPLIFICATION
	if(start_debug()) {
		cerr << "Topology: " << endl;
		topo *t, prev(-10,-10,-10);
		while(topology.read_item(&t) == NO_ERROR) {
			cerr << " " << *t;
			assert(prev.p <= t->p);
			assert(prev.c < t->c);
			if(is_level_line(t->c_z, granularity, e_granularity))
				cerr << "(simplifiable)";
			cerr << endl;
			prev = *t;
		}
		topology.seek(0);
	}
#endif

	ptime t_all, t_parent;
	t_all=t_parent=microsec_clock::local_time();

	//Prepare
	output.truncate(0);

	cerr << "------------ Done initial setup ------------ " << endl;
	double millis = (microsec_clock::local_time()-t_parent).total_milliseconds();
	cout << millis << " ms." << endl;
	t_parent=microsec_clock::local_time();

	// Queues:
	ami::queue<topo> q_topo; // enqueue(T), dequeue(**T), peek(**T)
	ami::queue<contour_point> q_segs; // Move queues together!

	// Current state:
	int parent = -1;
	map<int,topo> sibling_topos; // sibling(or self) -> topo
	map<int,contour*> contours; // contour label -> points.
	contour_point *cp = NULL; // input stream state.
	topo *top = NULL; // input strem state.

	// Put -1 children from stream to Q.
	q_topo.enqueue(topo(-1,-2,-1.85230002));
	loadChildrenFromStream(top, cp, -1, topology, sibling_topos, input_segments, contours);
	toQueues(NULL, *top, q_segs, q_topo, sibling_topos, contours, output);

	int cnt_contours_simplified = 0; // OK
//	int cnt_segs_all = 0;
	int cnt_segs_simplified = 0; // OK
	int cnt_segs_simplifiable = 0; // OK
	int linear_scans = 0, bfs_scans = 0, bfs_steps = 0;

	// read t => t.p.p and siblings on queue, t.p to be simplified, read t.c.
	while(!q_topo.is_empty()) { // handle all siblings in an iteration. Break when Q is empty.
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << "------------------------------------------------------" << endl;
#endif		

		// Clear all:
		sibling_topos.clear();
		for(map<int,contour*>::iterator it = contours.begin(); it != contours.end(); ++it) {
			delete it->second;
		}
		contours.clear();
		assert(contours.empty());
		// Load parent and siblings from topo queue:
		loadParentAndSiblings(parent, q_topo, q_segs, sibling_topos, contours);
#ifdef DEBUG_SIMPLIFICATION
		if(start_debug())
			cerr << endl;
#endif		
		millis = (microsec_clock::local_time()-t_parent).total_milliseconds();
		cout << endl << "Parent " << parent << ": pre: " << millis << "ms. ";
			
		ptime t_child=microsec_clock::local_time();						

		// handle all sibling:
		for(map<int,topo>::iterator it = sibling_topos.begin(); it != sibling_topos.end(); ++it) {
			int current = it->second.c;
			elev_t current_z = it->second.c_z;
			bool simplifyable = is_level_line(current_z, granularity, e_granularity) && parent != -1;
#ifdef DEBUG_SIMPLIFICATION
			if(start_debug()) {
				cerr << "Handling sibling " << it->first << " at topo:" << it->second << endl;
				if(simplifyable)
					cerr << "(Simplifiable)" << endl;
			}
#endif		
				
			// Load children contours from stream (to queue after simplified self):
			map<int,topo> children_topos;
			loadChildrenFromStream(top, cp, current, topology, children_topos, input_segments, contours);

			map<int,contour*>::iterator it3 = contours.find(current);
			assert(it3 != contours.end());
			contour *current_contour = it3->second;

			millis = (microsec_clock::local_time()-t_child).total_milliseconds();
			cout << endl << " Child " << current << " load: " << millis << "ms. ";
			t_child=microsec_clock::local_time();					


			if(simplifyable) {
				// Actually simplify t->p.
				decomposition d(parent, current, contours);
				cnt_segs_simplifiable += current_contour->size();
				millis = (microsec_clock::local_time()-t_child).total_milliseconds();
				cout << "decomp: " << millis << "ms. ";
				t_child=microsec_clock::local_time();					

				cdp(current_contour, e_simplify, &d, output); // current contour is changed.
#ifdef DEBUG_SIMPLIFICATION
				if(start_debug())
					cerr << "Simplified " << current << endl;
#endif		
				linear_scans += d.linear_scans;
				bfs_scans += d.bfs_scans;
				bfs_steps += d.bfs_steps;
				cnt_segs_simplified += current_contour->size();
				millis = (microsec_clock::local_time()-t_child).total_milliseconds();
				cout << "simp: " << millis << "ms. ";
				t_child=microsec_clock::local_time();					
				cnt_contours_simplified++;
			}
			else {
				for(contour::iterator it6 = current_contour->begin(); it6 != current_contour->end(); ++it6) {
					output.write_item(*it6);		
				}		
			}

			// INSERT (simplified) self and children INTO Queues (topo queue already done):
			toQueues(current_contour, it->second, q_segs, q_topo, children_topos, contours, output);
			millis = (microsec_clock::local_time()-t_child).total_milliseconds();
			cout << "write: " << millis << "ms. ";
			t_child=microsec_clock::local_time();					
		}
		millis = (microsec_clock::local_time()-t_parent).total_milliseconds();
		cout << endl << "Total: " << millis << "ms. ";
		t_parent=microsec_clock::local_time();					
	}

	// TODO: Update paper with BFS and selv in queue?			
	input_segments.seek(0);
	topology.seek(0);
	output.seek(0);
	cout << endl;
	cout << " Time usage for cdp in total: " << (microsec_clock::local_time()-t_all) << " ms." << endl;
	cout << "|topology| (stream len): " << topology.stream_len() << endl;
	cout << "#|All contours| (stream len): " << input_segments.stream_len()-topology.stream_len() << endl;
	cout << "#|simplifiable segments|: " << cnt_segs_simplifiable << endl;
	cout << "#|simplified segments|: " << cnt_segs_simplified << endl;
	cout << "#|output segments| (stream len): " << output.stream_len()-topology.stream_len() << endl;
	cout << "#Intersections: " << cnt_intersections << endl;
    cout << "#Max recursion for fixing crossings: " << cnt_rd << endl;
    cout << "#Extra linear scans: " << linear_scans << endl;
    cout << "#BFS scans and steps: " << bfs_scans << ", " << bfs_steps << endl;
    cout << "#bail for epsilon: " << cnt_bail_e << endl;
    cout << "#bail for decomposition: " << cnt_bail_d << endl;
}
