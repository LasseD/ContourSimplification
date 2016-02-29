// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_DECOMPOSITION_H__
#define __TEST_CONTOUR_SIMPLIFICATION_DECOMPOSITION_H__
#include <terrastream/common/common.h>
#include "io_contours/contour_types.h"
#include <tpie/array.h>
#include <set>
#include <vector>
#include <iomanip>

using namespace terrastream;
using namespace tpie;
using namespace std;

//#define DEBUG_SLC
//#define DEBUG_SLC_OLD

const float SQRT_EPSILON = 1e-7f; // sqrt EPSILON
const float EPSILON = SQRT_EPSILON*SQRT_EPSILON; // Only to compare when points are too close!

typedef vector<contour_point> contour;

namespace simplification {

bool start_debug(contour_point p);
bool start_debug();
bool is_level_line(elev_t candidate, elev_t granularity, float epsilon);

struct link_point {
	contour_point p;
	link_point *prev, *next, *up, *down;

	link_point(const link_point &lp);
	link_point(contour_point &p);

	bool operator<(const link_point l);
	bool operator==(const link_point l);
	
	bool is_starting();
	bool is_ending();
	bool is_transitioning();
	bool inside_up(int parent_contour);
	bool inside_up();
	bool inside_down(int parent_contour);
	bool inside_down();
	static void connect(link_point *prev, link_point *next);
	static void connect_vertical(link_point *u, link_point *d);

	inline friend std::ostream& operator << (std::ostream &ostr, const link_point &t) {
		ostr << t.p;

		link_point* points[4] = {t.prev,t.next,t.up,t.down};
		string desc[4] = {"prev", "next", "up", "down"};
		for(int i = 0; i < 4; i++) {
			link_point *point = points[i];
			ostr << desc[i] << ": ";
			if(point != NULL) {
				ostr << point->p << " ";
			}
			else {
				ostr << "- ";
			}
		}
		return ostr;
	}
};

struct lp_seg {
	lp_seg(link_point *lp);
	lp_seg(const lp_seg &s);
	lp_seg();
	lp_seg(contour_point &a);
	lp_seg(contour_point &a, contour_point &b);
	~lp_seg();
	contour_point p1() const;
	contour_point p2() const;
	link_point* lp1() const;
	link_point* lp2() const;
	link_point* lp_at(link_point *guide) const;

	bool intersects(contour_point &p1, contour_point &p2);
	lp_seg seg_above();
	lp_seg seg_below();
	bool contains_p_above(contour_point p,contour_point p_guide);
	lp_seg split_from_below(link_point *lp, int parent);
	lp_seg split_from_above(link_point *lp, int parent);
	bool inside_up(int parent_contour) const;
	xycoord_t eval(xycoord_t x) const; // Deprecate anything but use for creation! (ie. privatice)
	bool is_vertical() const;

	const bool operator==(const lp_seg &l) const { // for Sweep line
		return (lp1() == l.lp1() && lp2() == l.lp2()) ||
			(lp1() == l.lp2() && lp2() == l.lp1());
	}

/*	const bool operator<(const lp_seg &l) const { // for Sweep line
//		cerr << *this << " vs " << l << endl;
		if(*this == l)
			return false;		
		xycoord_t x = max(p1().x, l.p1().x);
		xycoord_t y = eval(x);
		xycoord_t ly = l.eval(x);
		if(y != ly) {
			return y < ly;
		}

		// start/end special cases:
		bool at_min = x == min(p1().x, p2().x);
		bool at_l_min = x == min(l.p1().x, l.p2().x);
		if(x == max(p1().x, p2().x) && at_l_min) {
			return true;
		}
		if(x == max(l.p1().x, l.p2().x) && at_min) {
			return false;
		}

		x--; // new comparing place. // TODO: DON'T COMPUTE NEW POINTS!
		if(at_min || at_l_min) {
			x +=2;
		}

		y = eval(x);
		ly = l.eval(x);
		if(y != ly)
			return y < ly;		
		return !inside_up(); // Warning: no inside_up cmp.
		}*/

	inline friend std::ostream& operator << (std::ostream &ostr, const lp_seg &t) {
		if(t.lp1() == NULL)
			return ostr << "Not set";
		ostr << setprecision(11) << t.p1() << "-" << setprecision(11) << t.p2();
//		ostr << " pointer: " << t.lp1();
		return ostr;
	}
private:
	bool intersects_line(contour_point &lp1, contour_point &lp2);
	bool inside_up() const; // of contour, ignoring label.
	link_point *split(xycoord_t x);
	link_point *lp;
};

int int_sign(float f);
int int_sign(float f, float g);

struct sl_cmp {
	int parent;
	xycoord_t *x;
	bool assertion_raised;
	contour_point asa1,asa2,asb1,asb2;

	sl_cmp(int p, xycoord_t *_x) : parent(p), x(_x), assertion_raised(false) {}

	bool has_assertion_raised(contour_point &a1, contour_point &a2, contour_point &b1, contour_point &b2) {
		if(!assertion_raised)
			return false;
		a1 = asa1;
		a2 = asa2;
		b1 = asb1;
		b2 = asb2;
		return true;
	}

	int det_cmp(contour_point const &p1,contour_point const &p2, 
				contour_point const &p3,contour_point const &p4) {
#ifdef DEBUG_SLC
		if(start_debug(p1)) {
			cerr << " det cmp " << p1 << "-" << p2 << " vs " << p3 << "-" << p4 << endl;
		}
#endif
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
		int s123 = int_sign(det123);
		int s124 = int_sign(det124);
		int s341 = int_sign(det341);
		int s342 = int_sign(det342);

		if(x1 == x2) { // First line vertical!
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " Vertical 12" << endl;
			}
#endif
			if(x3 == x4) { // Other line vertical!
#ifdef DEBUG_SLC
				if(start_debug(p1)) {
					cerr << " Vertical both!" << endl;
				}
#endif
				if(min(y1,y2) < max(y3,y4))
					return -1;
				if(min(y3,y4) < max(y1,y2))
					return 1;
				return 0;			
			}
			if(s341 == 0)
				return s342; // overlapping end point
			if(s342 == 0)
				return s341; // overlapping end point
			assertion_raised |= (s341 != s342);
			return s341;			
		}
		if(x3 == x4) { // Other line vertical!
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " Vertical 34" << endl;
			}
#endif
			if(s123 == 0)
				return -s124; // overlapping end point
			if(s124 == 0)
				return -s123; // overlapping end point
			assertion_raised |= (s123 != s124);
			return -s123;			
		}

#ifdef DEBUG_SLC
		if(start_debug(p1)) {
			cerr << " det123 " << det123 << ", det124 " << det124 << ", det341 " << det341 << ", det342 " << det342 << endl;
		}
#endif

		if(!(s123 == 0 && s124 == 0)) { // can be 0 because of small s12 length
			if(s123 == 0)
				return -s124; // overlapping end point
			if(s124 == 0)
				return -s123; // overlapping end point
			if(s123 == s124) { // If they don't agree, ask s34.
				return -s123;
			}
		}
		if(s342 == 0)
			return s341; // overlapping end point
		if(s341 == 0)
			return s342; // overlapping end point
		assertion_raised |= (s341 != s342);
		return s341;
	}

	int compare_y(contour_point p1,contour_point p2, 
				  contour_point p3,contour_point p4) {
		if(p1.x > p2.x)
			swap(p1,p2);
		if(p3.x > p4.x)
			swap(p3,p4);
#ifdef DEBUG_SLC
		if(start_debug(p1)) {
			cerr << " Comparing " << p1 << "-" << p2 << " vs " << p3 << "-" << p4 << endl;
		}
#endif
		// Overlapping points:
		assert(p1 != p2);
		assert(p3 != p4);
		if((p1 == p3 && p2 == p4) || (p1 == p4 && p2 == p3)) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " Same line" << endl;
			}
#endif
			return 0;
		}
		if(p2 == p3 && p1.x < p2.x && p3.x < p4.x) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " Meet 1: -1" << endl;
			}
#endif
			return -1;
		}
		if(p4 == p1 && p1.x < p2.x && p3.x < p4.x) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " Meet 2: 1" << endl;
			}
#endif
			return 1;
		}

		// Bounding box:
		xycoord_t miny12 = min(p1.y,p2.y);
		xycoord_t maxy12 = max(p1.y,p2.y);
		xycoord_t miny34 = min(p3.y,p4.y);
		xycoord_t maxy34 = max(p3.y,p4.y);
		if(maxy12 <= miny34) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " minmaxy1: -1" << endl;
			}
#endif
			return -1;
		}
		if(maxy34 <= miny12) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " minmaxy2: 1" << endl;
			}
#endif
			return 1;
		}

		xycoord_t minx12 = p1.x;
		xycoord_t maxx12 = p2.x;
		xycoord_t minx34 = p3.x;
		xycoord_t maxx34 = p4.x;
		assert(minx12 <= *x);
		assert(maxx12 >= *x);
		assert(minx34 <= *x);
		assert(maxx34 >= *x);
		assert(maxx12 >= minx34);
		assert(maxx34 >= minx12);

		if(minx12 < minx34) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " det1" << endl;
			}
#endif
			return det_cmp(p1,p2,p3,p4);
		}
 		if(minx12 > minx34) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " -det2" << endl;
			}
#endif
			return -det_cmp(p3,p4,p1,p2);
		}
		if(maxx12 < maxx34) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " det3" << endl;
			}
#endif
			return det_cmp(p1,p2,p3,p4);
		}
 		if(maxx12 > maxx34) {
#ifdef DEBUG_SLC
			if(start_debug(p1)) {
				cerr << " -det4" << endl;
			}
#endif
			return -det_cmp(p3,p4,p1,p2);
		}

#ifdef DEBUG_SLC
		if(start_debug(p1)) {
			cerr << " y split" << endl;
		}
#endif

		if(miny12 < miny34 || maxy12 < maxy34)
			return -1;
		if(miny12 > miny34 || maxy12 > maxy34)
			return 1;
		return 0;
	}

	int compare_y(lp_seg const &a, lp_seg const &b) {
		return compare_y(a.p1(),a.p2(), b.p1(),b.p2());
	}

	bool old_cmp(lp_seg const &a, lp_seg const &b) const {
#ifdef DEBUG_SLC_OLD
		if(start_debug(a.pe())) {
			cerr << "Old compare x:" << *x << ", " << a << " vs " << b << endl;
		}
#endif
		if(a == b)
			return false;		

		xycoord_t y = a.eval(*x);
		xycoord_t ly = b.eval(*x);
		if(y != ly) { // && abs(y-ly) > EPSILON) { // Safe distance.
#ifdef DEBUG_SLC_OLD
			if(start_debug(a.pe())) {
				cerr << "y:" << y << " ly:" << ly << endl;
			}
#endif
			return y < ly;
		}
		
		// start/end special cases:
		xycoord_t min_a = min(a.p1().x, a.p2().x);
		xycoord_t min_b = min(b.p1().x, b.p2().x);
		if(max(a.p1().x, a.p2().x) == min_b) {
			assert(*x == min_b);
#ifdef DEBUG_SLC_OLD
			if(start_debug(a.pe())) {
				cerr << "Special case 1: true" << endl;
			}
#endif
			return true;
		}
		if(min_a == max(b.p1().x, b.p2().x)) {
			assert(*x == min_a);
#ifdef DEBUG_SLC_OLD
			if(start_debug(a.pe())) {
				cerr << "Special case 2: false" << endl;
			}
#endif
			return false;
		}
		
		xycoord_t xx = *x-1; // new comparing place. // TODO: DON'T COMPUTE NEW POINTS! FAILS ON LARGE INPUT!
		if(*x <= min_a || *x <= min_b) {
			xx = *x+1;
		}
		
		y = a.eval(xx);
		ly = b.eval(xx);
		if(y != ly)
			return y < ly;		
/*		if(a.inside_up(parent) != b.inside_up(parent)) {
			cerr << " a " << a << " inside up: " << a.inside_up(parent) << " vs " << b << endl;
			return !a.inside_up(parent);
			}*/
		return false;
	}

	// a < b:
	bool operator() (lp_seg const &a, lp_seg const &b) {
		bool before = assertion_raised;
		int cy = compare_y(a,b);
		bool res = cy < 0;
#ifdef DEBUG_SLC
		if(start_debug()) {
			bool old_res = old_cmp(a,b);
			cerr << " cy: " << cy << ", res: " << res << ", old res: " << old_res << endl;
			if(assertion_raised)
				cerr << "ASSERTION RAISED!" << endl;
		}
#endif
//		assert(res == old_res);
		if(!before && assertion_raised) {
			asa1 = a.p1();
			asa2 = a.p2();
			asb1 = b.p1();
			asb2 = b.p2();
		}

		return res;
	}

	bool operator() (lp_seg const &a, link_point *b, bool left) {
		// Find real b to make seg from::
		if(left) {
			if(b->prev->p.x > b->p.x || (b->prev->p.x == b->p.x && b->prev->p.y > b->p.y))
				b = b->next;
		}
		else {
			if(b->prev->p.x < b->p.x || (b->prev->p.x == b->p.x && b->prev->p.y < b->p.y))
				b = b->next;
		}

		return (*this)(a, lp_seg(b));
	}
};

class decomposition {
public:
	decomposition(int parent, int skip, map<int,vector<contour_point>* > &contours);
	~decomposition();
    bool contains_line(vector<contour_point> &v, int p1, int p2);
	static link_point* addContourToPQ(set<link_point*,bool(*)(link_point*,link_point*)> &pq, contour &c, bool d);
	int linear_scans, bfs_scans, bfs_steps;
private:
	int parent_contour;
	vector<link_point*> points;
	link_point *guide;

	bool walk_contour(contour_point &p1, contour_point &p2, link_point* &start, vector<link_point*> &out);
	link_point* find(contour_point p,contour_point p_guide); // finds self and prev of seg below p.
	link_point* bfs_find(contour_point p,contour_point p_guide,link_point* guide); // finds self and prev of seg below p.
	void print();
};

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_DECOMPOSITION_H__*/
