// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
//
// Copyright 2010 SCALGO development

#include "io_contours/contour_types.h"
#include "decomposition.h"
#include <set>
#include <tpie/array.h>
#include "util.h"
#include <iomanip>

//#define DEBUG_DECOMPOSITION_CONSTRUCT
//#define DEBUG_DECOMPOSITION_WALK

using namespace simplification;

typedef vector<contour_point> contour;

namespace simplification {

const int DEBUG_START = 1115610; // 6619
bool START_DEBUG = false;
bool start_debug(contour_point p) {
	if(START_DEBUG || p.label >= DEBUG_START) {
		START_DEBUG = true;
		return true;
	}
	return false;
}

int int_sign(float f, float precision) {
	if(abs(f) < precision) {
		return 0;
	}
	if(f < 0)
		return -1;
	if(f > 0)
		return 1;
	return 0;
}
int int_sign(float f) {
	return int_sign(f,EPSILON);
}

bool start_debug() {
	if(DEBUG_START <= 0)
		START_DEBUG = true;
	return START_DEBUG;
}

bool is_level_line(elev_t candidate, elev_t granularity, float epsilon) {
	bool res = abs(round(candidate/granularity)*granularity - candidate) < epsilon/30;
	return res;
}

link_point::link_point(const link_point &lp) : p(lp.p), prev(lp.prev), next(lp.next), up(lp.up), down(lp.down) {
}

link_point::link_point(contour_point &cp) : p(cp), prev(NULL), next(NULL), up(NULL), down(NULL) {
}

inline xycoord_t get_y(contour_point &p1, contour_point &p2, xycoord_t x) {
	// TODO: Outphase this method!
	if(p1.x == p2.x && p1.x == x) { // vertical
		return (p1.y+p2.y)/2;
	}
	if(p1.x == x)
		return p1.y;
	if(p2.x == x)
		return p2.y;
	return p1.y + (p2.y - p1.y)*(x - p1.x)/(p2.x - p1.x);
}

// lex x, y, label, rank for PQ
bool link_point::operator<(const link_point l) {
	if(p.x != l.p.x)
		return l.p.x > p.x;
	if(p.y != l.p.y)			
		return l.p.y > p.y;
	if(p.label != l.p.label)			
		return l.p.label > p.label;
	return p.rank < l.p.rank;
}
bool link_point::operator==(const link_point l) { // for up/down shooting
	return l.p.x == p.x && l.p.y == p.y;
}

bool link_point::is_starting() {	
	xycoord_t x = min(next->p.x, prev->p.x);
	if(p.x == x) { // vertical.
		if(next->p.x == p.x && prev->p.x == x)
			return false; // transitioning.
		if(next->p.x == p.x)
			return next->p.y > p.y;
		assert(prev->p.x == p.x);
		return prev->p.y > p.y;
	}
	return p.x < x;
}

bool link_point::is_ending() {
	xycoord_t x = max(next->p.x, prev->p.x);
	if(p.x == x) { // vertical.
		if(next->p.x == p.x && prev->p.x == x)
			return false; // transitioning.
		if(next->p.x == p.x)
			return next->p.y < p.y;
		assert(prev->p.x == p.x);
		return prev->p.y < p.y;
	}
	return p.x > x;
}

bool link_point::is_transitioning() {
	xycoord_t ma = max(next->p.x, prev->p.x);
	xycoord_t mi = min(next->p.x, prev->p.x);
	/*cerr << "p: " << p << " ma: " << ma << " mi: " << mi << endl;
	cerr << "next: " << *next << endl;
	cerr << "prev: " << *prev << endl;//*/

	// Real transitioning:
	if((ma > p.x && mi < p.x) || (ma == p.x && mi == p.x)) {
		return true;
	}
	
	if(p.x == ma) { // vertical/left:
		if(next->p.x == p.x)
			return next->p.y > p.y;
		assert(prev->p.x == p.x);
		return prev->p.y > p.y;
	}
	if(p.x != mi) {
		return false;
	}
    if(next->p.x == p.x) {
	    return next->p.y < p.y;
	}
	assert(prev->p.x == p.x);
	return prev->p.y < p.y;
}

lp_seg::lp_seg(link_point *p) {
	assert(p != NULL);
	lp = p;
}
lp_seg::lp_seg(const lp_seg &s) {
	lp = s.lp;
	assert(lp != NULL);
}
lp_seg::lp_seg() {
	lp = NULL; // special constructor.
}

link_point* lp_seg::lp1() const {
	return lp;
}
link_point* lp_seg::lp2() const {
	return lp->prev;
}
contour_point lp_seg::p1() const {
	return lp1()->p;
}
contour_point lp_seg::p2() const {
	return lp2()->p;
}
link_point* lp_seg::lp_at(link_point *a) const {
	return p1() == a->p ? lp1() : lp2();
}

xycoord_t lp_seg::eval(xycoord_t x) const {
	contour_point a = p1();
	contour_point b = p2();
	return get_y(a, b, x);
}

bool lp_seg::inside_up(int parent_contour) const {
	assert(lp != NULL);
	bool swap = lp->p.label == parent_contour;	
	return swap == inside_up();
}

bool lp_seg::inside_up() const {
	if(p1().x != p2().x)
		return p1().x < p2().x;
	return p1().y > p2().y;
}

bool link_point::inside_up(int parent) {
	if(prev->p.x == p.x && next->p.x == p.x) { // vertical.
		return false;
	}
	if(prev->p.x == p.x && prev->p.y > p.y) { // vertical up
		return false;
	}
	if(next->p.x == p.x && next->p.y > p.y) { // vertical up
		return false;
	}
	bool swap = p.label == parent;	
	return swap != inside_up();
}
bool link_point::inside_down(int parent) {
	if(prev->p.x == p.x && next->p.x == p.x) { // vertical down
		return false;
	}
	if(prev->p.x == p.x && prev->p.y < p.y) { // vertical down
		return false;
	}
	if(next->p.x == p.x && next->p.y < p.y) { // vertical down
		return false;
	}

	bool swap = p.label == parent;	
	return swap != inside_down();
}

bool link_point::inside_down() {
	// Vertical:
	if(prev->p.x == p.x) {
		return next->p.x < p.x;
	}
	if(next->p.x == p.x) {
		return prev->p.x > p.x;
	}
	// left/right local extrema:
	if(prev->p.x < p.x && next->p.x < p.x) {
		xycoord_t x = p.x-1;
		return get_y(prev->p, p, x) > get_y(next->p, p, x);
	}
	if(prev->p.x > p.x && next->p.x > p.x) {
		xycoord_t x = p.x+1;
		return get_y(prev->p, p, x) < get_y(next->p, p, x);
	}
	return prev->p.x > p.x;
}

bool link_point::inside_up() {
	// Vertical:
	if(prev->p.x == p.x) {
		return next->p.x > p.x;
	}
	if(next->p.x == p.x) {
		return prev->p.x < p.x;
	}
	// left/right local extrema:
	if(prev->p.x < p.x && next->p.x < p.x) {
		xycoord_t x = p.x-1;
		return get_y(prev->p, p, x) > get_y(next->p, p, x);
	}
	if(prev->p.x > p.x && next->p.x > p.x) {
		xycoord_t x = p.x+1;
		return get_y(prev->p, p, x) < get_y(next->p, p, x);
	}
	return prev->p.x < p.x;
}

link_point* lp_seg::split(xycoord_t x) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug(p1())) {
		cerr << "  Splitting " << *this << " at " << x;
	}
#endif
	assert(x < max(p1().x,p2().x));
	assert(x > min(p1().x,p2().x));
	contour_point c(x, eval(x), -1, lp->p.label);
	link_point *np = new link_point(c);
	link_point::connect(lp->prev, np);
	link_point::connect(np, lp);
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug(p1())) {
		cerr << " into " << *this << " and new " << *np << endl;
	}
#endif
	return np;
}

lp_seg lp_seg::split_from_below(link_point *llp, int parent) {
	if((llp->next->p.x == llp->p.x && llp->next->p.y > llp->p.y) ||
	   (llp->prev->p.x == llp->p.x && llp->prev->p.y > llp->p.y)) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug(p1())) {
		cerr << "  Vertical self bail" << endl;
	}
#endif
		assert(false);
	}
	if(inside_up(parent)) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug(p1())) {
		cerr << "  Inside up bail " << endl;
	}
#endif
		return *this;
	}
	
	link_point *own_lp;
	bool no_split = true;
	if(llp->p.x == p1().x) {
		own_lp = lp1();
	}
	else if(llp->p.x == p2().x) {
		own_lp = lp2();
	}
	else {
		own_lp = split(llp->p.x);
		no_split = false;
	}
	while(no_split && own_lp->down != NULL && own_lp->down->p.y >= llp->p.y) {
		//assert(false); // TODO: THIS SHOULD NOT HAPPEN! -- see sweep line invariant.
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(p1())) {
			cerr << "  Walking down from " << *own_lp << " to " << *(own_lp->down) << endl;
		}
#endif
		own_lp = own_lp->down;
		if(*own_lp == *llp)
			return *this;//*/
	}
	link_point::connect_vertical(own_lp, llp);
	if(!no_split && own_lp->p.x > lp->p.x) {
		return lp_seg(own_lp);
	}
	return *this;
}

lp_seg::lp_seg(contour_point &a) { // TODO: Make factory contructor!
	lp = new link_point(a);
	lp->prev = (lp->next = lp);
}
lp_seg::lp_seg(contour_point &a,contour_point &b) { // TODO: Make factory contructor!
	lp = new link_point(a);
	link_point *lpp = new link_point(b);
	lp->prev = lpp;
	lp->next = lp;
}
lp_seg::~lp_seg() {
	if(lp != NULL && lp->next == lp) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(p1())) {
			cerr << "  Deleting special point lp_seg" << *lp << endl;
		}
#endif
		if(lp->prev != lp)
			delete lp->prev;
		delete lp;
	}
	lp = NULL;
}

lp_seg lp_seg::split_from_above(link_point *llp, int parent) {
	if((llp->next->p.x == llp->p.x && llp->next->p.y < llp->p.y) ||
	   (llp->prev->p.x == llp->p.x && llp->prev->p.y < llp->p.y)) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(p1())) {
			cerr << "  Vertical self bail" << endl;
		}
#endif
		assert(false);
	}
	if(!inside_up(parent)) {
		return *this;
	}

	link_point *own_lp;
	bool no_split = true;
	if(llp->p.x == p1().x) {
		own_lp = lp1();
	}
	else if(llp->p.x == p2().x) {
		own_lp = lp2();
	}
	else {
		own_lp = split(llp->p.x);
		no_split = false;
	}
	while(no_split && own_lp->up != NULL && own_lp->up->p.y <= llp->p.y) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(p1())) {
			cerr << "  Walking up closer from " << own_lp->p.rank << " to " << (own_lp->up)->p.rank << endl;
		}
#endif
//		assert(own_lp->up->up != own_lp);
		own_lp = own_lp->up;		
		if(*own_lp == *llp)
			return *this;
	}
	link_point::connect_vertical(llp, own_lp);
	if(!no_split && own_lp->p.x > lp->p.x) {
		return lp_seg(own_lp);
	}
	return *this;
}

void link_point::connect(link_point *prev, link_point *next) {
	assert(prev != NULL);
	assert(next != NULL);
	assert(prev != next);
	prev->next = next;
	next->prev = prev;
	if(prev->p.x == next->p.x && prev->p.y != next->p.y) { // Vertical handling here!
		if(prev->p.y < next->p.y) {
			connect_vertical(next, prev); // Vertical no longer need to be handled in PQ.
		}
		else {
			connect_vertical(prev, next); // Vertical no longer need to be handled in PQ.
		}
	}
}

void link_point::connect_vertical(link_point *up, link_point *down) {
	assert(up != NULL);
	assert(down != NULL);
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug(up->p)) {
		cerr << "  Connecting vertically, up/down:" << endl;
		cerr << "   " << *up << endl;
		cerr << "   " << *down << endl;
	}
#endif
	assert(up->p.x == down->p.x);
	assert(up->p.y >= down->p.y);

	while(up->down != NULL && (up->down->p.y > down->p.y || (up->down->p.y == down->p.y && up->down!=down))) {
		up = up->down;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(up->p)) {
			assert(up->down == NULL || up->down->down != up);
			cerr << "  Stepping down on up:" << *up << endl;
		}
#endif
	}
	while(down->up != NULL && (down->up->p.y < up->p.y || (down->up->p.y == up->p.y && up!=down->up))) {
		down = down->up;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug()) {
			assert(down->up == NULL || down->up->up != down);
			cerr << "  Stepping up on down:" << *down << endl;
		}
#endif
	}		

	if(down->up != NULL && down->up != up) { // merge in:
		link_point *over_up = up;
		while(over_up->up != NULL && over_up->up != down->up) {
			over_up = over_up->up;
			if(over_up->p.y > down->up->p.y)
				return; // Already set.
		}
		over_up->up = down->up;
		down->up->down = over_up;
	}
	if(up->down != NULL && up->down != down) { // merge in:
		link_point *over_down = down;
		while(over_down->down != NULL && over_down->down != up->down) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug()) {
				cerr << "  Downing:" << *over_down << endl;
			}
#endif
			over_down = over_down->down;
			if(over_down->p.y < up->down->p.y)
				return; // already set.
		}
		over_down->down = up->down;
		up->down->up = over_down;
	}

	up->down = down;
	down->up = up;
}

bool lp_seg::is_vertical() const {
	return p1().x == p2().x;
}

bool ptr_cmp(link_point* a, link_point* b) {
	return *a < *b;
}

typedef set<lp_seg,sl_cmp> lpsegpset;

lpsegpset::iterator scan_for(lp_seg &seg, lpsegpset sweep_line) {
	lpsegpset::iterator it; 
	for(it = sweep_line.begin(); it != sweep_line.end(); ++it) {
		if(*it == seg)
			return it;
	}
	return it;
}

inline void update_sweep_line(lpsegpset &sweep_line, lp_seg &seg, bool is_ending) {
	if(seg.is_vertical()) {
		return;
	}

	if(!is_ending) {
		pair<lpsegpset::iterator,bool> pair = sweep_line.insert(seg);
		assert(pair.second);
	}
	else {
		bool erased = sweep_line.erase(seg);
		assert(erased);
	}
}

// returns true if segment is left.
inline bool lr_insert(link_point *&ll, link_point *&lu, link_point *&rl, link_point *&ru, 
					  link_point *&vl, link_point *&vu, 
					  lp_seg &seg, link_point *point) {
	assert(point != NULL);
	link_point *copy = seg.lp1();
	link_point *other = seg.lp2();
	if(seg.p2() == point->p) {
		swap(other,copy);
	}
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug(seg.p1())) {
		cerr << "LR insert " << seg << " at " << *point << endl;
		cerr << "copy: " << *copy << ", other: " << *other << endl;
	}
#endif
	assert(copy->p == point->p);

	if(seg.is_vertical()) {
		if(other->p.y < point->p.y) { // left down:
			vl = copy;
			return true; // unused return value
		}
		else {// right up:
			vu = copy;
			return false; // unused return value
		}
	}
	else {
		xycoord_t x = point->p.x;
		sl_cmp cmp(-1, &x);
		if(other->p.x < x) { // left
			if(ll == NULL || cmp(seg, ll, true)) {
				ll = copy;
			}
			if(lu == NULL || !cmp(seg, lu, true)) {
				lu = copy;
			}
			return true;
		}
		else {// right
			if(rl == NULL || cmp(seg, rl, false)) {
				rl = copy;
			}
			if(ru == NULL || !cmp(seg, ru, false)) {
				ru = copy;
			}
			return false;
		}
	}
}

link_point* decomposition::addContourToPQ(set<link_point*,bool(*)(link_point*,link_point*)> &pq, contour &c, bool do_contractions) {
	bool first = true;
	link_point *first_point = NULL, *prev = NULL;
	for(contour::iterator it2 = c.begin(); it2 != c.end(); ++it2) {
		assert(first_point == NULL || first_point->prev == NULL);
		assert(prev == NULL || it2->rank == -1 || it2->rank != prev->p.rank);
		contour_point p = *it2;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(p)) {
			cerr << " Adding " << p << endl;
		}
#endif
		assert(prev == NULL || prev->p != p);
		link_point *lp;
		if(first) {
			lp = new link_point(p);
			first_point = lp;
			first = false;
		}
		else if(p == first_point->p) {
			link_point::connect(prev, first_point);
//				break; // ensures we stop when first self-loop to 0 occurs. 
		}
		else {
			lp = new link_point(p);
			link_point::connect(prev, lp);
		}
		prev = lp;
	}
#ifdef DEBUG_DECOMPOSITION_CONSTRUCTx
	if(start_debug(p)) {
		cerr << " Added all to lp-list" << endl;
	}
#endif				
	assert(first_point != NULL);
	assert(first_point->prev != NULL); // never got back to start?

	// Remove Double horizontal and overlapping segments: // TODO: Not necessary with more robust comparator!
	if(do_contractions) {
		prev = first_point;
		do {
			if(prev->prev->p.x == prev->next->p.x && prev->prev->p.y == prev->next->p.y) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
				if(start_debug()) {
					cerr << "Oh noes: " << *prev << endl;
					assert(false);
				}
#endif
			}
			if((prev->prev->p.x == prev->next->p.x && prev->p.x == prev->next->p.x) ||
			   (prev->prev->p.y == prev->next->p.y && prev->p.y == prev->next->p.y)) {
				// remove prev:
				prev->prev->next = prev->next;
				prev->next->prev = prev->prev;
				if(prev->up != NULL) {
					assert(prev->up->down == prev);
					prev->up->down = prev->down;
				}
				if(prev->down != NULL) {
					assert(prev->down->up == prev);
					prev->down->up = prev->up;
				}
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
				if(start_debug(prev->p)) {
					cerr << " Removed " << *prev << endl;
				}
#endif				
			}
			assert(prev->next != prev);
			prev = prev->next;
		}
		while(prev != first_point);
#ifdef DEBUG_DECOMPOSITION_CONSTRUCTx
		if(start_debug(prev->p)) {
			cerr << " Removed all dbl. horizontal and overlapping." << endl;
		}
#endif			
	}	
	
	// Insert to PQ
	prev = first_point;
	do {
		assert(prev->prev != NULL);
		assert(prev->next != NULL);
	    bool inserted = pq.insert(prev).second;
		assert(inserted);
		prev = prev->next;
	}
	while(prev != first_point);
	
	return first_point;
}

decomposition::decomposition(int parent, int skip, map<int,contour* > &contours) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug())
		cerr << "Constructing decomposition for parent " << parent << " without " << skip << endl;
#endif
	guide = NULL;
	parent_contour = parent;
	linear_scans = -1;
	bfs_scans = 0;
	bfs_steps = 0;
	// add all points:
	set<link_point*,bool(*)(link_point*,link_point*)> pq(ptr_cmp);

	for(map<int,vector<contour_point>* >::iterator it = contours.begin(); it != contours.end(); ++it) {
		if(it->first == skip || it->second->size() < 4) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug())
				cerr << "Ignoring contour " << it->first << " of size " << it->second->size() << endl;
#endif
			continue;
		}
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug())
			cerr << "Adding contour " << it->first << " of size " << it->second->size() << endl;
#endif
		// Add all points to pq and 'points':
		contour *cont = it->second;
		link_point *first_point = addContourToPQ(pq, *cont, true);
		assert(first_point != NULL);
		points.push_back(first_point);
	}
	// Do sweep for splitting and up/down setting.
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug()) {
		cerr << "Starting sweep on |PG| = " << pq.size()<< endl;
		print_set(pq);
		cerr << "---------------------------" << endl;
	}
#endif
	
	xycoord_t x;
	sl_cmp sweep_line_cmp(parent, &x);
	lpsegpset sweep_line(sweep_line_cmp);

	for(set<link_point*,bool(*)(link_point*,link_point*)>::iterator it_pq = pq.begin(); it_pq != pq.end();) {
		link_point *p = *it_pq;
		x = p->p.x;

		link_point *ll = NULL, *lu = NULL, *rl = NULL, *ru = NULL, *vu = NULL, *vl = NULL;
		lp_seg seg(p);
		while(true) { // Insert into left/right:
#ifdef DEBUG_DECOMPOSITION_CONSTRUCTx
			if(start_debug(seg.p1())) {
				cerr << " Updating l/r of " << *p << endl;
			}
#endif		
			// Update l/r:
			bool is_left = lr_insert(ll, lu, rl, ru, vl, vu, seg, p);
			update_sweep_line(sweep_line, seg, is_left);

			seg = lp_seg(p->next);
			is_left = lr_insert(ll, lu, rl, ru, vl, vu, seg, p);
			update_sweep_line(sweep_line, seg, is_left);

			++it_pq;
			if(it_pq == pq.end() || !(*p == **it_pq)) {
				break;
			}
			assert(*p == **it_pq);
			p = *it_pq;
			seg = lp_seg(p);
		}

#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
		if(start_debug(seg.p1())) {
			cerr << " Sweep line:" << endl;
			print_set(sweep_line);
		}

		if(!sweep_line.empty()) {
			lp_seg prev = *(sweep_line.begin());
			xycoord_t prev_y = prev.eval(x);
			for(lpsegpset::iterator it_sl = ++sweep_line.begin(); it_sl != sweep_line.end(); ++it_sl) {
				lp_seg sl_seg = *it_sl;
				if(prev_y > it_sl->eval(x)) {
					cerr << " SWEEP LINE INCONSISTENCY ON: " << *it_sl << " vs prev at y " << prev_y << endl;
					assert(false);
				}
				if(!sweep_line_cmp(prev,sl_seg)) {
					cerr << " SWEEP LINE INCONSISTENCY 2: " << endl;
					cerr << "  compare prev,sl_seg: " << sweep_line_cmp(prev,sl_seg) << endl;
					cerr << "  compare sl_seg,prev: " << sweep_line_cmp(sl_seg,prev) << endl;
					cerr << prev << endl;
					cerr << sl_seg << " at y " << prev_y << endl;
					assert(false);
				}
				if(sweep_line_cmp(sl_seg,prev)) {
					cerr << " SWEEP LINE INCONSISTENCY 3: " << endl;
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
		if(start_debug(seg.p1())) {
			cerr << "Handle points at " << *p << endl;
			if(ll != NULL)
				cerr << " LL: " <<  *ll << endl;
			if(lu != NULL)
				cerr << " LU: " <<  *lu << endl;
			if(rl != NULL)
				cerr << " RL: " <<  *rl << endl;
			if(ru != NULL)
				cerr << " RU: " <<  *ru << endl;
			if(vu != NULL)
				cerr << " VU: " <<  *vu << endl;
			if(vl != NULL)
				cerr << " VL: " <<  *vl << endl;
		}
#endif		

		/*
		  Link:
		  Top above:
		  Bottom below:
		  Left/Right up
		 */
		if(ll == NULL) { // only right:
			if(vl != NULL)
				rl = vl;
			if(vu != NULL)
				ru = vu;
			if(ru != rl)
				link_point::connect_vertical(ru, rl);						
		}
		else if(ru == NULL) { // only left:
			if(vl != NULL)
				ll = vl;
			if(vu != NULL)
				lu = vu;
			if(ll != lu)
				link_point::connect_vertical(lu, ll);						
		}
		else { // both above and below:
			if(ru != lu) { // above:
				if(vu != NULL) {
					if(vu != ru && ru->inside_up(parent)) {
						link_point::connect_vertical(vu, ru);
					}
					else if (vu != lu && lu->inside_up(parent)){
						link_point::connect_vertical(vu, lu);
					}
				}
				else if(ru->inside_up(parent)) { // normal transition above:
					if(ru->is_transitioning()) {
						link_point::connect_vertical(lu, ru);
					}
					else {
						link_point::connect_vertical(ru, lu);
					}
				}
			}
			if(rl != ll) { // below:
				if(vl != NULL) {
					if(vl != rl && rl->inside_down(parent)) {
						link_point::connect_vertical(rl, vl);
					}
					else if(vl != ll && ll->inside_down(parent)) {
						link_point::connect_vertical(ll, vl);
					}
				}
				else if(ll->inside_down(parent)) { 
					if(ll->is_transitioning()) {
						link_point::connect_vertical(ll, rl);
					}
					else {
						link_point::connect_vertical(rl, ll);
					}
				}
			}
		}

		// Connect up/down:
		lpsegpset::iterator up = sweep_line.upper_bound(seg), down = up;		
		while(up != sweep_line.end() && (up->p1() == p->p || up->p2() == p->p)) {
			++up;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug(seg.p1())) {
				cerr << " Pushing up up to ";
				if(up != sweep_line.end())
					cerr << *up << endl;
				else 
					cerr << "top" << endl;
			}
#endif
		}
		if(!sweep_line.empty()) {
			if(down != sweep_line.begin()) {
				--down;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
				if(start_debug(seg.p1())) {
					cerr << " Stepping down down to " << *down << endl;
				}
#endif
				while(down != sweep_line.begin() && ((p->p.x == down->p1().x && p->p.y <= down->p1().y) || (p->p.x == down->p2().x && p->p.y <= down->p2().y))) {
					--down;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
					if(start_debug(seg.p1())) {
						cerr << " Pushing down down to " << *down << endl;
					}
#endif
				}
				if(((p->p.x == down->p1().x && p->p.y <= down->p1().y) || (p->p.x == down->p2().x && p->p.y <= down->p2().y))) {
					down = sweep_line.end();
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
					if(start_debug(seg.p1())) {
						cerr << " Kill down" << endl;
					}
#endif
				}
			}
			else {
				down = sweep_line.end(); // Can't step down.
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
				if(start_debug(seg.p1())) {
					cerr << " No step down kill of down" << endl;
				}
#endif
			}
		}

		link_point *p_up = vu == NULL ? (ru == NULL ? lu : ru) : vu;
		link_point *p_down = vl == NULL ? (ll == NULL ? rl : ll) : vl;
		if(up != sweep_line.end() && p_up->inside_up(parent)) {
			lp_seg seg_up = *up;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug(seg.p1())) {
				cerr << " Shooting up. Removing " << *up << endl;
			}
#endif
			sweep_line.erase(up);
			lp_seg new_up = seg_up.split_from_below(p_up, parent);
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug(seg.p1())) {
				cerr << " Inserting " << new_up << endl;
			}
#endif
			pair<lpsegpset::iterator,bool> pair = sweep_line.insert(new_up);
			assert(pair.second);				
		}
		if(down != sweep_line.end() && p_down->inside_down(parent)) {
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug(seg.p1())) {
				cerr << " Shooting down " << endl;
			}
#endif
			lp_seg seg_down = *down;
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug(seg.p1())) {
				cerr << " Removing " << seg_down << endl;
			}
#endif
			sweep_line.erase(down);

			seg_down = seg_down.split_from_above(p_down, parent);
#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
			if(start_debug(seg.p1())) {
				cerr << " Inserting " << seg_down << endl;
			}
#endif
			pair<lpsegpset::iterator,bool> pair = sweep_line.insert(seg_down);
			assert(pair.second);
		}
	}

#ifdef DEBUG_DECOMPOSITION_CONSTRUCT
	if(start_debug()) {
		cerr << "Decomposition construction done " << endl;
		print();
	}
#endif
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug())
		print();
#endif
//	print();
}

decomposition::~decomposition() {
	for(vector<link_point*>::iterator it = points.begin(); it != points.end(); ++it) {
		link_point *first = *it;
		link_point *next = first->next;
		do {
			link_point *nn = next->next;
			delete next;
			next = nn;
		}
		while(next != first);
		delete first;
	}
	points.clear();
}

lp_seg lp_seg::seg_above() {
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug(p1())) {
		cerr << "  Seg above calculation for: " << *this << endl;
	}
#endif
	xycoord_t const x1 = p1().x;
	xycoord_t const x2 = p2().x;
	assert(x1 != x2);

	xycoord_t x;
	link_point *r = lp1();
	do {
		r = r->up;
	}
	while(r != NULL && r->next->p.x != x2 && r->prev->p.x != x2);
	link_point *candidate1 = r;

	r = lp2();
	do {
		r = r->up;
	}
	while(r != NULL && r->next->p.x != x1 && r->prev->p.x != x1);
	link_point *candidate2 = r;

	if(candidate1 == NULL) {
		if(candidate2 == NULL && lp2()->prev->p.x == x1 && lp2()->prev->p.y > p1().y) {
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1())) {
				cerr << "  Super special case bail: " << *this << endl;
			}
#endif
			return lp_seg(lp2());
		}

		assert(candidate2 != NULL);
		r = candidate2;
		x = x1;
	}
	else if(candidate2 == NULL) {
		if(candidate1 == NULL && lp1()->next->p.x == x2 && lp1()->next->p.y > p2().y) {
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1())) {
				cerr << "  Super special case bail2: " << *this << endl;
			}
#endif
			return lp_seg(lp1()->next);
		}

		assert(candidate1 != NULL);
		r = candidate1;
		x = x2;
	}
	else if(candidate1->prev == candidate2) {
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			cerr << " Perfect bail: " << lp_seg(candidate1) << endl;
		}
#endif
		return lp_seg(candidate1);
	}
	else if(candidate1->next == candidate2) {
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			cerr << " Perfect bail: " << lp_seg(candidate2) << endl;
		}
#endif
		return lp_seg(candidate2);		
	}
	else {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1())) {
			cerr << " WARNING ON UP! " << endl;
		}
#endif
		link_point *c1o = candidate1->next;
		if(c1o->p.x != x2)
			c1o = candidate1->prev;
		
		// Choose lowest!
		r = c1o->p.y < candidate2->p.y ? candidate1 : candidate2;
		x = c1o->p.y < candidate2->p.y ? x2 : x1;
	}
#ifdef DEBUG_DECOMPOSITION_WALKx
	if(start_debug(p1())) {
		cerr << " Above found: " << (r->prev->p.x == x ? lp_seg(r) : lp_seg(r->next)) << endl;
	}
#endif
	return r->prev->p.x == x ? lp_seg(r) : lp_seg(r->next);
}

lp_seg lp_seg::seg_below() {
#ifdef DEBUG_DECOMPOSITION_WALKx
	if(start_debug(p1())) {
		cerr << "  Seg below calculation for: " << *this << endl;
		cerr << " pa: " << *lp1() << endl;	
		cerr << " pb: " << *lp2() << endl;	
	}
#endif	
	xycoord_t x1 = p1().x;
	xycoord_t x2 = p2().x;
	assert(x1 != x2);

	link_point *r = lp1();
	xycoord_t x;
	do {
		r = r->down;
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			if(r != NULL)
				cerr << " down a: " << *r << endl;	
		}
#endif	
	}
	while(r != NULL && r->next->p.x != x2 && r->prev->p.x != x2);
	link_point *candidate1 = r;

	r = lp2();
	do {
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			if(r != NULL)
				cerr << " down b: " << *r << endl;	
		}
#endif	
		r = r->down;
	}
	while(r != NULL && r->next->p.x != x1 && r->prev->p.x != x1);
	link_point *candidate2 = r;

	if(candidate1 == NULL) {
		assert(candidate2 != NULL);
		r = candidate2;
		x = x1;
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			cerr << " Candidate 1 i NULL" << endl;
		}
#endif
	}
	else if(candidate2 == NULL) {
		assert(candidate1 != NULL);
		r = candidate1;
		x = x2;
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			cerr << " Candidate 2 i NULL" << endl;
		}
#endif
	}
	else if(candidate1->prev == candidate2) {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1())) {
			cerr << " Perfect bail: " << lp_seg(candidate1) << endl;
		}
#endif
		return lp_seg(candidate1);
	}
	else if(candidate1->next == candidate2) {
#ifdef DEBUG_DECOMPOSITION_WALKx
		if(start_debug(p1())) {
			cerr << " Perfect bail: " << lp_seg(candidate2) << endl;
		}
#endif
		return lp_seg(candidate2);		
	}
	else {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1())) {
			cerr << " WARNING ON DOWN! " << endl;
			cerr << " candidate 1: " << *candidate1 << endl;
			cerr << " candidate 2: " << *candidate2 << endl;
		}
#endif
		link_point *c1o = candidate1->next; // Candidate 1 - other
		if(c1o->p.x != x2)
			c1o = candidate1->prev;
		assert(c1o->p.x == x2);
		
		// Choose lowest!
		r = c1o->p.y > candidate2->p.y ? candidate2 : candidate1;
		x = c1o->p.y > candidate2->p.y ? x1 : x2;
	}
#ifdef DEBUG_DECOMPOSITION_WALKx
	if(start_debug(p1())) {
		cerr << " Below found: " << (r->prev->p.x == x ? lp_seg(r) : lp_seg(r->next)) << endl;
	}
#endif
	return r->prev->p.x == x ? lp_seg(r) : lp_seg(r->next);
}

/*  
  Returns true if p from the line ppm is inside of the trapezoid with this as the lower 
  segment. That is. ppm must have some part interior to the trapezoid.
 */
bool lp_seg::contains_p_above(contour_point p, contour_point pm) {
	contour_point left = p1();
	contour_point right = p2();
	if(left.x > right.x)
		swap(left,right);
	lp_seg above = seg_above();
	contour_point above_left = above.p1();
	contour_point above_right = above.p2();
	if(above_left.x > above_right.x)
		swap(above_left,above_right);

	if(left.x > p.x || right.x < p.x || eval(p.x) > p.y || above.eval(p.x) < p.y ||
	   (left.x == p.x && pm.x < p.x) || (right.x == p.x && pm.x > p.x)) {
		return false;
	}

	if(left.x == p.x) {
		if(pm.x == p.x) {
			return eval(pm.x) <= pm.y || above.eval(pm.x) >= pm.y;
		}
		if(left == p) {
			if(left.x < pm.x)
				return eval(pm.x) <= pm.y;
			else
				return right.y <= lp_seg(p,pm).eval(right.x);
		}
		if(above_left == p) {
			if(right.x < pm.x)
				return above.eval(pm.x) >= pm.y;
			else
				return above_right.y >= lp_seg(p,pm).eval(above_right.x);
		}
	}
	else if(right.x == p.x) {
		if(pm.x == p.x) {
			return eval(pm.x) <= pm.y || above.eval(pm.x) >= pm.y;
		}
		if(right == p) {
			if(left.x < pm.x)
				return eval(pm.x) <= pm.y;
			else
				return left.y <= lp_seg(p,pm).eval(left.x);
		}
		if(above_right == p) {
			if(left.x < pm.x)
				return above.eval(pm.x) >= pm.y;
			else
				return above_left.y >= lp_seg(p,pm).eval(above_left.x);
		}
	}
	return true;
}

link_point* decomposition::bfs_find(contour_point p,contour_point pm,link_point* guide) {
	bfs_scans++;
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug(p)) {
		cerr << "Warning: bfs scan " << bfs_scans << " for " << p << "-" << pm << " from " << *guide << endl;
	}
#endif
	lp_seg seg(guide);
	if(seg.is_vertical() || !seg.inside_up(parent_contour)) {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p)) {
			cerr << " Bad start: linear scan! " << endl;
		}
#endif
		return find(p,pm);
	}
	queue<link_point*> q;
	q.push(guide);
	set<link_point*> seen;

	int iterations = 0;
	while(!q.empty()) {
		link_point *llp = q.front();
		q.pop();
		++iterations;
		bfs_steps++;
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug()) {
			cerr << " BFS iteration " << *llp << endl;
		}
#endif
				
		lp_seg seg(llp);
		if(seg.contains_p_above(p,pm)) {
#ifdef DEBUG_DECOMPOSITION_WALK
			if(iterations >= 100) {
				cerr << " BFS found " << *llp << " in iteration " << iterations << ", BFS " << bfs_scans << endl;
			}
			if(start_debug()) {
				cerr << " BFS found " << *llp << " in iteration " << iterations << endl;
			}
#endif
			return llp;
		}

		link_point* above[2] = {seg.lp1(),seg.lp2()};
		for(int i = 0; i < 2; i++) {			
			link_point* lp = above[i];
			while(lp->up != NULL) {
				lp = lp->up;
			}
			
			while(lp != NULL) {
				lp_seg seg2(lp);
				if(lp != llp && !seg2.is_vertical() && seg2.inside_up(parent_contour) && seen.insert(lp).second) {
					q.push(lp);
				}
				seg2 = lp_seg(lp->next);
				if(lp->next != llp && !seg2.is_vertical() && seg2.inside_up(parent_contour) && seen.insert(lp->next).second) {
					q.push(lp->next);
				}
				lp = lp->down;
			}
		}
	}
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug()) {
		cerr << " Bad BFS: linear scan! " << endl;
	}
#endif
	return find(p,pm);
}

// finds self and prev of seg below p.
link_point* decomposition::find(contour_point p,contour_point pm) {
	linear_scans++;
#ifdef DEBUG_DECOMPOSITION_WALK
	if(linear_scans > 1) {
		cerr << "Error: too many linear scans: " << linear_scans << " for " << p << endl;
		assert(false);
	}
	if(start_debug(p)) {
		cerr << " Starting linear scan for p " << p << endl;
	}
#endif
	for(vector<link_point*>::iterator it = points.begin(); it != points.end(); ++it) {
		link_point *first = *it, *lp = first;
		do {
			lp_seg seg(lp);
			if(!seg.is_vertical() && seg.inside_up(parent_contour) && seg.contains_p_above(p,pm)) {
#ifdef DEBUG_DECOMPOSITION_WALK
				if(start_debug(p)) {
					cerr << "  Found: " << seg << endl;
					cerr << "  Found (lp): " << *lp << endl;
				}
#endif
				return lp;
			}
			lp = lp->next;
		}
		while(lp != first);
	}
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug(p)) {
		cerr << "Linear scan fail!" << endl;
	}
#endif
	return NULL;
}

bool lp_seg::intersects(contour_point &lp1, contour_point &lp2) {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1())) {
			cerr << " Test intersects of " << *this << " vs " << lp1 << "-" << lp2 << endl;
		}
#endif
	// top & bottom line:
	if(intersects_line(lp1,lp2)) {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1())) {
			cerr << " Intersecting a: " << lp1 << "-" << lp2 << endl;
		}
#endif
		return false;
	}	
	lp_seg above = seg_above();
	if(above.intersects_line(lp1,lp2)) {
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1())) {
			cerr << " Intersecting b: " << above << endl;
		}
#endif
		return false;
	}	
	// left:
	xycoord_t x_min = max(min(lp1.x,lp2.x),min(p1().x,p2().x));	
	xycoord_t y_left_min = eval(x_min);
	xycoord_t y_left_max = above.eval(x_min);
	xycoord_t p12y = get_y(lp1,lp2,x_min);
	if(y_left_min > p12y || y_left_max < p12y)
		return false; 
    // right:
	xycoord_t x_max = min(max(lp1.x,lp2.x),max(p1().x,p2().x));	
	xycoord_t y_right_min = eval(x_max);
	xycoord_t y_right_max = above.eval(x_max);
	p12y = get_y(lp1,lp2,x_max);
	if(y_right_min > p12y || y_right_max < p12y)
		return false; 
	return true;
}

bool lp_seg::intersects_line(contour_point &lp1, contour_point &lp2) {
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug(lp1)) {
		cerr << " Intersects line: " << *this << " vs " << lp1 << " - " << lp2 << endl;
	} 
#endif
	contour_point l = p1();
	contour_point r = p2();
	if(l.x > r.x)
		swap(l,r);
	contour_point ll = lp1;
	contour_point lr = lp2;
	if(ll.x > lr.x)
		swap(ll,lr);
	
	if(lr.x < l.x || ll.x > r.x) {
		return false; // out of x bound.
	}

	xycoord_t x1 = l.x;
	xycoord_t y1 = l.y;
	xycoord_t x2 = r.x;
	xycoord_t y2 = r.y;
	xycoord_t x3 = ll.x;
	xycoord_t y3 = ll.y;
	xycoord_t x4 = lr.x;
	xycoord_t y4 = lr.y;

	int order_l;
	if(x3 < x1) {
		float det341 = (x3 - x1)*(y4 - y1) - (x4 - x1)*(y3 - y1);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug()) {
			cerr << " det341: " << det341 << endl;
		}
#endif
		order_l = -simplification::int_sign(det341, SQRT_EPSILON);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(order_l > 0) {
			assert(y1 < lp_seg(ll,lr).eval(x1));
		}
		else if(order_l < 0) {
			assert(y1 > lp_seg(ll,lr).eval(x1));
		}
#endif
	}
	else {
		float det123 = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug()) {
			cerr << " det123: " << det123 << endl;
		}
#endif
		order_l = simplification::int_sign(det123, SQRT_EPSILON);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(order_l < 0) {
			assert(y3 < eval(x3));
		}
		else if(order_l > 0) {
			assert(y3 > eval(x3));
		}
#endif
	}

	int order_r;
	if(x4 > x2) {
		float det342 = (x3 - x2)*(y4 - y2) - (x4 - x2)*(y3 - y2);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug()) {
			cerr << " det342: " << det342 << endl;
		}
#endif
		order_r = -simplification::int_sign(det342, SQRT_EPSILON);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(order_r > 0) {
			assert(y2 < lp_seg(ll,lr).eval(x2));
		}
		else if(order_r < 0) {
			assert(y2 > lp_seg(ll,lr).eval(x2));
		}
#endif
	}
	else {
		float det124 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug()) {
			cerr << " det124: " << det124 << endl;
		}
#endif
		order_r = simplification::int_sign(det124, SQRT_EPSILON);
#ifdef DEBUG_DECOMPOSITION_WALK
		if(order_r < 0) {
			assert(y4 < eval(x4));
		}
		else if(order_r > 0) {
			assert(y4 > eval(x4));
		}
#endif
	}
	
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug()) {
		cerr << " orderl: " << order_l << ", order_r: " << order_r << endl;
	}
#endif
	return order_l == 0 || order_r == 0 || order_l != order_r;
}

bool decomposition::contains_line(vector<contour_point> &points, int p1i, int p2i) {
	if(points.empty())
		return true;
	contour_point p1 = points[p1i];
	contour_point pm = points[p1i+1];
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug()) {
		cerr << " Contains line for points between " << p1 << " and " << points[p2i] << endl;
		cerr << "  Guide is ";
		if(guide == NULL)
			cerr << "Null" << endl;
		else
			cerr << lp_seg(guide) << endl;
	}
#endif
	if(guide == NULL) {
		guide = find(p1,pm);
	}
	else if(!lp_seg(guide).contains_p_above(p1, pm)) {
		guide = bfs_find(p1,pm,guide);
	}
	link_point *lp = guide;
	if(guide == NULL)
		return false;

	// Build can. seq:
	vector<link_point*> canonical_sequence;
	canonical_sequence.push_back(guide);
	for(int i = p1i+1; i <= p2i; i++) {
		vector<link_point*> to_add;
		bool ok = walk_contour(points[i-1], points[i], lp, to_add);
		if(!ok)
			return false;
		for(vector<link_point*>::iterator it2 = to_add.begin(); it2 != to_add.end(); ++it2) {
			int cs = canonical_sequence.size();
			if(cs > 2 && canonical_sequence[cs-2] == *it2) {
				canonical_sequence.pop_back();
			}
			else {
				canonical_sequence.push_back(*it2);
			}
		}
	}
	
	// test CS:
	contour_point p2 = points[p2i];
	for(vector<link_point*>::iterator it = canonical_sequence.begin(); it != canonical_sequence.end(); ++it) {
		assert(*it != NULL);
		lp_seg ls(*it);
		if(!ls.intersects(p1,p2)) {
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug()) {
				cerr << " Intersecting c: " << ls << endl;
			}
#endif
			return false;
		}
	}   
	guide = lp;
	return true;
}

/* 
   Returns the (non-can.) seq. of the trapezoids visited (excluding start tz).
   Loop invariant: lp_seg(below,true) is seg below.
 */
bool decomposition::walk_contour(contour_point &p1, contour_point &p2, 
								 link_point* &start, vector<link_point*> &v) {
	assert(start != NULL);
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug(p1)) {
		cerr << "LINE: " << p1 << ", " << p2 << endl;
	}
#endif
	if(!lp_seg(start).contains_p_above(p1,p2)) {
		start = bfs_find(p1,p2,start);
		if(start == NULL) {
			return false;
		}
		v.push_back(start);
	}

	bool going_right = p2.x > p1.x || 
		(p1.x == p2.x && p1.x == max(start->p.x, start->prev->p.x));
	link_point *below = start;
	while(!lp_seg(below).contains_p_above(p2, p1)) {
		if(min(lp_seg(below).p1().x,lp_seg(below).p2().x) > max(p1.x,p2.x) ||
		   max(lp_seg(below).p1().x,lp_seg(below).p2().x) < min(p1.x,p2.x)) {
			cerr << "Error: Walking out of decomposition: Line " << p1 << "-" << p2 << " walking out of " << lp_seg(below) << endl;
			assert(false);
			return false;
		}
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1)) {
			cerr << " Walking from " << lp_seg(below) << " up " << lp_seg(below).seg_above() << endl;
		}
#endif		
		link_point *r = below;

		assert(r->prev->p.x != r->p.x); // not vertical.
		if(going_right == (r->prev->p.x > r->p.x)) {
			r = r->prev; // r is right of below.
		}
		xycoord_t rx = r->p.x;

		// Walk:
		// walk up:
		xycoord_t p12y = get_y(p1,p2,rx);
		do {
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1)) {
				cerr << "  Walking up from " << *r;
			}
			if(r->up == NULL) {
				cerr << "r->up == NULL Error: " << *r << endl;
			}
#endif		
			assert(r->up != NULL);
			r = r->up;
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1)) {
				cerr << " to " << *r << endl;
			}
#endif		
		} 
		while(r->p.y < p12y || 
			  (going_right && r->prev->p.x <= rx && r->next->p.x <= rx) || 
			  (!going_right && r->prev->p.x >= rx && r->next->p.x >= rx));
		if(r->p.y == p12y) {
#ifdef DEBUG_DECOMPOSITION_WALK
			cerr << "  ERROR: WHERE TO GO? " << endl;
#endif		
			return false; // Where to go?
		}
#ifdef DEBUG_DECOMPOSITION_WALK
		if(start_debug(p1)) {
			cerr << "  above of next: " << *r << " vs line_y: " << p12y << endl;
		}
#endif		
		// find r_down:		
		assert(!(r->prev->p.x == r->p.x && r->next->p.x == r->p.x));
		if(going_right) {
			lp_seg above(r->prev->p.x > rx ? r :r->next);
			if(above.inside_up(parent_contour)) {
				above = lp_seg(r->prev->p.x <= rx ? r : r->next);
			}
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1)) {
				cerr << "  Real above: " << above << endl;
			}
#endif		
			assert(above.p1().x != above.p2().x);
			lp_seg s_below = above.seg_below();
			
			below = s_below.lp1();//below = s_below.p1().x > s_below.p2().x ? s_below.lp1() : s_below.lp2();
			v.push_back(below);
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1)) {
				cerr << "  New below: " << s_below << endl;
			}
#endif		
			assert(s_below.inside_up(parent_contour));
		}
		else {
			lp_seg above(r->prev->p.x < rx ? r :r->next);
			if(above.inside_up(parent_contour))
				above = lp_seg(r->prev->p.x >= rx ? r : r->next);
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1)) {
				cerr << "  Real above: " << above << endl;
			}
#endif		
			assert(above.p1().x != above.p2().x);
			lp_seg s_below = above.seg_below();
			
			below = s_below.lp1();//s_below.p1().x < s_below.p2().x ? s_below.lp1() : s_below.lp2();
			v.push_back(below);
#ifdef DEBUG_DECOMPOSITION_WALK
			if(start_debug(p1)) {
				cerr << "  New below: " << s_below << endl;
			}
#endif		
			assert(s_below.inside_up(parent_contour));
		}
	}

/*	// No touching!
	s1 = lp_seg(below);
	s2 = s1.seg_above();
	if(s1.p1()==p2 || s1.p2()==p2 || s1.p1()==p2 || s2.p2()==p2) {
		return false;
	}*/

	start = below;
#ifdef DEBUG_DECOMPOSITION_WALK
	if(start_debug(p1)) {
		cerr << "  OK!: " << *start << endl;
	}
#endif		

	return true;
}

void decomposition::print() {
	cerr << "Decomposition:" << endl;
	for(vector<link_point*>::iterator it = points.begin(); it != points.end(); ++it) {
		cerr << " Contour with points: " << endl;
		link_point *first = *it;
		link_point *lp = first;
		do {
			assert(lp != NULL);
			cerr << "  " << *lp << endl;			
			assert(lp->prev != NULL);
			assert(lp->next != NULL);
			assert(lp->prev->next == lp);
			assert(lp->next->prev == lp);
			assert(lp->prev != lp);
			assert(lp->next != lp);
			assert(lp->up != lp);
			assert(lp->down != lp);
			if(lp->up != NULL) {
				assert(lp->up->down == lp);
				assert(lp->up->p.y >= lp->p.y);				
				assert(lp->up->p.x == lp->p.x);				
				assert(lp->up->up != lp);
			}
			if(lp->down != NULL) {
				assert(lp->down->up == lp);
				assert(lp->down->p.x == lp->p.x);
				assert(lp->down->p.y <= lp->p.y);
				assert(lp->down->down != lp);
			}
			lp = lp->next;
		}
		while(lp != first);
	}
	cerr << endl;
}

}

