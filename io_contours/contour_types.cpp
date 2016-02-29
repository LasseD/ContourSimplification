// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "contour_types.h"

using namespace std;

namespace terrastream {

contour_segment::contour_segment(xycoord_t _x1,xycoord_t _y1,
								 xycoord_t _x2,xycoord_t _y2,elev_t _z) :
  x1(_x1),y1(_y1),x2(_x2),y2(_y2),z(_z) {
	if ((x1>x2) || (x1==x2 && y1>y2)) {
	swap(x1,x2);
	swap(y1,y2);
  }
}

contour_point::contour_point(xycoord_t _x,xycoord_t _y,int r, int lbl) :
	segment_point(_x,_y),label(lbl),rank(r) {
}
contour_point::contour_point() : label(-1) {
}
contour_point::contour_point(const contour_point &p) : segment_point(p.x,p.y), label(p.label), rank(p.rank) {
}

const segment_point contour_segment::p1() const {
	return segment_point(x1,y1);
}
const segment_point contour_segment::p2() const {
	return segment_point(x2,y2);
}

bool p1_above_l(const segment_point &p1, 
				const segment_point &l1,
				const segment_point &l2) {
	if(p1 == l1 || p1 == l2)
		return false;
	
	xycoord_t x3 = p1.x;
	xycoord_t y3 = p1.y;
	xycoord_t x1 = l1.x;
	xycoord_t y1 = l1.y;
	xycoord_t x2 = l2.x;
	xycoord_t y2 = l2.y;

	float det123 = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
	return det123 > 0;
}

contour_point operator-(const contour_point& a, const contour_point& b) {
	return contour_point(a.x-b.x,a.y-b.y,a.rank,a.label);
}
contour_point operator+(const contour_point& a, const contour_point& b) {
	return contour_point(a.x+b.x,a.y+b.y,a.rank,a.label);
}
double dot(const contour_point& a, const contour_point& b) {
	return a.x*b.x + a.y*b.y;
}
contour_point operator*(const contour_point& a, const double scale) {
	contour_point b = a;
	b.x*=scale;
	b.y*=scale;
	return b;
}
contour_point operator*(const double scale, const contour_point& a) {
	contour_point b = a;
	b.x*=scale;
	b.y*=scale;
	return b;
}
double contour_point::length() const {
	return (double)sqrt(x*x + y*y);
}


signed_contour_segment::signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,
								xycoord_t y2,elev_t z,bool _sign) :
  contour_segment(x1,y1,x2,y2,z),sign(_sign) {
}

labelling_signed_contour_segment::labelling_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t _z,bool _sign) : z(_z), sign(_sign) {
	if (x1>x2 || (x1==x2 && y1>y2)) {
		swap(x1,x2);
		swap(y1,y2);
	}
	p[0] = segment_point(x1, y1);
	p[1] = segment_point(x2, y2);
}

labelled_signed_contour_segment::labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign,int lbl) :
	signed_contour_segment(x1,y1,x2,y2,z,sign),label(lbl) {
}

ranked_labelled_signed_contour_segment::ranked_labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign,int lbl,int rnk) :
  labelled_signed_contour_segment(x1,y1,x2,y2,z,sign,lbl),rank(rnk) {
}

labelled_signed_contour_segment::labelled_signed_contour_segment(labelling_signed_contour_segment &s,int lbl) :
	signed_contour_segment(s.p[0].x,s.p[0].y,s.p[1].x,s.p[1].y,s.z,s.sign),label(lbl) {
}

triangle_point::triangle_point(xycoord_t _x,xycoord_t _y,elev_t _z) : x(_x),y(_y),z(_z) {}

triangle::triangle(triangle_point p1,triangle_point p2,triangle_point p3) {
  points[0]=p1;
  points[1]=p2;
  points[2]=p3;
}

triangulated_ranked_labelled_signed_contour_segment::triangulated_ranked_labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign,int lbl,int rnk,decomposition_triangle* tri) :
  ranked_labelled_signed_contour_segment(x1,y1,x2,y2,z,sign,lbl,rnk),triangle(tri) {
}

decomposition_triangle::decomposition_triangle(triangle_point p1,triangle_point p2,triangle_point p3) {
  points[0]=p1;
  points[1]=p2;
  points[2]=p3;
  neighbours[0]=neighbours[1]=neighbours[2]=0;
}

} //namespace terrastream
