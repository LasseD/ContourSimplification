// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CONTOUR_TYPES_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CONTOUR_TYPES_H__
#include <terrastream/common/common.h>
#include <terrastream/common/tflow_types.h>

namespace terrastream{

  template<class T>
  struct ptr_cmp {
	bool operator() (const T* a, const T* b) const {
	  return (*a) < (*b);
	}
  };

struct segment_point {
	public:
	xycoord_t x, y;
	
	segment_point() {}
	segment_point(xycoord_t _x, xycoord_t _y) : x(_x), y(_y) {}
	segment_point(const segment_point &l) : x(l.x), y(l.y) {}
	
	bool operator==(const segment_point &s) const {
		return x == s.x && y == s.y;
	}
};

struct contour_point : segment_point {
	int rank, label;

	contour_point();
	contour_point(xycoord_t _x,xycoord_t _y,int r, int lbl);
	contour_point(const contour_point &p);

	friend std::ostream& operator << (std::ostream &ostr, const contour_point &l) {
		return ostr << "(" << l.rank << "@" << l.label << ":" << l.x << "," << l.y << ")";
	}
	bool operator==(const contour_point &s) const {
		return x == s.x && y == s.y;
	}
	bool operator!=(const contour_point &s) const {
		return !(*this == s);
	}
    // For Priority Queue of contains_intersections.
	bool operator<(const contour_point &s) const { 
		if(x != s.x)
			return x < s.x;
		if(y != s.y)
			return y < s.y;
		if(label != s.label)
			return label < s.label;
		return rank < s.rank;
	}
	double length() const;
};

contour_point operator-(const contour_point& a, const contour_point& b);
double dot(const contour_point& a, const contour_point& b);
contour_point operator*(const contour_point& a, const double scale);
contour_point operator*(const double scale, const contour_point& a);
contour_point operator+(const contour_point& a, const contour_point& b);
bool p1_above_l(const segment_point &p1, const segment_point &l1, const segment_point &l2);

  //The endpoints of the contour_segment are sorted such that x1<x2, or in case of tie, such that y1<y2.
  //The z-value is the height of the contour
  struct contour_segment {
	xycoord_t x1,y1,x2,y2;
	elev_t z;
	contour_segment() {}
	contour_segment(xycoord_t _x1,xycoord_t _y1,xycoord_t _x2,xycoord_t _y2,elev_t _z);

 const segment_point p1() const;
 const segment_point p2() const;

	bool operator<(const contour_segment &l) const{
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

	bool operator<=(const contour_segment &l) const{
		return operator<(l) || operator==(l);
	}

	bool operator==(const contour_segment &l) const{
		return x1==l.x1 && y1==l.y1 && x2==l.x2 && y2==l.y2;
	}

	friend std::ostream& operator << (std::ostream &ostr, const contour_segment &l) {
	  return ostr << "(" << l.x1 << "," << l.y1 << ") -> (" << l.x2 << "," << l.y2 << ") z=" << l.z;
	}
  };

  //This type of contour has an extra boolean for storing information
  struct signed_contour_segment : contour_segment {
	bool sign;
	signed_contour_segment() {};
	signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign=false);

	friend std::ostream& operator << (std::ostream &ostr, const signed_contour_segment &l) {
		ostr << "(" << l.x1 << "," << l.y1 << ") -> (" ;
		if (l.x1 == l.x2) {
			ostr << "-||-,";
		} else {
			ostr << l.x2 << ",";
		}
		if (l.y1 == l.y2) {
			ostr << "-||-";
		} else {
			ostr << l.y2;
		}
		return ostr << ") z=" << l.z << " sign=" << l.sign;
	}
  };

  struct labelling_signed_contour_segment{
	  segment_point p[2];
	  elev_t z;
	  bool sign;
	  labelling_signed_contour_segment() {};
	  labelling_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t _z,bool _sign=false);
	  friend std::ostream& operator << (std::ostream &ostr, const labelling_signed_contour_segment &l) {
		return ostr << "(" << l.p[0].x << "," << l.p[0].y << ") -> (" << l.p[1].x << "," << l.p[1].y << ") z=" << l.z << " sign=" << l.sign;
	  }
  };


  struct labelled_signed_contour_segment : signed_contour_segment{
	  int label;
	  labelled_signed_contour_segment() {};
	  labelled_signed_contour_segment(labelling_signed_contour_segment &s,int lbl);
	  labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign,int lbl);
	  friend std::ostream& operator << (std::ostream &ostr, const labelled_signed_contour_segment &l) {
		return ostr << "(" << l.x1 << "," << l.y1 << ") -> (" << l.x2 << "," << l.y2 << ") z=" << l.z << " sign=" << l.sign << " label=" << l.label;
	  }
  };

  struct ranked_labelled_signed_contour_segment : labelled_signed_contour_segment{
	int rank;
	ranked_labelled_signed_contour_segment() {};
	ranked_labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign,int label,int rnk);
	friend std::ostream& operator << (std::ostream &ostr, const ranked_labelled_signed_contour_segment &l) {
		ostr << "(" << l.x1 << "," << l.y1 << ") -> (" ;
		if (l.x1 == l.x2) {
			ostr << "-||-,";
		} else {
			ostr << l.x2 << ",";
		}
		if (l.y1 == l.y2) {
			ostr << "-||-";
		} else {
			ostr << l.y2;
		}
		return ostr << ") z=" << l.z << " sign=" << l.sign << " rank=" << l.rank << " label=" << l.label;
	}
  };

  struct triangle_point{
	xycoord_t x,y;
	elev_t z;

	//Sorts the triangle points by x, using y as tiebreaker
	bool operator<(const triangle_point &r) const{
		return x != r.x ? x < r.x : y < r.y;
	}

	triangle_point() {}
	triangle_point(xycoord_t _x,xycoord_t _y,elev_t _z);
	  
	  friend std::ostream& operator << (std::ostream &ostr, const triangle_point &t) {
		  return ostr << "(" << t.x << "," << t.y << "," << t.z << ")";
	  }
  };

  struct triangle{
	triangle_point points[3];
	triangle() {}
	triangle(triangle_point p1,triangle_point p2,triangle_point p3);
  };

  struct map_info{
	xycoord_t minX,maxX,minY,maxY;
	elev_t minZ,maxZ;
  };

  struct topology_edge{
	  int c, p;
	  elev_t c_z;//, p_z;
	  topology_edge() {}
	  topology_edge(int child,int parent,elev_t cz) : c(child), p(parent), c_z(cz) {}
	  friend std::ostream& operator << (std::ostream &ostr, const topology_edge &l) {
		  return ostr << "c:" << l.c << "@" << l.c_z << ", p:" << l.p;// << "@" << l.p_z;
	  }
  };

  struct decomposition_triangle : triangle {
	decomposition_triangle* neighbours[3];
	decomposition_triangle() {}
	decomposition_triangle(triangle_point p1,triangle_point p2,triangle_point p3);
  };

  struct triangulated_ranked_labelled_signed_contour_segment : ranked_labelled_signed_contour_segment{
	decomposition_triangle const* triangle;
	triangulated_ranked_labelled_signed_contour_segment() {};
	triangulated_ranked_labelled_signed_contour_segment(xycoord_t x1,xycoord_t y1,xycoord_t x2,xycoord_t y2,elev_t z,bool sign,int label,int rnk,decomposition_triangle* triangle);
	friend std::ostream& operator << (std::ostream &ostr, const triangulated_ranked_labelled_signed_contour_segment &l) {
	  return ostr << "(" << l.x1 << "," << l.y1 << ") -> (" << l.x2 << "," << l.y2 << ") z=" << l.z <<
		  " sign=" << l.sign << " label=" << l.label << " rank=" << l.rank << " triangle=" << l.triangle;
	}
  };

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_CONTOUR_TYPES_H__*/
