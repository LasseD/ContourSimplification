// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_TOPOLOGY_H__
#define __TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_TOPOLOGY_H__
#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include "contour_types.h"
#include <algorithm>
#include <vector>

namespace terrastream{

template<typename T>
void print(stream<T> &s) {
	std::cerr << "|segs| = " << s.stream_len() << std::endl;	
	T *t, prev;	
	while (s.read_item(&t) == NO_ERROR) {
		if(t->label != prev.label)
			std::cerr << *t << std::endl;
		prev = *t;
	}
	s.seek(0);
}

void relabel(stream<ranked_labelled_signed_contour_segment> &ls,
			 stream<ranked_labelled_signed_contour_segment> &out);

void build_topology(stream<ranked_labelled_signed_contour_segment> &ls,stream<topology_edge> &topology);

void order_for_simplification(stream<topology_edge> &topology,
							  stream<contour_point> &contours_in,
							  stream<contour_point> &contours_out);

struct pq_entry {
	int p, lv, c;
	
	pq_entry(int _p, int _lv, int _c);
	pq_entry();
	pq_entry(const pq_entry &p);
	
	friend std::ostream& operator << (std::ostream &ostr, const pq_entry &l) {
		return ostr << "P " << l.p << ", lv " << l.lv << ", c " << l.c;
	}

	bool operator<(const pq_entry &l) const{
		if(p != l.p)
			return p < l.p;
		return c < l.c;
	}
};

template<typename T>
class expand_set {
	std::vector<T> v;	// bool gone_io; if iov is not NULL.
	stream<T> *iov;
	size_t capacity;
	bool sorted;
	bool first; // for next
	typename std::vector<T>::iterator itv;

public:	
	expand_set(size_t cap) {
		capacity = cap;
		sorted = false;
		iov = NULL;
		first = true;
	}

	~expand_set() {
		if(iov != NULL) {
			iov->truncate(0);
			delete iov;
		}
		else
			v.clear();
	}

	void insert(T const &t) {
		assert(!sorted);
		if(v.size() == capacity) { // Switch tech:
			iov = new stream<T>();
			typename std::vector<T>::iterator it;
			for(it = v.begin(); it != v.end(); ++it) {
				iov->write_item(*it);
			}
			v.clear();
		}
		if(iov != NULL) {
			iov->write_item(t);
		}
		else {
			v.push_back(t);
		}
	}

	template<typename Compare>
	void sort(Compare comp) {
		if(iov == NULL) {
			std::sort(v.begin(), v.end(), comp);
		}
		else {
			ts_sort(iov,&comp);
		}
		sorted = true;
	}
	
	bool next(T &t) {
		assert(sorted);
		if(iov == NULL) {
			if(first) {
				itv = v.begin();
				first = false;
			}
			if(itv == v.end())
				return false;
			t = *itv;
			++itv;
			return true;
		}
		else {
			T *item;
			err res = iov->read_item(&item);
			if(res != NO_ERROR)
				return false;
			t = *item;
			return true;
		}
	}
};

}
#endif /*__TEST_CONTOUR_SIMPLIFICATION_IO_CONTOURS_TOPOLOGY_H__*/
