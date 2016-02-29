// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#ifndef __TEST_CONTOUR_SIMPLIFICATION_UTIL_H__
#define __TEST_CONTOUR_SIMPLIFICATION_UTIL_H__

#include <set>
using namespace std;

template<typename T>
void print_list(std::list<T*> &s, T* ioi, int to_show) {
  cerr << "|list<>| = " << s.size() << endl;

  bool found = false;
  typename list<T*>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	  found |= (*it == ioi);
	  if(!found)
		  continue;
	  if(to_show-- == 0)
		  return;
	  cerr << "| " << *(*it) << endl;
  }
}

template<typename T>
void print_set(std::set<T*> &s) {
  cerr << "|set<>| = " << s.size() << endl;
//  if(s.size() > 20)
//	  return;
  typename set<T*>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	cerr << " " << *(*it) << endl;
  }
}

template<typename T, typename D>
void print_set(std::set<T*,D> &s) {
  cerr << "|set<>| = " << s.size() << endl;
//  if(s.size() > 20)
//	  return;
  typename set<T*,D>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	cerr << " " << *(*it) << endl;
  }
}

template<typename T, typename D>
void print_set(std::set<T,D> &s) {
  cerr << "|set<>| = " << s.size() << endl;
//  if(s.size() > 20)
//	  return;
  typename set<T,D>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	cerr << " " << *it << endl;
  }
}

template<typename T>
void print_set(std::set<T> &s) {
  cerr << "|set<>| = " << s.size() << endl;
//  if(s.size() > 20)
//	  return;
  typename set<T>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	cerr << " " << (*it) << endl;
  }
}

template<typename T>
void print_set(std::set<T*,bool(*)(T*,T*)> &s) {
  cerr << "|set<>| = " << s.size() << endl;
//  if(s.size() > 20)
//	  return;
  typename set<T*,bool(*)(T*,T*)>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	cerr << " " << *(*it) << endl;
  }
}

template<typename T>
void print_set(std::set<T,bool(*)(T,T)> &s) {
  cerr << "|set<>| = " << s.size() << endl;
//  if(s.size() > 20)
//	  return;
  typename set<T,bool(*)(T,T)>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
	cerr << " " << *it << endl;
  }
}
	
template<typename T>
void print_set(std::set<T> &s, T item_of_interest) {
  cerr << "|set<>| = " << s.size() << endl;
  typename set<T>::iterator it_ioi = s.upper_bound(item_of_interest);
  typename set<T>::iterator it;
  int i = 0;
  if(it_ioi == s.end()) {
	  if(s.empty()) {
		  // i = 0;
	  }
	  else {
		  i = s.size()-1;
	  }
  }
  else {
	  item_of_interest = *it_ioi;
	  for (it = s.begin(); it != s.end(); ++it) {
		  if(!(*it < item_of_interest || item_of_interest < *it))
			  break;
		  i++;
	  }
  }

  int j = 0;
  for (it = s.begin(); it != s.end(); ++it) {
	  if(abs(j-i) < 10) 
		  cerr << " " << (*it) << endl;
	  j++;
  }
}

template<typename T>
void print_set(std::set<T*> &s, T *item_of_interest) {
	std::set<T> t;
	typename std::set<T*>::iterator it;
	for(it = s.begin(); it != s.end(); ++it) {
		t.insert(**it);
	}
	print_set(t, *item_of_interest);
}

#endif /*__TEST_CONTOUR_SIMPLIFICATION_UTIL_H__*/
