// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "cw_ordering.h"
#include <terrastream/common/sort.h>
#include <vector>
#include <queue>
#include <map>
#include <utility>

using namespace terrastream;
using namespace tpie;
using namespace std;

//#define DEBUG_CW_ORDERING

typedef labelled_signed_contour_segment lss;
typedef ranked_labelled_signed_contour_segment rlss;
typedef pair<xycoord_t,xycoord_t> pt;
typedef unsigned int CW_SIZE_T; // The size_t type used throughout the code -- kept as 32-bit to avoid doubling memory usage

/* Order segments by label, and secondarely such that the first segment from each contour
 * has the leftmost endpoint and lies on the upper hull.
 */
struct order_for_augment{
  int compare(const lss &l1,const lss &l2) const {
	if (l1.label!=l2.label) return l1.label-l2.label;
	if (l1.x1!=l2.x1) return l1.x1<l2.x1 ? -1 : 1;
	if (l1.y1!=l2.y1) return l1.y1<l2.y1 ? -1 : 1;
	//(p1.x-p3.x)*(p2.y-p3.y)-(p2.x-p3.x)*(p1.y-p3.y);
	xycoord_t turn_sign = (l1.x1-l2.x2)*(l1.y2-l2.y2)-(l1.x2-l2.x2)*(l1.y1-l2.y2);
	return turn_sign ==0 ? 0 : (turn_sign<0 ? -1 : 1);
  }
};

struct ep{
	pt p;
	vector<pair<CW_SIZE_T,CW_SIZE_T> > ns; // adjacent point, seg index in contour
	ep(pt _p) : p(_p) {}
};

/*void write_contour(vector<lss> &contour,vector<CW_SIZE_T> &segs, int label, stream<rlss> &out) {
  order_for_augment cmp;
  //Check for real segments and start the contours at the leftmost segment
  bool has_real=false;
  CW_SIZE_T first=0;
  for (CW_SIZE_T i=0;i<segs.size();i++) {
	has_real |= contour[segs[i]].sign==0;
	if (cmp.compare(contour[segs[i]],contour[segs[first]])<0)
	  first=i;
  }
  if (!has_real) return;
  //Output the contour
  //cerr << "Outputting contour " << label << std::endl;
  int ranks = 0;
  for (CW_SIZE_T i=first;i<segs.size();i++) {
	lss &l = contour[segs[i]];
	out.write_item(rlss(l.x1,l.y1,l.x2,l.y2,l.z,l.sign,label,ranks++));
//	cerr << " " << l << std::endl;
  }
  for (CW_SIZE_T i=0;i<first;i++) {
	lss &l = contour[segs[i]];
	out.write_item(rlss(l.x1,l.y1,l.x2,l.y2,l.z,l.sign,label,ranks++));
//	cerr << " " << l << std::endl;
  }
  }*/

const xycoord_t PI = acos(-1.0);

void build_eps(vector<ep> &eps, vector<lss> &contour) {
	map<pt,CW_SIZE_T> ep_map; // pt: pair<xycoord_t,xycoord_t>, ep_map: point->index
	for (CW_SIZE_T i=0;i<contour.size();i++) {
		lss &l = contour[i];
		pt p1 = pt(l.x1,l.y1);
		pt p2 = pt(l.x2,l.y2);
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
}

//Sort neighbors of nodes of degree > 2 (in clockwise order)
void order_points(vector<ep> &eps) {
	for (CW_SIZE_T i=0;i<eps.size();i++) {
		if (eps[i].ns.size()>2) {
			vector<pair<xycoord_t,CW_SIZE_T> > tmp; // cos, index.
			for (CW_SIZE_T j=0;j<eps[i].ns.size();j++) {
				pair<CW_SIZE_T,CW_SIZE_T> n = eps[i].ns[j];
				xycoord_t dx = eps[n.first].p.first-eps[i].p.first;
				xycoord_t dy = eps[n.first].p.second-eps[i].p.second;
				xycoord_t l = sqrt(dx*dx+dy*dy);
				dx/=l;
				xycoord_t acos_r = acos(dx);
				if (acos_r<0) acos_r += 2*PI;
				if (dy>0 && acos_r>PI || dy<0 && acos_r<PI)
					acos_r = 2*PI-acos_r;
				//For cw order, we negate the angle before sorting
				acos_r = -acos_r;
				tmp.push_back(make_pair(acos_r,j));
			}
			sort(tmp.begin(),tmp.end());
			vector<pair<CW_SIZE_T,CW_SIZE_T> > new_ns;
			for (CW_SIZE_T j=0;j<tmp.size();j++) {
				new_ns.push_back(eps[i].ns[tmp[j].second]);
			}
			eps[i].ns=new_ns; // TODO: Actually updated?
		}
	}
}

/*void order_contour(vector<lss> &contour,int &label, stream<rlss> &out) {
#ifdef DEBUG_CW_ORDERING
	if(!contour.empty()) {
		cerr << "Ordering contour at " << contour[0].z << " of size " << contour.size() << endl;
		for(vector<lss>::iterator it = contour.begin(); it != contour.end(); ++it) {
			cerr << " " << *it << endl;
		}
	}
#endif
	//Build the graph
	vector<ep> eps; // ep: p=pt point, ns=vector(unref neighbour index,seg index in contour)
	build_eps(eps, contour);

	//Sort neighbors of nodes of degree > 2 (in clockwise order)
	order_points(eps);

	//Traverse the contour
	stack<CW_SIZE_T> multi_fanout;
	vector<CW_SIZE_T> segs;
	vector<CW_SIZE_T> pts;
	segs.push_back(0);
	pts.push_back(0);
	pts.push_back(1);
	CW_SIZE_T at = 1; // Where we are in pts
	CW_SIZE_T from = 0; // where we came from in pts
	vector<bool> vis(contour.size(),false); // visited segments
	vector<bool> pt_vis(eps.size(),false); // visited points
	vis[0]=true;
	pt_vis[0]=true;

	while (true) {
		//if(pt_vis[at]) {
		if (at==pts[0]) {
			//We finished the current contour-cycle
			write_contour(contour,segs,label++,out);
			segs.clear();
			pts.clear();
			//Backtract to point with multi_fanout // I smell hack... / LD
			while (true) {
				//If none exists, we are done
				if (multi_fanout.empty()) 
					return;
				//Start a new cycle from here, if possible
				CW_SIZE_T back_to = multi_fanout.top();
				ep &e = eps[back_to];
#ifdef DEBUG_CW_ORDERING
				cerr << "Handle multiple fanout point (" << e.p.first << "," << e.p.second << ") at " << back_to << endl;
#endif
				CW_SIZE_T j=0;
				for (;j<e.ns.size();j++)
					if (vis[e.ns[j].second]) 
						break;
				CW_SIZE_T k;
				bool found=false;
				for (CW_SIZE_T i=0;i<e.ns.size();i++) {
					k = (j+i)%e.ns.size();
					if (!vis[e.ns[k].second]) {
						//Found next edge to leave through
						found=true;
						break;
					}
				}
				if (!found) {
					//This multifanout point has no exits left
					//Look for another multifanout point to continue from
#ifdef DEBUG_CW_ORDERING
					cerr << "Multiple fanout without exits left:" << back_to << endl;
#endif
					multi_fanout.pop();
					continue;
				}
				//We found an edge to exit from
				from = back_to;
				at = e.ns[k].first;
				segs.push_back(e.ns[k].second);
				pts.push_back(from);
				pts.push_back(at);
				vis[e.ns[k].second]=true;
				break;
			}
		} // if at 0

		pt_vis[at]=true;
		//Continue the traversal
		CW_SIZE_T next;
		if (eps[at].ns.size()!=2) {
			//Find next neighbor to exit through
			CW_SIZE_T cur;
			for (CW_SIZE_T i=0;i<eps[at].ns.size();i++)
				if (eps[at].ns[i].first==from) {
					cur = i;
					break;
				}
			//Attempt to end the current loop
			next = (cur+1+eps[at].ns.size())%eps[at].ns.size(); // MUST FAIL! - -> +
			//Remember that we saw a multifanout point
			multi_fanout.push(at);
		} else {
			//This is a normal node with two neighbors
			next = vis[eps[at].ns[0].second] ? 1 : 0;
		}
		pair<CW_SIZE_T,CW_SIZE_T> des = eps[at].ns[next];
		vis[des.second]=true;
		segs.push_back(des.second);
		pts.push_back(des.first);
		from=at;
		at=des.first;
	}
	}*/

void write_contour(vector<lss> &contour, int &label, stream<rlss> &out) {
	//Check for real segments and start the contours at the leftmost segment
	bool has_real=false;
	int first=0;
	for(CW_SIZE_T i=0;i<contour.size();i++) {
		has_real |= contour[i].sign==0;
		if(contour[i].x1 < contour[first].x1)
			first=i;
	}
#ifdef DEBUG_CW_ORDERING
	cerr << " Leftmost segment of contour " << label << ": " << contour[first] << endl;
#endif

	if (!has_real) {
#ifdef DEBUG_CW_ORDERING
		cerr << "Contour " << label << " without real segment. Not output." << endl;
#endif
		return;
	}
	
	lss first_segment = contour[first];
	lss prev_segment = contour[(first-1+contour.size())%contour.size()];
#ifdef DEBUG_CW_ORDERING
	cerr << "prev seg: " << prev_segment << endl;
	cerr << "first seg: " << first_segment << endl;
#endif
	assert(prev_segment.x2 == first_segment.x1 && prev_segment.y2 == first_segment.y1);
	bool reverse = false;
	if(first_segment.x1 == first_segment.x2) {
		reverse = first_segment.y1 > first_segment.y2;
	}
	else if(prev_segment.x1 == prev_segment.x2) {
		reverse = prev_segment.y1 > prev_segment.y2;
	}
	else {
		lss l2 = first_segment;
		lss l1 = prev_segment;
		xycoord_t turn_sign = (l1.x1-l2.x2)*(l1.y2-l2.y2)-(l1.x2-l2.x2)*(l1.y1-l2.y2);
		reverse = turn_sign > 0;
	}

	//Output the contour
	if(reverse) {
#ifdef DEBUG_CW_ORDERING
		cerr << "Outputting reversed contour " << label << " of size " << contour.size() << std::endl;
#endif
		int ranks = 0;
		for (int i=first;i>=0;--i) {
			lss &l = contour[i];
			out.write_item(rlss(l.x2,l.y2,l.x1,l.y1,l.z,l.sign,label,ranks++));
#ifdef DEBUG_CW_ORDERING
			cerr << " " << l << endl;
#endif
		}
		for (int i=contour.size()-1;i>first;--i) {
			lss &l = contour[i];
			out.write_item(rlss(l.x2,l.y2,l.x1,l.y1,l.z,l.sign,label,ranks++));
			cerr << " " << l << std::endl;
		}
	}
	else {
#ifdef DEBUG_CW_ORDERING
		cerr << "Outputting non-reversed contour " << label << " of size " << contour.size() << std::endl;
#endif
		int ranks = 0;
		for (CW_SIZE_T i=first;i<contour.size();i++) {
			lss &l = contour[i];
#ifdef DEBUG_CW_ORDERING
			cerr << " " << l << std::endl;
#endif
			out.write_item(rlss(l.x1,l.y1,l.x2,l.y2,l.z,l.sign,label,ranks++));
		}
		for (CW_SIZE_T i=0;i<first;i++) {
			lss &l = contour[i];
			out.write_item(rlss(l.x1,l.y1,l.x2,l.y2,l.z,l.sign,label,ranks++));
#ifdef DEBUG_CW_ORDERING
			cerr << " " << l << std::endl;
#endif
		}
	}

	label++;
}

void order_contour2(vector<lss> &contour,int &label, stream<rlss> &out) {
#ifdef DEBUG_CW_ORDERING
	if(!contour.empty()) {
		cerr << "Ordering contour at " << contour[0].z << " of size " << contour.size() << endl;
		for(vector<lss>::iterator it = contour.begin(); it != contour.end(); ++it) {
			cerr << " " << *it << endl;
		}
	}
#endif
	//Build the graph
	vector<ep> eps; // ep: p=pt point, ns=vector(unref neighbour index,seg index in contour)
	build_eps(eps, contour);

	//Sort neighbors of nodes of degree > 2 (in clockwise order)
	order_points(eps);

	// Graph needed for simple alg.:
	map<pt,vector<CW_SIZE_T> > neighbour_info; // pt -> neighbours
	for(vector<ep>::iterator it = eps.begin(); it != eps.end(); ++it) {
		ep point_info = *it;
		vector<CW_SIZE_T> ns2;
		vector<pair<CW_SIZE_T,CW_SIZE_T> >::iterator it2;
		for(it2 = point_info.ns.begin(); it2 != point_info.ns.end(); ++it2) {
			ns2.push_back(it2->second);
		}
		if(ns2.size() % 2 == 1) {
			cerr << "Error: Input error to cw_ordering. Contour " << label << " with node of degree " << ns2.size() << endl;						
			cerr << "Warning: Contour removed from input. " << endl;			
			return;
		}
//		assert(ns2.size() % 2 == 0);
		neighbour_info[point_info.p] = ns2;
	}
	eps.clear(); // clear space.

#ifdef DEBUG_CW_ORDERING
	cerr << "Data structure neighbour_info built ||=" << neighbour_info.size() << endl;
	for(map<pt,vector<CW_SIZE_T> >::iterator it = neighbour_info.begin(); it != neighbour_info.end(); ++it) {
			if(it->second.size() == 2)
				continue;
			cerr << " " << it->first.first << "," << it->first.second << ": ";
			for(vector<CW_SIZE_T>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
				cerr << "," << *it2;
			}
			cerr << endl;
	}
#endif

	// Segment state:
	vector<bool> seg_vis(contour.size(),false); // visited segs
	stack<CW_SIZE_T> multi_fanout;

	int to_handle = 0; 
	vector<lss> new_contour;
	while(to_handle < contour.size()) {
		int contour_index = to_handle;
		bool current_oriented = true; // from p1 to p2.
#ifdef DEBUG_CW_ORDERING
		cerr << "Starting to output new contour from " << contour_index << " (" << contour[contour_index] << ")" << endl;
#endif
		// handle contour from "to_handle" by walking until any orig encountered:
		do {
			assert(!seg_vis[contour_index]);
			// "output" seg:
			lss seg(contour[contour_index]);
			if(!current_oriented) {
				swap(seg.x1,seg.x2);
				swap(seg.y1,seg.y2);
			}
			new_contour.push_back(seg);
			seg_vis[contour_index] = true;
#ifdef DEBUG_CW_ORDERING
			cerr << "Adding " << seg << endl;
#endif
			// find next:
			CW_SIZE_T next; // next seg.
			pt p2 = make_pair(seg.x2,seg.y2);
			vector<CW_SIZE_T> p2n = neighbour_info[p2];
			if (p2n.size()!=2) {
				//Find self
				CW_SIZE_T cur;
				for(CW_SIZE_T i = 0; i < p2n.size(); i++) {
					if(p2n[i] == contour_index) {
						cur = i;
						break;
					}
				}
				next = p2n[(cur+1)%p2n.size()];
				// Update multiple fanout node:
				vector<CW_SIZE_T> np2n;
				for(CW_SIZE_T i = 0; i < p2n.size(); i++) {
					if(!(i == cur || i == (cur+1)%p2n.size())) {
						np2n.push_back(p2n[i]);
					}
				}
				assert(np2n.size() == p2n.size()-2);
				neighbour_info[p2] = np2n;
				assert(neighbour_info[p2].size() == np2n.size());
#ifdef DEBUG_CW_ORDERING
				cerr << "Point " << p2.first << "," << p2.second << " now with neighbours: ";
				for(vector<CW_SIZE_T>::iterator it3 = np2n.begin(); it3 != np2n.end(); ++it3) {
					cerr << "," << *it3;
				}
				cerr << endl;
#endif
			} else { //This is a normal node with two neighbors
				next = p2n[0] == contour_index ? p2n[1] : p2n[0];				
			}
#ifdef DEBUG_CW_ORDERING
			cerr << "Going to " << next << " (" << contour[next] << ")" << endl;
#endif
			// update orientation:
			current_oriented = make_pair(contour[next].x1,contour[next].y1) == p2;
			contour_index = next;
		}
		while(contour_index != to_handle);

		// Output new_contour
#ifdef DEBUG_CW_ORDERING
		cerr << " Writing contour for contour " << contour[0].label << endl;
#endif
		write_contour(new_contour, label, out);
		new_contour.clear();

		// update to_handle:
		while(seg_vis[to_handle]) {
			++to_handle;
			if(to_handle >= contour.size())
				return;
		}
	}
}

void terrastream::print_stream_status(stream<lss> &in) {
  cout << "Stream status:" << endl;

  lss* t;
  err ae;
  int label = -1;
  while (true) {
	  int contour_size = 1;
	  while ((ae=in.read_item(&t))==NO_ERROR && t->label==label) {
		  contour_size++;
	  }
	  if(label != -1)
		  cerr << "contour " << label << " of size " << contour_size << endl;
	  label = t->label;
	  if (ae!=NO_ERROR) break;
  }
  in.seek(0);
}

void terrastream::cw_order(stream<lss> &in,stream<rlss> &out) {
  if (in.stream_len()==0) return;
  cout << "Sorting contours into clockwise order\n";
  cout << "Assuming a single contour fits in memory\n";

  //Sort
  order_for_augment cmp;
  in.seek(0);
  ts_sort(&in,&cmp);
  in.seek(0);
//  cout << "Done sorting.\n";

  //print_stream_status(in);

  //Retrieve single contours
  vector<lss> contour;
  lss* t;
  err ae = in.read_item(&t);
  int label = 0;
  tflow_progress progress("Performing clockwise ordering of segments", "Performing clockwise ordering of segments", 0, in.stream_len(), 1);
  while(ae==NO_ERROR) {
	//Extract the next contour and check that it contains a real segment
	contour.clear();
	contour.push_back(*t);

	int contour_size = 1;
	while ((ae=in.read_item(&t))==NO_ERROR && t->label==contour[0].label) {
	  progress.step();
	  contour.push_back(*t);
	  contour_size++;
	}
	
	if (contour.size() > std::numeric_limits<CW_SIZE_T>::max())
		throw std::runtime_error("Too many segments in contour, redefine CW_SIZE_T to larger type and recompile");
	//Process it if it has a real segment
	order_contour2(contour,label,out);
	
	//Advance to next contour...
  }
  progress.done();
  cout << "Done cw-ordering.\n";
}
