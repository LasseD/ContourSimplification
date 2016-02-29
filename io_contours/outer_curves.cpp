// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "outer_curves.h"
#include <terrastream/common/sort.h>

#include <iostream>

using namespace std;
using namespace tpie;
using namespace terrastream;

xycoord_t foo_sign(triangle_point p1,triangle_point p2,triangle_point p3) {
	return (p1.x-p3.x)*(p2.y-p3.y)-(p2.x-p3.x)*(p1.y-p3.y);
}

struct endpoint_segment{
	triangle_point p1,p2;
	endpoint_segment() {}
	endpoint_segment(triangle_point _p1,triangle_point _p2) {
		if (_p2.x<_p1.x || _p2.x==_p1.x && _p2.y<_p1.y) {
			swap(_p1,_p2);
		}
		p1=_p1; p2=_p2;
	}
	bool operator<(const endpoint_segment &s) const{
		if (p1.x!=s.p1.x) return p1.x<s.p1.x;
		if (p1.y!=s.p1.y) return p1.y<s.p1.y;
		if (p2.x==s.p2.x && p2.y==s.p2.y) return 0;
		//return 0;
		xycoord_t sgn = foo_sign(p1,p2,s.p2); //(p1.x-s.p2.x)*(p2.y-s.p2.y)-(p2.x-s.p2.x)*(p1.y-s.p2.y);
		//cout << "sign is " << sign << "\n";
		return sgn<0;
	}
};

struct cmp_segments{
	inline int compare(const endpoint_segment &i1, const endpoint_segment &i2) {
		if (i1.p1.x!=i2.p1.x) return i1.p1.x<i2.p1.x ? -1 : 1;
		if (i1.p1.y!=i2.p1.y) return i1.p1.y<i2.p1.y ? -1 : 1;
		if (i1.p2.x!=i2.p2.x) return i1.p2.x<i2.p2.x ? -1 : 1;
		if (i1.p2.y!=i2.p2.y) return i1.p2.y<i2.p2.y ? -1 : 1;
		return 0;
	}
};

void write_tangent(triangle_point p1,triangle_point p2,triangle_point p3,elev_t gran,map_info &info,stream<endpoint_segment> &out,bool upper) {
	xycoord_t dx1 = p2.x-p1.x;
	xycoord_t dx2 = p2.x-p3.x;
	xycoord_t dy1 = p2.y-p1.y;
	xycoord_t dy2 = p2.y-p3.y;
	xycoord_t len1 = static_cast<xycoord_t>(sqrt(dx1*dx1+dy1*dy1));
	xycoord_t len2 = static_cast<xycoord_t>(sqrt(dx2*dx2+dy2*dy2));
	dx1/=len1; dy1/=len1;
	dx2/=len2; dy2/=len2;
	xycoord_t dir_x = dx1+dx2;
	xycoord_t dir_y = dy1+dy2;
	if (dir_x==0 && dir_y==0) {
		//Special case, colinear
		if (upper) {
			dir_x = -dy1;
			dir_y = dx1;
		} else {
			dir_x = dy1;
			dir_y = -dx1;
		}
	}
	const elev_t slope = (info.minZ-info.maxZ);
	xycoord_t tan_len = static_cast<xycoord_t>(((info.minZ-gran)-p2.z)/slope)*4;
	xycoord_t dir_len = static_cast<xycoord_t>(sqrt(dir_x*dir_x+dir_y*dir_y));
	dir_x *= tan_len/dir_len;
	dir_y *= tan_len/dir_len;
	triangle_point p4(p2.x+dir_x,p2.y+dir_y,info.minZ-gran);
	endpoint_segment tangent(p2,p4);
	//cout << "Tan: " << p2.x << "," << p2.y << "->" << p4.x << "," << p4.y << "\n";
	out.write_item(tangent);
}

void intersect(endpoint_segment *cur,endpoint_segment *next,elev_t gran,stream<labelling_signed_contour_segment> &out) {
	xycoord_t minz_cur = min(cur->p1.z,cur->p2.z);
	xycoord_t maxz_cur = max(cur->p1.z,cur->p2.z);
	xycoord_t minz_next = min(next->p1.z,next->p2.z);
	xycoord_t maxz_next = max(next->p1.z,next->p2.z);
	elev_t minz = min(minz_cur,minz_next);
	elev_t maxz = max(maxz_cur,maxz_next);
	elev_t pz;
	int hs = int(ceil(minz/gran))-1;
	while ((pz=++hs*gran)<maxz) {
		if (pz < minz)
			continue;
		//Linear interpolation
		if (pz>min(maxz_cur,maxz_next)) {
			//Intersection is between tangent and the two boundary vertices
			xycoord_t cur_boundX = cur->p1.z<cur->p2.z ? cur->p2.x : cur->p1.x;
			xycoord_t cur_boundY = cur->p1.z<cur->p2.z ? cur->p2.y : cur->p1.y;
			xycoord_t next_boundX = next->p1.z<next->p2.z ? next->p2.x : next->p1.x;
			xycoord_t next_boundY = next->p1.z<next->p2.z ? next->p2.y : next->p1.y;
			xycoord_t zcur = maxz_cur, znext = maxz_next;
			if (next_boundX < cur_boundX || (next_boundX == cur_boundX && next_boundY < cur_boundY)) {
				swap(zcur, znext);
				swap(cur_boundX, next_boundX);
				swap(cur_boundY, next_boundY);
			}
			xycoord_t t1 = (pz-zcur)/(znext-zcur);
			xycoord_t nx1 = cur_boundX+t1*(next_boundX-cur_boundX);
			xycoord_t ny1 = cur_boundY+t1*(next_boundY-cur_boundY);
			endpoint_segment *tang = maxz_cur>maxz_next ? cur : next;
			xycoord_t t2 = (pz-tang->p1.z)/(tang->p2.z-tang->p1.z);
			xycoord_t nx2 = tang->p1.x + t2*(tang->p2.x-tang->p1.x);
			xycoord_t ny2 = tang->p1.y + t2*(tang->p2.y-tang->p1.y);
			labelling_signed_contour_segment s(nx1,ny1,nx2,ny2,pz,true);
			out.write_item(s);
		} else {
			//Intersection is on the tangents
			xycoord_t t1 = (pz-cur->p1.z)/(cur->p2.z-cur->p1.z);
			xycoord_t t2 = (pz-next->p1.z)/(next->p2.z-next->p1.z);
			xycoord_t nx1 = cur->p1.x+t1*(cur->p2.x-cur->p1.x);
			xycoord_t ny1 = cur->p1.y+t1*(cur->p2.y-cur->p1.y);
			xycoord_t nx2 = next->p1.x+t2*(next->p2.x-next->p1.x);
			xycoord_t ny2 = next->p1.y+t2*(next->p2.y-next->p1.y);
			labelling_signed_contour_segment s(nx1,ny1,nx2,ny2,pz,true);
			out.write_item(s);
		}
	}
}

void terrastream::add_outer_curves(stream<triangle> &tris,elev_t gran,map_info &info,stream<labelling_signed_contour_segment> &out) {
	if (tris.stream_len()==0) return;
	tris.seek(0);
	//Create a stream of the edges in the triangulation
	tflow_progress create_progress("Creating edge stream from triangles", "Creating edge stream from triangles", 0, tris.stream_len(), 1);
	stream<endpoint_segment> edges;
	triangle *t;
	while (tris.read_item(&t)==NO_ERROR) {
		create_progress.step();
		for (int i=0;i<3;i++) {
			endpoint_segment s(t->points[i],t->points[(i+1)%3]);
			edges.write_item(s);
		}
	}
	create_progress.done();
	//Sort the edges by left endpoint
	tflow_progress sort_progress("Sorting segments", "Sorting segments", 0, edges.stream_len(), 1);
	cmp_segments cmp;
	edges.seek(0);
	ts_sort(&edges, &cmp,&sort_progress);
	edges.seek(0);
	//Boundary edges have no duplicates, others have
	//Sweep left to right, maintaining upper and lower hull intersection with sweep line.
	//Write tangents to stream
	stream<endpoint_segment> upper_tan;
	stream<endpoint_segment> lower_tan;
	endpoint_segment cur, *next;
	triangle_point up,down,up_prev,down_prev;
	//Find topmost segment
	edges.read_item(&next);
	cur=*next;
	endpoint_segment *tmp_edge, first_edge;
	while ((edges.read_item(&tmp_edge))==NO_ERROR && tmp_edge->p1.x == cur.p1.x && tmp_edge->p1.y == cur.p1.y)
		if (*tmp_edge < cur) cur=*tmp_edge;
	first_edge=cur;
	edges.seek(0);
	down_prev = up = cur.p2;
	up_prev = down = cur.p1;
	//cout << "up is " << up.x << " " << up.y << " " << up.z << "\n";
	//cout << "down is " << down.x << " " << down.y << " " << down.z << "\n";
	tflow_progress generate_progress("Generating tangents", "Generating tangents", 0, edges.stream_len(), 1);
	while (true) {
		bool is_internal=false;
		ami::err er;
		while ((er=edges.read_item(&next))==NO_ERROR && cmp.compare(cur, *next) == 0) {
			generate_progress.step();
			is_internal=true;
		}
		if (!is_internal && cmp.compare(cur,first_edge)!=0) {
			//We have a boundary segment
			if (cur.p1.x==up.x && cur.p1.y==up.y) {
				//Upper hull
				triangle_point up_new = cur.p2;
				write_tangent(up_prev,up,up_new,gran,info,upper_tan,true);
				//Prepare for next
				up_prev=up;
				up=up_new;
				//cout << "up is " << up.x << " " << up.y << " " << up.z << "\n";
			} else {
				//Lower hull
				triangle_point down_new = cur.p2;
				write_tangent(down_prev,down,down_new,gran,info,lower_tan,false);
				//Prepare for next
				down_prev=down;
				down=down_new;
				//cout << "down is " << down.x << " " << down.y << " " << down.z << " and prev " << down_prev.x << " " << down_prev.y << "\n";
			}
		}
		if (er!=NO_ERROR) {
			//Reached end of stream
			break;
		}
		//Prepare for next iteration
		cur=*next;
	}
	//Write the last tangent line
	write_tangent(up_prev,up,down_prev,gran,info,upper_tan,true);
	generate_progress.done();
	//Output segments
	upper_tan.seek(0);
	lower_tan.seek(0);
	//Begin with upper hull segments
	tflow_progress outer_progress("Creating outer contours", "Creating outer contours", 0, upper_tan.stream_len() + lower_tan.stream_len(), 1);
	upper_tan.read_item(&next);
	cur = *next;
	endpoint_segment first_upper = cur;
	while (upper_tan.read_item(&next)==NO_ERROR) {
		outer_progress.step();
		//Intersect
		intersect(&cur,next,gran,out);
		cur=*next;
	}
	endpoint_segment last_upper = cur;
	//Do lower hull segments
	lower_tan.read_item(&next);
	cur = *next;
	endpoint_segment first_lower = cur;
	while (lower_tan.read_item(&next)==NO_ERROR) {
		outer_progress.step();
		intersect(&cur,next,gran,out);
		cur=*next;
	}
	endpoint_segment last_lower = cur;
	outer_progress.done();
	//Do last intersections
	intersect(&first_upper,&first_lower,gran,out);
	intersect(&last_upper,&last_lower,gran,out);
}
