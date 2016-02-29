// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "intersect.h"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace tpie;
using namespace terrastream;

inline void update_map_info(map_info &m,triangle *t) {
	for (int i=0;i<3;i++) {
		m.minX = min(t->points[i].x,m.minX);
		m.maxX = max(t->points[i].x,m.maxX);
		m.minY = min(t->points[i].y,m.minY);
		m.maxY = max(t->points[i].y,m.maxY);
		m.minZ = min(t->points[i].z,m.minZ);
		m.maxZ = max(t->points[i].z,m.maxZ);
	}
}

map_info terrastream::intersect(stream<triangle> &in,elev_t gran,
								stream<signed_contour_segment> &out) {
	return intersect(in,gran,0.0f,out);
}

/*
  Intersect the triangle t at the elevation pz.
*/
inline void intersect_once(elev_t pz, elev_t minz, elev_t maxz,triangle* t,elev_t *zs,
						   stream<signed_contour_segment> &out) {
	if (pz < minz || pz > maxz) {
		return;
	}
	//Count number of vertices hit by the plane
	int hit_vertex[2];
	int hit_cnt=0;
	for (int i=0;i<3;i++) if (zs[i]==pz) hit_vertex[hit_cnt++]=i;
	if (hit_cnt==2) {
		//The contour lies on an edge.
		//Find the vertex not hit.
		int not_hit;
		for (int i=0;i<3;i++)
			if (i!=hit_vertex[0] && i!=hit_vertex[1]) not_hit=i;
		//Create the contour segment
		signed_contour_segment s(t->points[hit_vertex[0]].x,t->points[hit_vertex[0]].y,
								 t->points[hit_vertex[1]].x,t->points[hit_vertex[1]].y,
								 pz,t->points[not_hit].z<pz);
		out.write_item(s);
	} else if (hit_cnt==1) {
		if (pz==minz || pz==maxz) return; //Just a single point contour
		//The middle point got hit.
		//Find between which two points a and b, the triangle got hit as well.
		int a=0;
		int b;
		for (;a<3;a++) {
			b = (a+1)%3;
			if (zs[a]==pz || zs[b]==pz) continue;
			break;
		}
		//Compute hitting point by linear interpolation
		int min_p = a, max_p = b;
		if (t->points[max_p] < t->points[min_p])
			swap(max_p, min_p);
		double h = (double(pz-zs[min_p])/double(zs[max_p]-zs[min_p]));
		xycoord_t hx = t->points[min_p].x+h*(t->points[max_p].x - t->points[min_p].x);
		xycoord_t hy = t->points[min_p].y+h*(t->points[max_p].y - t->points[min_p].y);
		//Create the contour segment
		signed_contour_segment s(t->points[hit_vertex[0]].x,t->points[hit_vertex[0]].y,hx,hy,pz);
		out.write_item(s);
	} else {
		//No points got hit
		xycoord_t hxs[2];
		xycoord_t hys[2];
		int hits=0;
		//Check all edges of the triangle
		for (int a=0;a<3;a++) {
			int b = (a+1)%3;
			if ((pz-zs[a])*(pz-zs[b])>0)
				continue; //Both points are on same side of contour
			int min_p = a, max_p = b;
			if (t->points[max_p] < t->points[min_p])
				swap(max_p, min_p);
			double h = (double(pz-zs[min_p])/double(zs[max_p]-zs[min_p]));
			hxs[hits] = t->points[min_p].x+h*(t->points[max_p].x - t->points[min_p].x);
			hys[hits] = t->points[min_p].y+h*(t->points[max_p].y - t->points[min_p].y);
			hits++;
		}
		//Create the contour segment
		signed_contour_segment s(hxs[0],hys[0],hxs[1],hys[1],pz);
		out.write_item(s);
	}
}

map_info terrastream::intersect(stream<triangle> &in,elev_t gran,float z_diff,
								stream<signed_contour_segment> &out) {
	cerr << "Intersect: " << gran << "," << z_diff << endl;
	if (in.stream_len()==0) return map_info();
	in.seek(0);

	map_info res;
	triangle* t;
	in.read_item(&t);
	res.minX=res.maxX=t->points[0].x;
	res.minY=res.maxY=t->points[0].y;
	res.minZ=res.maxZ=t->points[0].z;
	tflow_progress progress("Intersecting triangles", "Intersecting triangles", 0, in.stream_len(), 1);
	do{
		update_map_info(res,t);
		//Start intersecting
		elev_t zs[3];
		for (int i=0;i<3;i++) zs[i]=t->points[i].z;
		elev_t minz = min(min(zs[0],zs[1]),zs[2]);
		elev_t maxz = max(max(zs[0],zs[1]),zs[2]);
		if (minz==maxz) continue; //Flat triangles gives no contours
		elev_t pz;
		int hs = int(ceil(minz/gran))-2;
		while ((pz=++hs*gran)<=maxz+gran) {
			if (z_diff > 0 && 2*z_diff <= gran) {
				if(2*z_diff < gran) {
					for (int i=-1;i<2;i++) {
						intersect_once(pz+i*z_diff,minz,maxz,t,zs,out);
					}
				}
				else {
					intersect_once(pz,minz,maxz,t,zs,out);
					intersect_once(pz-z_diff,minz,maxz,t,zs,out);
				}
			} 
			else {
				intersect_once(pz,minz,maxz,t,zs,out);
			}
		}
		progress.step();
	}while(in.read_item(&t)==ami::NO_ERROR);
	progress.done();
	return res;
}
