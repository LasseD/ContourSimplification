// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
//
// Copyright 2010, 2011 SCALGO development

#include "contour_to_shape.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include "io_contours/contour_types.h"
#include <stdlib.h>
#include <string.h>
#include <shapefil.h>

#include <terrastream/common/common.h>
#include <tpie/portability.h>
#include <tpie/stream.h>
#include <set>
#include <terrastream/common/tflow_types.h>
#include <terrastream/common/wlabel.h>
#include "decomposition.h"

using namespace std;
using namespace terrastream;
//using namespace tpie;
using namespace tpie::ami;
using namespace simplification;

/*
  SHPT_POLYGON (2D polygon without measure)
 */
static int z_field_index;
static int id_field_index;
static int contour_index = -1;
static int contour_index_o = -1;
static int contour_index_b = -1;
static int contour_index_s = -1;
bool outputPoints(SHPHandle &shpHandle, DBFHandle &dbfHandle,
				  SHPHandle &shpHandle_o, DBFHandle &dbfHandle_o,
				  SHPHandle &shpHandle_b, DBFHandle &dbfHandle_b,
				  SHPHandle &shpHandle_s, DBFHandle &dbfHandle_s,
				  vector<contour_point> &points, elev_t z,
				  float contour_interval, float e_z,
				  bool is_simplified) {
	bool is_level = is_level_line(z, contour_interval, e_z);
	int npoints = points.size();

	if (npoints < 3) {
		return false;
	}
//	cerr << "Outputting " << points.size() << " points for " << points[0].label << " at " << z << endl;
	double *x = new double[npoints];
	double *y = new double[npoints];
  
	int i = 0;
	for(vector<contour_point>::iterator it = points.begin(); it != points.end(); ++it) {
		x[i] = it->x;
		y[i] = -it->y;
//		cerr << it->label << ":" << x[i] << "," << y[i] << endl;
		i++;
	}

	int id = points[0].label;
	SHPObject *o = SHPCreateObject(SHPT_POLYGON, ++contour_index, 0, NULL, NULL,
								   npoints, x, y, NULL, NULL);
	SHPWriteObject(shpHandle, -1, o);
	SHPDestroyObject(o);
//  cerr << "Write att. " << contour_index << "," << id_field_index << ":" << id << endl;
	int ok = DBFWriteIntegerAttribute(dbfHandle, contour_index, id_field_index, id);
	assert(ok);
	ok = DBFWriteDoubleAttribute(dbfHandle, contour_index, z_field_index, z);
	assert(ok);
	
	if (is_level && is_simplified) {
		o = SHPCreateObject(SHPT_POLYGON, ++contour_index_s, 0, NULL, NULL,
							npoints, x, y, NULL, NULL);
		SHPWriteObject(shpHandle_s, -1, o);
		SHPDestroyObject(o);
		ok = DBFWriteIntegerAttribute(dbfHandle_s, contour_index_s, id_field_index, id);
		assert(ok);
		ok = DBFWriteDoubleAttribute(dbfHandle_s, contour_index_s, z_field_index, z);
		assert(ok);
	} else if (!is_level) {
		if(is_simplified) {
			cerr << "Error: contour " << contour_index << ", id counter " << id_field_index << ":" << id << " at z:" << z << ", gran:" << contour_interval << ", e_z:" << e_z << endl;
		}
		assert(!is_simplified);
		o = SHPCreateObject(SHPT_POLYGON, ++contour_index_b, 0, NULL, NULL,
							npoints, x, y, NULL, NULL);
		SHPWriteObject(shpHandle_b, -1, o);
		SHPDestroyObject(o);
		ok = DBFWriteIntegerAttribute(dbfHandle_b, contour_index_b, id_field_index, id);
		assert(ok);
		ok = DBFWriteDoubleAttribute(dbfHandle_b, contour_index_b, z_field_index, z);
		assert(ok);
	} else {
		o = SHPCreateObject(SHPT_POLYGON, ++contour_index_o, 0, NULL, NULL,
							npoints, x, y, NULL, NULL);
		SHPWriteObject(shpHandle_o, -1, o);
		SHPDestroyObject(o);
		ok = DBFWriteIntegerAttribute(dbfHandle_o, contour_index_o, id_field_index, id);
		assert(ok);
		ok = DBFWriteDoubleAttribute(dbfHandle_o, contour_index_o, z_field_index, z);
		assert(ok);
	}

	delete[] x;
	delete[] y;
	return true;
}


void outputStreams(SHPHandle &shpHandle, DBFHandle &dbfHandle,
				  SHPHandle &shpHandle_o, DBFHandle &dbfHandle_o,
				  SHPHandle &shpHandle_b, DBFHandle &dbfHandle_b,
				  SHPHandle &shpHandle_s, DBFHandle &dbfHandle_s,
				  stream<contour_point> &uns,
				  stream<contour_point> &sim,
				  stream<topology_edge> &topology,
				  float contour_interval, float e_z) {
	contour_point *pu, *ps;
	topology_edge *topo;
	bool u_ok = uns.read_item(&pu) == NO_ERROR;
	bool s_ok = sim.read_item(&ps) == NO_ERROR;
	topology.seek(0);
	while(topology.read_item(&topo) == NO_ERROR) {
		if(u_ok && pu->label == topo->c) {
			vector<contour_point> pts;
			do {
				pts.push_back(*pu);			
			}
			while((u_ok = (uns.read_item(&pu) == NO_ERROR)) && pu->label == topo->c);
			outputPoints(shpHandle, dbfHandle, shpHandle_o, dbfHandle_o,
						 shpHandle_b, dbfHandle_b, shpHandle_s, dbfHandle_s,
						 pts, topo->c_z,
						 contour_interval, e_z, false);
		}
		if(s_ok && ps->label == topo->c) {
			vector<contour_point> pts;
			do {
				pts.push_back(*ps);			
			}
			while((s_ok = (sim.read_item(&ps) == NO_ERROR)) && ps->label == topo->c);
			outputPoints(shpHandle, dbfHandle, shpHandle_o, dbfHandle_o,
						 shpHandle_b, dbfHandle_b, shpHandle_s, dbfHandle_s,
						 pts, topo->c_z,
						 contour_interval, e_z, true);
		}
	}
}

void shape::to_shape(char const* const file_suffix,
					 stream<contour_point>& unsimplified_stream,
					 stream<topology_edge> &topology,				   
					 stream<contour_point>* const simplified_stream,
					 float contour_interval, float e_z) {
	SHPHandle shpHandle = SHPCreate(file_suffix, SHPT_POLYGON);
	string s = ((string)file_suffix) + "_o";
	SHPHandle shpHandle_o = SHPCreate(s.c_str(), SHPT_POLYGON);
	s = ((string)file_suffix) + "_b";
	SHPHandle shpHandle_b = SHPCreate(s.c_str(), SHPT_POLYGON);
	s = ((string)file_suffix) + "_s";
	SHPHandle shpHandle_s = SHPCreate(s.c_str(), SHPT_POLYGON);

	s = ((string)file_suffix) + ".dbf";
	DBFHandle dbfHandle = DBFCreate(s.c_str());
	s = ((string)file_suffix) + "_o.dbf";
	DBFHandle dbfHandle_o = DBFCreate(s.c_str());
	s = ((string)file_suffix) + "_b.dbf";
	DBFHandle dbfHandle_b = DBFCreate(s.c_str());
	s = ((string)file_suffix) + "_s.dbf";
	DBFHandle dbfHandle_s = DBFCreate(s.c_str());
	
	string z_name = "z";
	z_field_index = DBFAddField(dbfHandle, z_name.c_str(), FTDouble, 8, 4);
	DBFAddField(dbfHandle_o, z_name.c_str(), FTDouble, 8, 4);
	DBFAddField(dbfHandle_b, z_name.c_str(), FTDouble, 8, 4);
	DBFAddField(dbfHandle_s, z_name.c_str(), FTDouble, 8, 4);
//  cerr << " added z field on " << z_field_index << endl;
	string id_name = "label";
	id_field_index = DBFAddField(dbfHandle, id_name.c_str(), FTInteger, 8, 0);
	DBFAddField(dbfHandle, id_name.c_str(), FTInteger, 8, 0);
	DBFAddField(dbfHandle_o, id_name.c_str(), FTInteger, 8, 0);
	DBFAddField(dbfHandle_b, id_name.c_str(), FTInteger, 8, 0);
	DBFAddField(dbfHandle_s, id_name.c_str(), FTInteger, 8, 0);
//  cerr << " added label field on " << id_field_index << endl;
	

	outputStreams(shpHandle, dbfHandle,
				  shpHandle_o, dbfHandle_o,
				  shpHandle_b, dbfHandle_b,
				  shpHandle_s, dbfHandle_s,
				  unsimplified_stream, *simplified_stream, topology, contour_interval, e_z);
	
	if (simplified_stream != NULL) {
		simplified_stream->seek(0);
	}
	unsimplified_stream.seek(0);
	topology.seek(0);
	
	DBFClose(dbfHandle);
	DBFClose(dbfHandle_o);
	DBFClose(dbfHandle_b);
	DBFClose(dbfHandle_s);
	SHPClose(shpHandle);
	SHPClose(shpHandle_o);
	SHPClose(shpHandle_b);
	SHPClose(shpHandle_s);
}
