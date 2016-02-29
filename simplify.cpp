// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
//
// Copyright 2010 SCALGO development

#include <stdio.h>
#include "io_contours/contour_types.h"
#include "io_contours/mif_outputter.h"
#include <tpie/portability.h>
#include <cstdlib>
#include <tpie/stream.h>
#include "contour_reader.h"
#include "contour_to_shape.h"
#include "contour_simplification.h"
#include <tpie/persist.h>
#include <set>
#include <time.h>

using namespace tpie;
using namespace ami;

typedef triangulated_ranked_labelled_signed_contour_segment ts;
typedef ranked_labelled_signed_contour_segment rlss;
typedef elev_t height_type;
typedef contour_point cp;

namespace terrastream {
namespace simplification {

// TODO: Include output file.
void run(grid_reader<height_type> &reader, float contour_interval, float e_z, float e_dp, bool cdp) {
	time_t sec = time(NULL);

	ifstream ifile("topo.tpie");
	bool streams_made = (ifile);
	ifile.close();

	stream<topology_edge> topo_stream("topo.tpie", streams_made ? READ_STREAM : WRITE_STREAM);
	stream<contour_point> unsimplified_stream("unsimplified.tpie", streams_made ? READ_STREAM : WRITE_STREAM);
	stream<contour_point> simplified_stream;

	assert(unsimplified_stream.is_valid());
	assert(topo_stream.is_valid());

	if (!streams_made) { 
		if(contour_interval == 0) {
			cerr << " Creating test." << endl;
			// parent:
			int r = 0;
			int l = 1;
			unsimplified_stream.write_item(cp(0,0,r++,l));
			unsimplified_stream.write_item(cp(0,9,r++,l));
			unsimplified_stream.write_item(cp(16,9,r++,l));
			unsimplified_stream.write_item(cp(16,5,r++,l));
			unsimplified_stream.write_item(cp(16,0,r++,l));
			unsimplified_stream.write_item(cp(0,0,r++,l));
			// to simplify:
			r = 0;
			l++;
			unsimplified_stream.write_item(cp(1,1,r++,l));
			unsimplified_stream.write_item(cp(1,6,r++,l));
			unsimplified_stream.write_item(cp(3,7,r++,l));
			unsimplified_stream.write_item(cp(5,8,r++,l));
			unsimplified_stream.write_item(cp(7,7,r++,l));
			unsimplified_stream.write_item(cp(9,8,r++,l));
			unsimplified_stream.write_item(cp(15,8,r++,l));
			unsimplified_stream.write_item(cp(14,1,r++,l));
			unsimplified_stream.write_item(cp(8,1,r++,l));
			unsimplified_stream.write_item(cp(9,4,r++,l));
			unsimplified_stream.write_item(cp(6,4,r++,l));
			unsimplified_stream.write_item(cp(7,1,r++,l));
			unsimplified_stream.write_item(cp(1,1,r++,l));
			// child 1:
			r = 0;
			l++;
			unsimplified_stream.write_item(cp(2,2,r++,l));
			unsimplified_stream.write_item(cp(4,5,r++,l));
			unsimplified_stream.write_item(cp(6,2,r++,l));
			unsimplified_stream.write_item(cp(2,2,r++,l));
			// child 2:
			r = 0;
			l++;
			unsimplified_stream.write_item(cp(6,5,r++,l));
			unsimplified_stream.write_item(cp(8,7,r++,l));
			unsimplified_stream.write_item(cp(13,7,r++,l));
			unsimplified_stream.write_item(cp(13,5,r++,l));
			unsimplified_stream.write_item(cp(13,2,r++,l));
			unsimplified_stream.write_item(cp(9,2,r++,l));
			unsimplified_stream.write_item(cp(11,5,r++,l));
			unsimplified_stream.write_item(cp(6,5,r++,l));
			// topo:
			topo_stream.write_item(topology_edge(1, -1, 0.9));
			topo_stream.write_item(topology_edge(2, 1, 1));
			topo_stream.write_item(topology_edge(3, 2, 1.1));
			topo_stream.write_item(topology_edge(4, 2, 1.1));
			return;
		}
		contour_reader::read_grid(reader,contour_interval,e_z,topo_stream,unsimplified_stream);
		cerr << " Input streams created. Exiting." << endl;
		cerr << (time(NULL)-sec) << " seconds" << endl;
		return;
	}
	else {
		cerr << " Input streams already created" << endl;
	}

	cerr << "------------ Done read ------------ " << endl;
	cerr << (time(NULL)-sec) << " seconds" << endl;
	sec = time(NULL);

	if (!cdp) {
		cerr << "Running normal Douglas Peucker" << endl;
		cerr << "FATAL ERROR: NO NORMAL DP ANYMORE!";
		//douglas_peucker(e_dp, unsimplified_stream, simplified_stream);
		e_z = contour_interval /1000;
	} else {
		if(contour_interval == 0) {
			contour_interval = 1;
			e_z = 0.1;
		}
		cerr << "Running constrained Douglas Peucker for e=" << e_dp << endl;
		constrained_dp(e_dp, 
					   unsimplified_stream, 
					   topo_stream, 
					   contour_interval, e_z,
					   simplified_stream);
	}

	shape::to_shape("test", unsimplified_stream, topo_stream, &simplified_stream, contour_interval, e_z);
	// transform stream to output:                    
	topo_stream.seek(0);                    
	simplified_stream.seek(0);                      
	stream<ranked_labelled_signed_contour_segment> out_stream;                  
                         
	contour_point *point, prev;                     
	topology_edge *edge = NULL;
	prev.label = -10; // no such label.                    
	elev_t z = 0;
	int rank = 0;
	while(simplified_stream.read_item(&point) == NO_ERROR) {                  
		if(prev.label == point->label) {                    
			out_stream.write_item(ranked_labelled_signed_contour_segment(prev.x,prev.y,              
																		 point->x,point->y,              
																		 z,false,               
																		 prev.label,rank++));             
		}                        
		else {
			ami::err er;
			do {
				er = topo_stream.read_item(&edge);
				assert(er == NO_ERROR);
			} 
			while(edge->c != point->label);
			z = edge->c_z;
			rank = 0;
		}                        
		prev = *point;                      
	}                         
	topo_stream.seek(0);                                        
	simplified_stream.seek(0);                      
//	terrastream::output_mif(out_stream,topo_stream,"mifout");

	cerr << "------------ Done run ------------ " << endl;
	cerr << (time(NULL)-sec) << " additional seconds" << endl;
	sec = time(NULL);
}

}
}
