// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "contour_reader.h"
#include <math.h>
#include "io_contours/contour_types.h"
#include "io_contours/contours.h"
#include <terrastream/common/nodata.h>
#include "../../terrastream/common/raster_drivers/gdal_grid_reader.h"

//#define DEBUG_CONTOUR_READER

using namespace terrastream;
using namespace tpie::ami;
using namespace std;

typedef ranked_labelled_signed_contour_segment rlss;

static int BORDER_ELEV = -1;

inline height_type row(height_type *row, int const x, int const width) {
	assert(x >=0);
	assert(x <= width+1);
	if (x == 0 || x == width+1)
		return BORDER_ELEV;
	height_type res = row[x-1];
	if(is_nodata(res)) // no-val. TODO FIXME: Use real no val
		return BORDER_ELEV;
	return row[x-1];
}

// width: non-bordered width.
void make_triangles(stream<triangle> &tris, 
					height_type *row2, height_type *row1,
					int const width, int const y) {
	for (int x = 0; x < width+1; x++) {
		triangle_point share1(x,y+1,row(row2, x, width));
		triangle_point share2(x+1,y,row(row1, x+1, width));
		// lower triangle: diagonal, cw. (opposite)
		tris.write_item(triangle(triangle_point(x,y,row(row1, x, width)),
								 share1,share2));
		// upper triangle: diagonal, cw.
		tris.write_item(triangle(triangle_point(x+1,y+1,row(row2, x+1, width)),
								 share2,share1));
	}
}

// void contour_reader::read_grid(char *file,
// 							   float const contour_interval,
// 							   float const e_z,
// 							   stream<topology_edge>& os_topo,
// 							   stream<rlss>& os_segs) {
// 	boost::program_options::variables_map varmap;
// 	grid_reader<height_type> reader;
// 	reader.open(file, -9999, varmap);
// 	read_grid(reader, contour_interval, e_z, os_topo, os_segs);
// }

void contour_reader::read_grid(grid_reader<height_type> &reader,
							   float const contour_interval,
							   float const e_z,
							   stream<topology_edge>& os_topo,
							   stream<contour_point>& os_segs) {
#ifdef DEBUG_CONTOUR_READER
	cerr << "Starting to read grid from file" << endl;
#endif
	os_segs.truncate(0);
	os_topo.truncate(0);
	os_segs.seek(0);
	os_topo.seek(0);
//	  cerr << "Options: " << reader.get_options() << endl;
	int width = reader.get_ncols();
#ifdef DEBUG_CONTOUR_READER
	cerr << "Width: " << width << endl;
#endif

	height_type* row1 = new height_type[width];
	height_type* row2 = new height_type[width];
	for (int i = 0; i < width; i++)
		row2[i]=BORDER_ELEV; // border row.

	stream<triangle> tris;

	int y = 0;
	bool b = true;
	while (reader.next_row(b ? row1 : row2)) {
		if (b)
			make_triangles(tris, row1, row2, width, y);
		else
			make_triangles(tris, row2, row1, width, y);
		// Prepare next:
		y++;
		b = !b;
	}
	for (int i = 0; i < width; i++) {
		if (b)
			row1[i]=BORDER_ELEV; // border row.
		else
			row2[i]=BORDER_ELEV; // border row.
	}
	if (b)
		make_triangles(tris, row1, row2, width, y);
	else
		make_triangles(tris, row2, row1, width, y);

	tris.seek(0);
	delete []row1;
	delete []row2;

#ifdef DEBUG_CONTOUR_READER
	cerr << "Triangles constructed. Moving on to contour lines"  << endl;
#endif
	compute_contours(tris,contour_interval, e_z, os_segs, os_topo); 
#ifdef DEBUG_CONTOUR_READER
	cerr << "Done contour lines"  << endl;
#endif
}
