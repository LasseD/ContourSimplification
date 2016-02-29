// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "mif_outputter.h"
#include <math.h>
#include <terrastream/common/sort.h>

using namespace std;
using namespace terrastream;
typedef segment_point point;
typedef ranked_labelled_signed_contour_segment con_segment;

struct cmp_topo_edge_child{
	inline int compare(const topology_edge &t1, const topology_edge &t2) {
		return t1.c == t2.c ? 0 : (t1.c < t2.c ? -1 : 1);
	}
};

///////////////////////////////////////////////////////////////
/// Figures out at which end the new segment fits the previous
/// point, and adds the right point to the input vector
///////////////////////////////////////////////////////////////
inline void add_next_point(vector<point>& points, con_segment* segment) {
	point p1(segment->x1, segment->y1), p2(segment->x2, segment->y2);
	if (points.back() == p1)
		points.push_back(p2);
	else {
		assert(points.back() == p2);
		points.push_back(p1);
	}
}

inline void write_contour(vector<point>& points, ofstream& mif) {
	mif << "Pline " << points.size() << "\n";
	for (vector<point>::iterator it = points.begin(); it != points.end(); ++it)
		mif << it->x << " " << it->y << "\n";
	mif << "    Pen (1,2,0)\n";
	mif << "    Smooth\n";
	points.clear();
}


void terrastream::output_mif(stream<con_segment>& stream_in, stream<topology_edge> &topos, string basename) {

	cmp_topo_edge_child cmp;
	topos.seek(0);
//	ts_sort(&topos, &cmp);
	topos.seek(0);

	log_info() << "Outputting .mif and .mid file\n";
	ofstream mif((basename + ".mif").c_str());
	ofstream mid((basename + ".mid").c_str());
	//Write header
	mif << "Version 300\n";
	mif << "Charset \"WindowsLatin1\"\n";
	mif << "Delimiter \" \"\n";
	mif << "Columns 3 \n";
	mif << "  Kote Float\n";
	mif << "  Label Integer\n";
	mif << "  Parent Integer\n";
	mif << "Data\n\n";
	mif << fixed << setprecision(8);
	mif.imbue(locale("C"));
	mid.imbue(locale("C"));

	con_segment* segment;
	vector<point> points;
	tpie::ami::err err;
	stream_in.seek(0);
	while ((err = stream_in.read_item(&segment)) == tpie::ami::NO_ERROR) {
		if (segment->rank ==0) { //Reached end of previous contour
			if (segment->label != 1) { //Nothing to output before first contour is read
				write_contour(points, mif);
			}
			//Add both points for starting segment of new contour
			points.push_back(point(segment->x1, segment->y1));
			points.push_back(point(segment->x2, segment->y2));
			//Write data to .mid file.
			topology_edge* top_edge;
			while(topos.read_item(&top_edge) == tpie::ami::NO_ERROR && top_edge->c != segment->label) {
//				throw std::runtime_error("An error occured while reading from topology edge stream");
			}
			assert(top_edge->c == segment->label);
			mid << segment->z << " " << segment->label << " " << top_edge->p << "\n";
		} else {//Figure out which point from the segment to add
			add_next_point(points, segment);
		}
	}
	if (err != tpie::ami::END_OF_STREAM) {
		stringstream ss;
		ss << "An error occured while reading from segment stream\n"
		   << "Errorcode was: " << int(err) << "\n";
		throw std::runtime_error(ss.str());
	}
	//Write last contour
	write_contour(points, mif);
	stream_in.seek(0);
}
