// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

#include "tin_to_triangle.h"
#include <terrastream/common/tin_io.h>
#include <terrastream/common/tin_reader.h>


using namespace std;
using namespace terrastream;


template <typename Id_type>
struct cmp_triangle_ids{
	inline int compare(const Id_type &i1, const Id_type &i2) {
		return i1.next_id < i2.next_id ? -1 : (i1.next_id > i2.next_id ? 1 : 0);
	}
};

struct two_points_one_id {
	triangle_point p1, p2;
	tin_vertex_id next_id;

	triangle create_next(triangle_point p3) {return triangle(p1, p2, p3);}
};

struct one_point_two_ids {
	triangle_point p1;
	tin_vertex_id next_id, id3;

	two_points_one_id create_next(triangle_point p2) {
		two_points_one_id tmp ={p1, p2, id3};
		return tmp;
	}
};

struct three_ids {
	tin_vertex_id next_id, id2, id3;

	one_point_two_ids create_next(triangle_point p1) {
		one_point_two_ids temp={p1, id2, id3};
			return temp;
	}
};

template <typename StreamT,typename Cmp>
void sort(stream<StreamT> &stream) {
	tflow_progress progress("Sorting", "Sorting", 0, stream.stream_len(), 1);
	Cmp cmp;
	stream.seek(0);
	ts_sort(&stream, &cmp, &progress);
	stream.seek(0);
}

template <typename In_t, typename Out_t>
void assign_points(stream<In_t> &in, stream<Out_t> &out, stream<tin_node<elev_t> > &nodes) {
	nodes.seek(0);
	tin_node<elev_t>* node;
	In_t* in_t;
	tpie::ami::err err;
	tpie::ami::err err_node;
	tflow_progress progress("Assigning", "Assigning", 0, in.stream_len(), 1);
	while ((err = in.read_item(&in_t)) == tpie::ami::NO_ERROR) {
		while (in_t->next_id > nodes.tell()-1) {
			if ((err_node = nodes.read_item(&node)) != tpie::ami::NO_ERROR) {
				stringstream ss;
				ss << "Couldn't read node nr: "
				   << nodes.tell() << " Error nr: "
				   << int(err_node) << "\n";
				throw std::runtime_error(ss.str());
			}
		}
		Out_t output = in_t->create_next(triangle_point(node->x, node->y, node->z));
		if (out.write_item(output) != tpie::ami::NO_ERROR)
			throw std::runtime_error("An error occured while writing to stream!");
		progress.step();
	}
	if (err != tpie::ami::END_OF_STREAM) {
		stringstream ss;
		ss << "Error while reading from temporary stream.\n Error code is: "
			 << int(err) << "\n";
		throw std::runtime_error(ss.str());
	}
	progress.done();
}

void terrastream::tin_to_triangle(tin_reader<elev_t> &in, stream<triangle> &out) {

	log_info() << "Converting triangles into stream\n";
	stream<three_ids> three_ids_stream;
	tflow_progress progress("Converting", "Converting", 0, in.get_triangle_upperbound(), 1);
	tin_triangle triangle_in;
	while (in.next_triangle(triangle_in)) {
		three_ids tmp = {triangle_in.nodes[0], triangle_in.nodes[1], triangle_in.nodes[2]};
		if (three_ids_stream.write_item(tmp) != tpie::ami::NO_ERROR)
			throw std::runtime_error("An error occured while writing to stream!");
		progress.step();
	}
	progress.done();

	log_info() << "Converting nodes into stream\n";
	stream<tin_node<elev_t> > node_stream;
	tflow_progress progress2("Converting", "Converting", 0, in.get_node_count(), 1);
	tin_node<elev_t> node_in;
	while (in.next_node(node_in)) {
		if (node_stream.write_item(node_in) != tpie::ami::NO_ERROR)
			throw std::runtime_error("An error occured while writing to stream!");
		progress2.step();
	}
	progress2.done();

	log_info() << "Sorting by first id\n";
	sort<three_ids, cmp_triangle_ids<three_ids> >(three_ids_stream);

	log_info() << "Assigning first point to triangles\n";
	stream<one_point_two_ids> one_point_stream;
	assign_points(three_ids_stream, one_point_stream, node_stream);

	log_info() << "Sorting by second id\n";
	sort<one_point_two_ids, cmp_triangle_ids<one_point_two_ids> >(one_point_stream);

	log_info() << "Assigning second point to triangles\n";
	stream<two_points_one_id> two_points_stream;
	assign_points(one_point_stream, two_points_stream, node_stream);

	log_info() << "Sorting by third id\n";
	sort<two_points_one_id, cmp_triangle_ids<two_points_one_id> >(two_points_stream);

	log_info() << "Assigning third point to triangles\n";
	assign_points(two_points_stream, out, node_stream);
}
