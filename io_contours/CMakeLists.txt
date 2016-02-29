set(GRID_SRCS)

set(GRID_HDRS)

set(TIN_SRCS)

set(TIN_HDRS)
	
set(COMMON_HDRS
	intersect.h
	contour_types.h
	ridge_removal.h
	outer_curves.h
	cw_ordering.h
	topology.h
	contours.h
	tin_to_triangle.h

	mif_outputter.h
)

set(COMMON_SRCS
	intersect.cpp
	contour_types.cpp	
	ridge_removal.cpp
	outer_curves.cpp
	cw_ordering.cpp
	topology.cpp
	contours.cpp
	tin_to_triangle.cpp

	mif_outputter.cpp
)

add_terrastream_library(
	io_contours_test
	COMMON_SRCS
	COMMON_HDRS
	TIN_SRCS
	TIN_HDRS
	GRID_SRCS
	GRID_HDRS
)

