add_subdirectory(io_contours) 

set(GRID_SRCS)

set(GRID_HDRS)

set(COMMON_HDRS
  decomposition.h
  contour_simplification.h
  contour_reader.h
  contour_to_shape.h
  util.h
  simplify.h
)

set(COMMON_SRCS
	decomposition.cpp
	contour_simplification.cpp	
	contour_reader.cpp
	contour_to_shape.cpp
	simplify.cpp
)

set(TIN_SRCS
)

set(TIN_HDRS
)

add_terrastream_library(
  io_contour_simplification
	COMMON_SRCS
	COMMON_HDRS
	TIN_SRCS
	TIN_HDRS
	GRID_SRCS
	GRID_HDRS
)

