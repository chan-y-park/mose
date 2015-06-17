typedef struct {double x; double y;} coordinate;

int cgal();

extern "C" int find_intersections_of_segments(
    coordinate *segment_1, long segment_1_size,
    coordinate *segment_2, long segment_2_size,
    coordinate *intersections, int max_num_of_intersections);
