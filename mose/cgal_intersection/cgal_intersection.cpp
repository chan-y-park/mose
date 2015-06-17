// Computing intersection points among curves using the sweep line.
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <list>

#include <fstream>
#include "cgal_intersection.h"

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Curve_2                               Segment_2;
int cgal(){
        // Construct the input segments.
        Segment_2 segments[] = {Segment_2 (Point_2 (1, 5), Point_2 (8, 5)),
                              Segment_2 (Point_2 (1, 1), Point_2 (8, 8)),
                              Segment_2 (Point_2 (3, 1), Point_2 (3, 8)),
                              Segment_2 (Point_2 (8, 5), Point_2 (8, 8))};

        // Compute all intersection points.
        std::list<Point_2>     pts;
        CGAL::compute_intersection_points (segments, segments + 4,
                                         std::back_inserter (pts));

        // Print the result.
        std::cout << "Found " << pts.size() << " intersection points: " 
            << std::endl; 
        std::copy (pts.begin(), pts.end(), 
            std::ostream_iterator<Point_2>(std::cout, "\n"));
        // Compute the non-intersecting sub-segments 
        // induced by the input segments.
        std::list<Segment_2>   sub_segs;
        CGAL::compute_subcurves(segments, segments + 4, 
            std::back_inserter(sub_segs));
        std::cout << "Found " << sub_segs.size()
                << " interior-disjoint sub-segments." << std::endl;
        CGAL_assertion (CGAL::do_curves_intersect (segments, segments + 4));
        return 0;
    }

extern "C" int find_intersections_of_segments(
    coordinate *segment_1, long segment_1_size,
    coordinate *segment_2, long segment_2_size,
    coordinate *intersections, int max_num_of_intersections)
{

    long i_1, i_2;
    int num_of_intersections = 5;

    std::ofstream fout("segment_coordinates.txt");

    fout << "Segment 1:\n";
    for (i_1 = 0; i_1 < segment_1_size; i_1++)
        fout << segment_1[i_1].x << ", " << segment_1[i_1].y << "\n";
    
    fout << "========================================\n";

    fout << "Segment 2:\n";
    for (i_2 = 0; i_2 < segment_2_size; i_2++)
        fout << segment_2[i_2].x << ", " << segment_1[i_2].y << "\n";

    fout.close();

    for (int j = 0; j < num_of_intersections; j++){
        if (j == max_num_of_intersections)
            break;
        intersections[j].x = (j / 2.0);
        intersections[j].y = (j / 3.0);
    }

    // If the buffer size is too small, return the number of intersections.
    if (num_of_intersections > max_num_of_intersections)
        return num_of_intersections;
    else
        return 0;
}
