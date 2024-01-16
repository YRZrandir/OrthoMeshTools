#include <iostream>
#include <argparse/argparse.hpp>
#include <CGAL/boost/graph/io.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Vector_3.h>
#include "../Polyhedron.h"

int main(int argc, char* argv[])
{
    using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;
    
    argparse::ArgumentParser parser;
    parser.add_argument("--input_file", "-i").required().help("");
    parser.add_argument("--input_label", "-l").required().help("");
    parser.add_argument("--output_file", "-o").required().help("");
    parser.add_argument("--output_label", "-ol").required().help("");
    parser.parse_args(argc, argv);

    Polyhedron mesh;
    CGAL::IO::read_polygon_mesh(parser.get("-i"), mesh, CGAL::parameters::verbose(true));
    mesh.LoadLabels(parser.get("-l"));

    bool upper = true;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label >= 31 && hv->_label <= 49)
        {
            upper = false;
            break;
        }
    }
    
    KernelEpick::Point_3 centroid;
    KernelEpick::Plane_3 plane;
    CGAL::linear_least_squares_fitting_3( mesh.points_begin(), mesh.points_end(), plane, centroid, CGAL::Dimension_tag<0>());
    if(upper && CGAL::angle(plane.orthogonal_vector(), KernelEpick::Vector_3(0, 0, 1)) == CGAL::Angle::ACUTE)
    {
        plane = plane.opposite();
    }
    if(!upper && CGAL::angle(plane.orthogonal_vector(), KernelEpick::Vector_3(0, 0, -1)) == CGAL::Angle::ACUTE)
    {
        plane = plane.opposite();
    }

    double proj_min = std::numeric_limits<double>::max();
    double proj_max = std::numeric_limits<double>::lowest();
    for(auto hv : CGAL::vertices(mesh))
    {
        double proj = CGAL::scalar_product(plane.orthogonal_vector(), hv->point() - centroid);
        proj_min = std::min(proj_min, proj);
        proj_max = std::max(proj_max, proj);
    }
    KernelEpick::Plane_3 clip_plane(centroid + plane.orthogonal_vector() / std::sqrt(plane.orthogonal_vector().squared_length()) * proj_min * 0.5, plane.orthogonal_vector());
    CGAL::Polygon_mesh_processing::clip(mesh, clip_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::allow_self_intersections(true));

    mesh.WriteOBJ(parser.get("-o"));
    return 0;
}