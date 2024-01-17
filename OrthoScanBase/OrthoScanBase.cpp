#include <iostream>
#include <argparse/argparse.hpp>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Vector_3.h>
#include "../Polyhedron.h"
#include "../MeshFix/MeshFix.h"

using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;

void GenerateBase1(Polyhedron &mesh)
{
    bool upper = true;
    for (auto hv : CGAL::vertices(mesh))
    {
        if (hv->_label >= 31 && hv->_label <= 49)
        {
            upper = false;
            break;
        }
    }

    KernelEpick::Point_3 centroid;
    KernelEpick::Plane_3 plane;
    CGAL::linear_least_squares_fitting_3(mesh.points_begin(), mesh.points_end(), plane, centroid, CGAL::Dimension_tag<0>());
    if (upper && CGAL::angle(plane.orthogonal_vector(), KernelEpick::Vector_3(0, 0, 1)) == CGAL::Angle::ACUTE)
    {
        plane = plane.opposite();
    }
    if (!upper && CGAL::angle(plane.orthogonal_vector(), KernelEpick::Vector_3(0, 0, -1)) == CGAL::Angle::ACUTE)
    {
        plane = plane.opposite();
    }

    double proj_min = std::numeric_limits<double>::max();
    double proj_max = std::numeric_limits<double>::lowest();
    for (auto hv : CGAL::vertices(mesh))
    {
        double proj = CGAL::scalar_product(plane.orthogonal_vector(), hv->point() - centroid);
        proj_min = std::min(proj_min, proj);
        proj_max = std::max(proj_max, proj);
    }
    KernelEpick::Plane_3 clip_plane(centroid + plane.orthogonal_vector() / std::sqrt(plane.orthogonal_vector().squared_length()) * proj_min * 0.5, plane.orthogonal_vector());
    CGAL::Polygon_mesh_processing::clip(mesh, clip_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::allow_self_intersections(true));

    // Add base mesh
    std::vector<Polyhedron::Halfedge_handle> borders;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(borders));

    int max_size = 0;
    Polyhedron::Halfedge_handle hole_handle = nullptr;
    for (auto hole : borders)
    {
        int size = 0;
        for (auto hh : CGAL::halfedges_around_face(hole, mesh))
        {
            size++;
        }
        if (size > max_size)
        {
            hole_handle = hole;
            max_size = size;
        }
    }

    std::vector<Polyhedron::Halfedge_handle> hole_edges;
    for (auto hh : CGAL::halfedges_around_face(hole_handle, mesh))
    {
        hole_edges.push_back(hh);
    }
    std::vector<Polyhedron::Halfedge_handle> new_edges;
    auto h = hole_edges[0];
    auto g = hole_edges[1];
    for (size_t i = 0; i < hole_edges.size(); i++)
    {
        auto ret = mesh.add_vertex_and_facet_to_border(h, g);
        ret->vertex()->point() = h->vertex()->point() - plane.orthogonal_vector() * 3.0;
        ret->vertex()->_label = 1;
        h = ret->opposite();
        g = h->next();
        new_edges.push_back(ret->next()->opposite());
    }

    h = new_edges.front();
    g = h->next()->next();
    for (size_t i = 0; i < new_edges.size(); i++)
    {
        h = mesh.add_facet_to_border(h, g)->opposite();
        g = h->next()->next();
    }
    std::vector<Polyhedron::Facet_handle> patch_faces;
    std::vector<Polyhedron::Vertex_handle> patch_vertex;
    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(mesh, h, std::back_inserter(patch_faces), std::back_inserter(patch_vertex));

    for (auto hv : patch_vertex)
    {
        hv->_label = 1;
    }
}

void GenerateBase2(Polyhedron &mesh)
{
    bool upper = true;
    for (auto hv : CGAL::vertices(mesh))
    {
        if (hv->_label >= 31 && hv->_label <= 49)
        {
            upper = false;
            break;
        }
    }
    KernelEpick::Point_3 centroid;
    KernelEpick::Plane_3 plane;
    CGAL::linear_least_squares_fitting_3(mesh.points_begin(), mesh.points_end(), plane, centroid, CGAL::Dimension_tag<0>());
    if (upper && CGAL::angle(plane.orthogonal_vector(), KernelEpick::Vector_3(0, 0, 1)) == CGAL::Angle::ACUTE)
    {
        plane = plane.opposite();
    }
    if (!upper && CGAL::angle(plane.orthogonal_vector(), KernelEpick::Vector_3(0, 0, -1)) == CGAL::Angle::ACUTE)
    {
        plane = plane.opposite();
    }

    double proj_min = std::numeric_limits<double>::max();
    double proj_max = std::numeric_limits<double>::lowest();
    for (auto hv : CGAL::vertices(mesh))
    {
        double proj = CGAL::scalar_product(plane.orthogonal_vector(), hv->point() - centroid);
        proj_min = std::min(proj_min, proj);
        proj_max = std::max(proj_max, proj);
    }
    using LA = CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<typename CGAL::Eigen_sparse_matrix<double>::EigenType>>;
    using HeatMethod = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Polyhedron, CGAL::Heat_method_3::Intrinsic_Delaunay,
                                                                              typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::const_type, LA, Polyhedron::Traits::Kernel>;
    std::unordered_map<typename Polyhedron::Vertex_handle, double> distances;
    HeatMethod heat_method(mesh);

    std::vector<typename Polyhedron::Vertex_handle> sources;
    for (auto hh : CGAL::halfedges(mesh))
    {
        if (hh->facet() != nullptr && hh->opposite() != nullptr && hh->opposite()->facet() != nullptr)
        {
            if (hh->facet()->_label == 0 && hh->opposite()->facet()->_label != 0)
            {
                sources.push_back(hh->vertex());
            }
        }
    }
    heat_method.add_sources(sources);
    heat_method.estimate_geodesic_distances(boost::make_assoc_property_map(distances));

    std::unordered_set<typename Polyhedron::Facet_handle> face_to_remove;
    for (auto hv : CGAL::vertices(mesh))
    {
        int l = hv->_label;
        if (!(l >= 11 && l <= 29 || l >= 31 && l <= 49))
        {
            if (distances[hv] > (proj_max - proj_min) * 0.2)
            {
                for (auto hf : CGAL::faces_around_target(hv->halfedge(), mesh))
                {
                    face_to_remove.insert(hf);
                }
            }
        }
    }

    for (auto hf : face_to_remove)
    {
        mesh.erase_facet(hf->halfedge());
    }

    auto [vertices, faces] = mesh.ToVerticesTriangles();
    FixMeshWithLabel(vertices, faces, mesh.WriteLabels(), mesh, true, 1000, false, true, 0, 0, false, 10);

    // Add base mesh
    std::vector<Polyhedron::Halfedge_handle> borders;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(borders));

    int max_size = 0;
    Polyhedron::Halfedge_handle hole_handle = nullptr;
    for (auto hole : borders)
    {
        int size = 0;
        for (auto hh : CGAL::halfedges_around_face(hole, mesh))
        {
            size++;
        }
        if (size > max_size)
        {
            hole_handle = hole;
            max_size = size;
        }
    }

    std::vector<Polyhedron::Halfedge_handle> hole_edges;
    for (auto hh : CGAL::halfedges_around_face(hole_handle, mesh))
    {
        hole_edges.push_back(hh);
    }
    std::vector<Polyhedron::Halfedge_handle> new_edges;
    auto h = hole_edges[0];
    auto g = hole_edges[1];
    for (size_t i = 0; i < hole_edges.size(); i++)
    {
        auto ret = mesh.add_vertex_and_facet_to_border(h, g);
        double proj = CGAL::scalar_product(plane.orthogonal_vector(), h->vertex()->point() - centroid);
        ret->vertex()->point() = h->vertex()->point() - plane.orthogonal_vector() * (proj - proj_min);
        ret->vertex()->_label = 1;
        h = ret->opposite();
        g = h->next();
        new_edges.push_back(ret->next()->opposite());
    }

    h = new_edges.front();
    g = h->next()->next();
    for (size_t i = 0; i < new_edges.size(); i++)
    {
        h = mesh.add_facet_to_border(h, g)->opposite();
        g = h->next()->next();
    }
    std::vector<Polyhedron::Facet_handle> patch_faces;
    std::vector<Polyhedron::Vertex_handle> patch_vertex;
    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(mesh, h, std::back_inserter(patch_faces), std::back_inserter(patch_vertex));

    for (auto hv : patch_vertex)
    {
        hv->_label = 1;
    }

    borders.clear();
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(borders));
    for (auto hh : borders)
    {
        std::vector<Polyhedron::Facet_handle> patch_faces;
        std::vector<Polyhedron::Vertex_handle> patch_vertex;
        CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(mesh, hh, std::back_inserter(patch_faces), std::back_inserter(patch_vertex));
    }
}

void Optimize(Polyhedron &mesh)
{
    struct Visitor : public CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<Polyhedron>
    {
        std::vector<typename Polyhedron::Facet_handle> new_faces;
        void before_subface_creations(face_descriptor fd)
        {
        }
        void after_subface_created(face_descriptor fd)
        {
            new_faces.push_back(fd);
        }
        void after_subface_creations()
        {
        }
    };
    double avg_area = 0.0;
    for(auto hf : CGAL::faces(mesh))
    {
        avg_area += std::sqrt(CGAL::squared_area(hf->halfedge()->vertex()->point(), hf->halfedge()->next()->vertex()->point(), hf->halfedge()->next()->next()->vertex()->point()));
    }
    avg_area /= mesh.size_of_facets();

    std::vector<typename Polyhedron::Facet_handle> faces_to_split;
    for(auto hf : CGAL::faces(mesh))
    {
        if(hf->halfedge()->vertex()->_label == 0 && hf->halfedge()->next()->vertex()->_label == 0 && hf->halfedge()->next()->next()->vertex()->_label == 0)
        {
            double area = std::sqrt(CGAL::squared_area(hf->halfedge()->vertex()->point(), hf->halfedge()->next()->vertex()->point(), hf->halfedge()->next()->next()->vertex()->point()));
            if(area > avg_area * 4)
            {
                faces_to_split.push_back(hf);
            }
        }
    }

    std::vector<typename Polyhedron::Facet_handle> new_faces;
    std::vector<typename Polyhedron::Vertex_handle> new_vertices;
    CGAL::Polygon_mesh_processing::refine(mesh, faces_to_split, std::back_inserter(new_faces), std::back_inserter(new_vertices));
    for(auto hv : new_vertices)
    {
        hv->_label = 0;
    }
}

int main(int argc, char *argv[])
{

    argparse::ArgumentParser parser;
    parser.add_argument("--input_file", "-i").required().help("");
    parser.add_argument("--input_label", "-l").required().help("");
    parser.add_argument("--output_file", "-o").required().help("");
    parser.add_argument("--output_label", "-ol").required().help("");
    parser.parse_args(argc, argv);

    Polyhedron mesh;
    CGAL::IO::read_polygon_mesh(parser.get("-i"), mesh, CGAL::parameters::verbose(true));
    mesh.LoadLabels(parser.get("-l"));
    Optimize(mesh);
    GenerateBase2(mesh);

    mesh.WriteOBJ(parser.get("-o"));
    mesh.WriteLabels(parser.get("-ol"));
    return 0;
}