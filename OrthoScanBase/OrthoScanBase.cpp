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
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
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
        if(hv->_label >= 11 && hv->_label <= 29)
        {
            upper = true;
            break;
        }
    }
    std::cout << "Computing principle axis...";
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

    std::cout << "Computing geodesic distance...";
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
    
    std::cout << "Erasing faces...";
    std::unordered_set<typename Polyhedron::Facet_handle> face_to_remove;
    for (auto hv : CGAL::vertices(mesh))
    {
        int l = hv->_label;
        if (!(l >= 11 && l <= 29 || l >= 31 && l <= 49))
        {
            if (distances[hv] > (proj_max - proj_min) * 0.05)
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
        if(hf != nullptr && hf->halfedge() != nullptr && !hf->halfedge()->is_border())
        {
            mesh.erase_facet(hf->halfedge());
        }
    }

    auto [vertices, faces] = mesh.ToVerticesTriangles();
    FixMeshWithLabel(vertices, faces, mesh.WriteLabels(), mesh, true, 1000, false, true, 0, 0, false, 10);

    // Add base mesh
    std::cout << "Building mesh...";
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
    std::vector<Polyhedron::Facet_handle> new_faces;
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
        new_faces.push_back(ret->facet());
    }

    h = new_edges.front();
    g = h->next()->next();
    for (size_t i = 0; i < new_edges.size(); i++)
    {
        auto ret = mesh.add_facet_to_border(h, g);
        h = ret->opposite();
        g = h->next()->next();
        new_faces.push_back(ret->facet());
    }

    std::vector<Polyhedron::Vertex_handle> hole_vertices;
    for(auto hh : CGAL::halfedges_around_face(h, mesh))
    {
        hole_vertices.push_back(hh->vertex());
    }
    {
        for(int ite = 0; ite < 3; ite++)
        {
            std::vector<Polyhedron::Traits::Kernel::Point_3> new_points;
            for(size_t i = 0; i < hole_vertices.size(); i++)
            {
                size_t prev = i == 0 ? hole_vertices.size() - 1 : i - 1;
                size_t next = (i + 1) % hole_vertices.size();
                new_points.push_back(CGAL::midpoint(hole_vertices[prev]->point(), hole_vertices[next]->point()));
            }
            for(size_t i = 0; i < hole_vertices.size(); i++)
            {
                hole_vertices[i]->point() = new_points[i];
            }
        }
    }

    double target_len = 0.0;
    int cnt = 0;
    for(auto hh : CGAL::halfedges_around_face(h, mesh))
    {
        target_len += std::sqrt(CGAL::squared_distance(hh->vertex()->point(), hh->next()->vertex()->point()));
        cnt++;
    }
    target_len = target_len / cnt * 4.0;
    std::unordered_set<Polyhedron::Facet_handle> old_faces;
    for(auto hf : CGAL::faces(mesh))
        old_faces.insert(hf);
    
    std::cout << "Remeshing..." << std::endl;
    CGAL::Polygon_mesh_processing::isotropic_remeshing(new_faces, target_len, mesh, CGAL::parameters::relax_constraints(true).number_of_relaxation_steps(3).number_of_iterations(3));
    std::unordered_set<Polyhedron::Vertex_handle> vertex_to_fair;
    new_faces.clear();
    for(auto hf : CGAL::faces(mesh))
    {
        if(old_faces.count(hf) == 0)
        {
            new_faces.push_back(hf);
            for(auto hh : CGAL::halfedges_around_face(hf->halfedge(), mesh))
                vertex_to_fair.insert(hh->vertex());
        }
    }
    for(auto hf : new_faces)
    {
        for(auto hh : CGAL::halfedges_around_face(hf->halfedge(), mesh))
        {
            if(hh->opposite()->is_border())
                vertex_to_fair.erase(hh->vertex());
        }
    }
    CGAL::Polygon_mesh_processing::fair(mesh, std::vector<Polyhedron::Vertex_handle>(vertex_to_fair.begin(), vertex_to_fair.end()), CGAL::parameters::number_of_iterations(1));

    Polyhedron::Halfedge_handle hole_hh = nullptr;
    for(auto hf : new_faces)
    {
        for(auto hh : CGAL::halfedges_around_face(hf->halfedge(), mesh))
        {
            if(hh->opposite()->is_border())
            {
                hole_hh = hh->opposite();
                break;
            }
        }
        if(hole_hh != nullptr)
            break;
    }
    std::vector<Polyhedron::Facet_handle> patch_faces;
    std::vector<Polyhedron::Vertex_handle> patch_vertex;
    CGAL::Polygon_mesh_processing::triangulate_hole(mesh, hole_hh, std::back_inserter(patch_faces));
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
    CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces(mesh, CGAL::parameters::needle_threshold(100).cap_threshold(std::cos(3.14159 * 0.9)));
    double avg_area = 0.0;
    double avg_len = 0.0;
    for(auto hf : CGAL::faces(mesh))
    {
        auto p0 = hf->halfedge()->vertex()->point();
        auto p1 = hf->halfedge()->next()->vertex()->point();
        auto p2 = hf->halfedge()->next()->next()->vertex()->point();
        avg_area += std::sqrt(CGAL::squared_area(p0, p1, p2));
        avg_len += std::sqrt(CGAL::squared_distance(p0, p1));
        avg_len += std::sqrt(CGAL::squared_distance(p0, p2));
        avg_len += std::sqrt(CGAL::squared_distance(p1, p2));
    }
    avg_area /= mesh.size_of_facets();
    avg_len /= mesh.size_of_halfedges();

    std::unordered_set<typename Polyhedron::Facet_handle> faces_to_split;
    double threshold = avg_len * avg_len * 16;
    for(auto hf : CGAL::faces(mesh))
    {
        int cnt = 0;
        if(hf->halfedge()->vertex()->_label == 0) cnt++;
        if(hf->halfedge()->next()->vertex()->_label == 0) cnt++;
        if(hf->halfedge()->next()->next()->vertex()->_label == 0) cnt++;
        if(cnt >= 1)
        {
            auto p0 = hf->halfedge()->vertex()->point();
            auto p1 = hf->halfedge()->next()->vertex()->point();
            auto p2 = hf->halfedge()->next()->next()->vertex()->point();
            if(CGAL::squared_distance(p0, p1) > threshold || CGAL::squared_distance(p0, p2) > threshold || CGAL::squared_distance(p1, p2) > threshold)
            {
                faces_to_split.insert(hf);
            }
        }
    }

    std::vector<std::vector<typename Polyhedron::Facet_handle>> face_patches;
    while(!faces_to_split.empty())
    {
        auto ff = *faces_to_split.begin();
        faces_to_split.erase(ff);
        face_patches.emplace_back();
        std::queue<typename Polyhedron::Facet_handle> q;
        q.push(ff);
        while(!q.empty())
        {
            auto hf = q.front();
            q.pop();
            face_patches.back().push_back(hf);
            for(auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
            {
                if(faces_to_split.count(nei) != 0)
                {
                    q.push(nei);
                    faces_to_split.erase(nei);
                }
            }
        }
    }
    std::cout << "processing " << face_patches.size() << " patches.";
    for(auto& patch : face_patches)
    {
        CGAL::Polygon_mesh_processing::isotropic_remeshing(patch, avg_len * 2, mesh, CGAL::parameters::number_of_iterations(1));
    }
}

int main(int argc, char *argv[])
{
    argparse::ArgumentParser parser;
    parser.add_argument("--input_file", "-i").required().help("specify the input mesh.");
    parser.add_argument("--input_label", "-l").required().help("specify the input labels.");
    parser.add_argument("--output_file", "-o").required().help("specify the output file.");
    parser.add_argument("--output_label", "-ol").required().help("specify the output labels.");
    parser.parse_args(argc, argv);

    Polyhedron mesh;
    CGAL::IO::read_polygon_mesh(parser.get("-i"), mesh, CGAL::parameters::verbose(true));
    mesh.LoadLabels(parser.get("-l"));
    try
    {
        std::cout << "Optimizing...";
        Optimize(mesh);
        std::cout << "Done." << std::endl;
        std::cout << "Generating...";
        GenerateBase2(mesh);
        std::cout << "Done." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    mesh.WriteOBJ(parser.get("-o"));
    mesh.WriteLabels(parser.get("-ol"));
    return 0;
}