#include <iostream>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Surface_mesh_deformation.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Vector_3.h>
#include "../Polyhedron.h"
#include "../MeshFix/MeshFix.h"
#include "../EasyOBJ.h"
#include "../SegClean/SegClean.h"

#ifndef FOUND_PYBIND11
#include <argparse/argparse.hpp>
#endif

//#define DEBUG_ORTHOSCANBASE

using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;

void LabelProcessing(Polyhedron& mesh)
{
    using hVertex = typename Polyhedron::Vertex_handle;
    auto aabb = CGAL::bbox_3(mesh.points_begin(), mesh.points_end());
    double threshold = std::max(aabb.x_span(), std::max(aabb.y_span(), aabb.z_span())) / 150.0;
    
    std::unordered_map<hVertex, int> new_label_set;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label != 0)
        {
            continue;
        }
        std::unordered_set<hVertex> neighbors;
        std::unordered_set<int> labels;
        std::queue<hVertex> q;
        q.push(hv);
        neighbors.insert(hv);
        labels.insert(hv->_label);
        while(!q.empty())
        {
            hVertex curr = q.front();
            q.pop();
            for(auto nei : CGAL::vertices_around_target(curr, mesh))
            {
                if(neighbors.count(nei) == 0 && CGAL::squared_distance(nei->point(), hv->point()) < threshold * threshold)
                {
                    neighbors.insert(nei);
                    labels.insert(nei->_label);
                    q.push(nei);
                }
            }
        }
        if(labels.size() >= 3)
        {
            double nearest_dist = std::numeric_limits<double>::max();
            hVertex nearest_hv = nullptr;
            for(auto nei : neighbors)
            {
                if(nei->_label == hv->_label)
                {
                    continue;
                }
                double d = CGAL::squared_distance(nei->point(), hv->point());
                if(d < nearest_dist)
                {
                    nearest_dist = d;
                    nearest_hv = nei;
                }
            }
            new_label_set.insert({hv, nearest_hv->_label});
        }
    }

    for(auto [hv, label] : new_label_set)
    {
        hv->_label = label;
    }
}

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
    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(mesh, h,
        CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)).vertex_output_iterator(std::back_inserter(patch_vertex)));
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
    using HeatMethod = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Polyhedron, CGAL::Heat_method_3::Direct,
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
            double dist = std::numeric_limits<double>::max();
            if(distances.find(hv) != distances.end())
            {
                dist = distances[hv];
            }
            if (dist > (proj_max - proj_min) * 0.05)
            {
                for (auto hf : CGAL::faces_around_target(hv->halfedge(), mesh))
                {
                    face_to_remove.insert(hf);
                }
            }
        }
    }

#ifdef DEBUG_ORTHOSCANBASE
    {
        std::vector<Polyhedron::Point_3> points;
        std::vector<Eigen::Vector3f> colors;
        std::vector<TTriangle<size_t>> faces;
        int nv = 0;
        int nt = 0;
        for(auto hf : CGAL::faces(mesh))
        {
            auto p0 = hf->halfedge()->vertex()->point();
            auto p1 = hf->halfedge()->next()->vertex()->point();
            auto p2 = hf->halfedge()->prev()->vertex()->point();
            points.push_back(p0);
            points.push_back(p1);
            points.push_back(p2);
            faces.emplace_back(nv, nv + 1, nv + 2);
            nv += 3;
            if(face_to_remove.count(hf) != 0)
            {
                colors.push_back(Eigen::Vector3f(1, 0, 0));
                colors.push_back(Eigen::Vector3f(1, 0, 0));
                colors.push_back(Eigen::Vector3f(1, 0, 0));
            }
            else
            {
                colors.push_back(Eigen::Vector3f(0.7f, 0.7f, 0.7f));
                colors.push_back(Eigen::Vector3f(0.7f, 0.7f, 0.7f));
                colors.push_back(Eigen::Vector3f(0.7f, 0.7f, 0.7f));
            }
        }
        WriteVCFAssimp<Polyhedron::Traits::Kernel, size_t>("removed1.obj", points, colors, faces);
    }
#endif

    // After removing, we want the remaining part to be one single connected component. So we try to merge them here.
    std::cout << "Computing connected components..." << std::endl;
    std::unordered_set<typename Polyhedron::Facet_handle> remain_faces;
    for(auto hf : CGAL::faces(mesh))
        if(face_to_remove.count(hf) == 0)
            remain_faces.insert(hf);
    std::vector<std::unordered_set<typename Polyhedron::Facet_handle>> conn_comp_remain_faces;
    while(!remain_faces.empty())
    {
        std::unordered_set<typename Polyhedron::Facet_handle> component;
        std::queue<typename Polyhedron::Facet_handle> q;
        q.push(*remain_faces.begin());
        remain_faces.erase(q.front());

        while(!q.empty())
        {
            auto hf = q.front();
            q.pop();
            component.insert(hf);
            for(auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
            {
                if(nei != nullptr && remain_faces.count(nei) != 0)
                {
                    q.push(nei);
                    remain_faces.erase(nei);
                }
            }
        }

        if(component.size() <= 10)
            continue;
        
        bool component_contain_tooth = false;
        for(auto hf : component)
        {
            if(hf == nullptr)
                continue;
            for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
            {
                if(hv->_label >= 11 && hv->_label <= 49)
                {
                    component_contain_tooth = true;
                    break;
                }
            }
            if(component_contain_tooth)
                break;
        }
        if(!component_contain_tooth)
        {
            continue;
        }
        conn_comp_remain_faces.push_back(std::move(component));
    }
    std::sort(conn_comp_remain_faces.begin(), conn_comp_remain_faces.end(), [&mesh](auto comp0, auto comp1){
        int maxl0 = 0;
        int minl0 = 100;
        int maxl1 = 0;
        int minl1 = 100;
        for(auto hf : comp0)
        {
            for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
            {
                if(hv->_label > 10 && hv->_label < 50)
                {
                    maxl0 = std::max(maxl0, hv->_label);
                    minl0 = std::min(minl0, hv->_label);
                }
            }
        }
        for(auto hf : comp1)
        {
            for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
            {
                if(hv->_label > 10 && hv->_label < 50)
                {
                    maxl1 = std::max(maxl1, hv->_label);
                    minl1 = std::min(minl1, hv->_label);
                }
            }
        }
        if(maxl0 / 10 == 1 || maxl0 / 10 == 3)
            maxl0 = 37 - maxl0;
        if(minl0 / 10 == 1 || minl0 / 10 == 3)
            minl0 = 37 - minl0;
        if(maxl1 / 10 == 1 || maxl1 / 10 == 3)
            maxl1 = 37 - maxl1;
        if(minl1 / 10 == 1 || minl1 / 10 == 3)
            minl1 = 37 - minl1;
        return minl0 < minl1;
    });
    std::cout << "Found " << conn_comp_remain_faces.size() << " remaining parts." << std::endl;
    typedef CGAL::Search_traits_3<typename Polyhedron::Traits::Kernel>      Traits_base;
    typedef boost::property_map<Polyhedron,CGAL::vertex_point_t>::type            Vertex_point_pmap;
    typedef CGAL::Search_traits_adapter<typename Polyhedron::Vertex_handle, Vertex_point_pmap, Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
    typedef K_neighbor_search::Tree                                         Tree;
    typedef Tree::Splitter                                                  Splitter;
    typedef K_neighbor_search::Distance                                     Distance;
    typedef CGAL::Surface_mesh_shortest_path_traits<typename Polyhedron::Traits::Kernel, Polyhedron>    ShortestPathTraits;
    typedef boost::property_map<Polyhedron, boost::vertex_external_index_t>::const_type   Vertex_index_map;
    typedef boost::property_map<Polyhedron, CGAL::halfedge_external_index_t>::const_type  Halfedge_index_map;
    typedef boost::property_map<Polyhedron, CGAL::face_external_index_t>::const_type      Face_index_map;
    typedef CGAL::Surface_mesh_shortest_path<ShortestPathTraits,
                                         Vertex_index_map,
                                         Halfedge_index_map,
                                         Face_index_map>  Surface_mesh_shortest_path;
    struct Sequence_collector
    {
        // Luckily we don't need an accurate sequence, so we just record all 'related' vertices.
        std::unordered_set<typename Polyhedron::Facet_handle> faces;
        std::unordered_set<typename Polyhedron::Vertex_handle> vertices;
        Polyhedron& m;
        Sequence_collector(Polyhedron& _m) : m(_m) {}
        void operator()(typename Polyhedron::Halfedge_handle he, double alpha)
        {
            vertices.insert(he->vertex());
            vertices.insert(he->prev()->vertex());
            for(auto nei : CGAL::faces_around_target(he, m))
            {
                if(nei != nullptr)
                    faces.insert(nei);
            }
            for(auto nei : CGAL::faces_around_target(he->prev(), m))
            {
                if(nei != nullptr)
                    faces.insert(nei);
            }
        }
        void operator()(typename Polyhedron::Vertex_handle v)
        {
            for(auto nei : CGAL::faces_around_target(v->halfedge(), m))
            {
                if(nei != nullptr)
                    faces.insert(nei);
            }
            vertices.insert(v);
        }
        void operator()(typename Polyhedron::Facet_handle f, ShortestPathTraits::Barycentric_coordinates alpha)
        {
            if(f != nullptr)
            {
                faces.insert(f);
                vertices.insert(f->halfedge()->vertex());
                vertices.insert(f->halfedge()->next()->vertex());
                vertices.insert(f->halfedge()->prev()->vertex());
            }   
        }
    };
    
    for(int i = 1; i < conn_comp_remain_faces.size(); i++)
    {
        std::cout << "computing nearest parts..";
        auto& part1 = conn_comp_remain_faces[0];
        auto& part2 = conn_comp_remain_faces[i];
        std::unordered_set<typename Polyhedron::Vertex_handle> part1_vertices;
        std::unordered_set<typename Polyhedron::Vertex_handle> part2_vertices;
        for(auto hf : part1)
        {
            if(hf != nullptr)
            {
                for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
                    part1_vertices.insert(hv);
            }
        }
        for(auto hf : part2)
        {
            if(hf != nullptr)
            {
                for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh)) 
                    part2_vertices.insert(hv);
            }
        }
        
        std::cout << "searching nearest pair..";
        Vertex_point_pmap vppmap = get(CGAL::vertex_point, mesh);
        // Insert vertices in the tree
        Tree tree(part1_vertices.begin(), part1_vertices.end(), Splitter(), Traits(vppmap));
        std::pair<typename Polyhedron::Vertex_handle, typename Polyhedron::Vertex_handle> nearest_pair;
        double nearest_dist = std::numeric_limits<double>::max();
        for(auto hv : part2_vertices)
        {
            K_neighbor_search search(tree, hv->point(), 1, 0.0, true, Distance(vppmap));
            if(search.begin()->second < nearest_dist)
            {
                nearest_pair = {search.begin()->first, hv};
                nearest_dist = search.begin()->second;
            }
        }

        std::cout << "computing shortest path..";
        Surface_mesh_shortest_path shortest_paths(mesh,
                                            get(boost::vertex_external_index, mesh),
                                            get(CGAL::halfedge_external_index, mesh),
                                            get(CGAL::face_external_index, mesh),
                                            get(CGAL::vertex_point, mesh));
        std::cout << "query..";
        // TODO: this accurate sequence is very slow. maybe we can 1. use a sub mesh 2. try other methods e.g. Dijkstra?
        shortest_paths.add_source_point(nearest_pair.first);
        Sequence_collector sequence(mesh);
        auto [dist, it] = shortest_paths.shortest_path_sequence_to_source_points(nearest_pair.second, sequence);
        std::cout << "query end..";

        if(dist < 0)
        {
            std::cout << "ERROR: components are not reachable" << std::endl;
        }
        // make the connecting part wider
        std::unordered_set<typename Polyhedron::Vertex_handle> connector_vertices;
        std::queue<typename Polyhedron::Vertex_handle> queue;
        for(auto hv : sequence.vertices)
        {
            queue.push(hv);
        }
        while(!queue.empty())
        {
            auto hv = queue.front();
            queue.pop();
            connector_vertices.insert(hv);
            for(auto nei : CGAL::vertices_around_target(hv, mesh))
            {
                if(connector_vertices.count(nei) != 0 || nei->_label > 10)
                {
                    continue;
                }
                double min_dist = std::numeric_limits<double>::max();
                for(auto vv : sequence.vertices)
                {
                    min_dist = std::min(min_dist, CGAL::squared_distance(vv->point(), nei->point()));
                }
                if(min_dist < 4.0)
                {
                    queue.push(nei);
                }
            }
        }

        std::unordered_set<typename Polyhedron::Facet_handle> connector_faces;
        for(auto hv : connector_vertices)
        {
            for(auto hf : CGAL::faces_around_target(hv->halfedge(), mesh))
            {
                connector_faces.insert(hf);
            }
        }
        for(auto hf : sequence.faces)
        {
            connector_faces.insert(hf);
            for(auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                connector_faces.insert(nei);
        }

        for(auto hf : part2)
            part1.insert(hf);
        for(auto hf : connector_faces)
            part1.insert(hf);
        std::cout << "done." << std::endl;
    }
    // conn_comp_remain_faces[0] should contains all faces now
    std::erase_if(face_to_remove, [&conn_comp_remain_faces](auto hf) { return conn_comp_remain_faces[0].count(hf) != 0; });

#ifdef DEBUG_ORTHOSCANBASE
    {
        std::vector<Polyhedron::Point_3> points;
        std::vector<Eigen::Vector3f> colors;
        std::vector<TTriangle<size_t>> faces;
        int nv = 0;
        int nt = 0;
        for(auto hf : CGAL::faces(mesh))
        {
            auto p0 = hf->halfedge()->vertex()->point();
            auto p1 = hf->halfedge()->next()->vertex()->point();
            auto p2 = hf->halfedge()->prev()->vertex()->point();
            points.push_back(p0);
            points.push_back(p1);
            points.push_back(p2);
            faces.emplace_back(nv, nv + 1, nv + 2);
            nv += 3;
            if(face_to_remove.count(hf) != 0)
            {
                colors.push_back(Eigen::Vector3f(1, 0, 0));
                colors.push_back(Eigen::Vector3f(1, 0, 0));
                colors.push_back(Eigen::Vector3f(1, 0, 0));
            }
            else
            {
                colors.push_back(Eigen::Vector3f(0.7f, 0.7f, 0.7f));
                colors.push_back(Eigen::Vector3f(0.7f, 0.7f, 0.7f));
                colors.push_back(Eigen::Vector3f(0.7f, 0.7f, 0.7f));
            }
        }
        WriteVCFAssimp<Polyhedron::Traits::Kernel, size_t>("removed2.obj", points, colors, faces);
    }
#endif

    for (auto hf : face_to_remove)
    {
        if(hf != nullptr && hf->halfedge() != nullptr && !hf->halfedge()->is_border())
        {
            mesh.erase_facet(hf->halfedge());
        }
    }
#ifdef DEBUG_ORTHOSCANBASE
    mesh.WriteOBJ("cleaned.obj");
#endif
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
    CGAL::Polygon_mesh_processing::fair(mesh, std::vector<Polyhedron::Vertex_handle>(vertex_to_fair.begin(), vertex_to_fair.end()));

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
    CGAL::Polygon_mesh_processing::triangulate_hole(mesh, hole_hh,
        CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)).vertex_output_iterator(std::back_inserter(patch_vertex)));
    for (auto hv : patch_vertex)
    {
        hv->_label = 1;
    }

    borders.clear();
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(borders));
    for (auto hh : borders)
    {
        CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(mesh, hh);
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
        int l0 = hf->halfedge()->vertex()->_label;
        int l1 = hf->halfedge()->next()->vertex()->_label;
        int l2 = hf->halfedge()->next()->next()->vertex()->_label;
        if(l0 == 0 || l0 % 10 == 8) cnt++;
        if(l1 == 0 || l1 % 10 == 8) cnt++;
        if(l2 == 0 || l2 % 10 == 8) cnt++;
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
    mesh.UpdateFaceLabels();
}

std::vector<std::array<std::pair<int, float>, 3>> ComputeDeformWeights(const Polyhedron& gum_mesh, const CrownFrames<KernelEpick>& crown_frames, bool upper)
{
    int label_start = upper ? 11 : 31;
    int label_end = upper ? 28 : 48;
    std::vector<std::array<std::pair<int, float>, 3>> weight_table;
    for(auto hv : CGAL::vertices(gum_mesh))
    {
        auto p = hv->point();
        std::vector<std::pair<int, float>> weights;
        for(int label = label_start; label <= label_end; label++)
        {
            if(crown_frames.Frames().count(label) == 0)
                continue;
            auto frame_center = crown_frames.GetFrame(label).pos;
            auto frame_dir = crown_frames.GetFrame(label).up;
            double dist_to_frame = CGAL::squared_distance( KernelEpick::Line_3{ frame_center, frame_dir }, p );
            weights.emplace_back(label, dist_to_frame);
        }

        std::sort(weights.begin(), weights.end(), [](auto& lh, auto& rh){ return lh.second < rh.second; });
        std::array<std::pair<int, float>, 3> weight_table_this;
        weight_table_this.fill({0, 0.f});
        if(weights[0].second < 1e-8f)
        {
            weight_table_this[0] = {weights[0].first, 1.0f};
        }
        else
        {
            float w_sum = 0.f;
            for(int i = 0; i < 3 && i < weights.size(); i++)
            {
                float w = weights[i].second;
                w = std::exp(-w * w / 25.0);
                weight_table_this[i] = {weights[i].first, w};
                w_sum += weight_table_this[i].second;
            }
            for(auto& [label, w] : weight_table_this)
                if(label != 0)
                    w /= w_sum;
        }
        weight_table.push_back(weight_table_this);
    }

    return weight_table;
}

std::vector<std::array<std::pair<int, float>, 3>> ComputeDeformWeightsGeoDist(const Polyhedron& gum_mesh, const CrownFrames<KernelEpick>& crown_frames, bool upper)
{
    std::array<bool, 50> has_label;
    std::fill(has_label.begin(), has_label.end(), false);
    for(auto hv : CGAL::vertices(gum_mesh))
    {
        if(hv->_label > 10 || hv->_label == 1)
        {
            has_label[hv->_label] = true;
        }
    }
    std::vector<int> all_labels;
    for(int label = 1; label < 50; label++)
    {
        if(has_label[label])
        {
            all_labels.push_back(label);
        }
    }

    using LA = CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<typename CGAL::Eigen_sparse_matrix<double>::EigenType>>;
    using HeatMethod = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Polyhedron, CGAL::Heat_method_3::Direct,
                                                                              typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::const_type, LA, Polyhedron::Traits::Kernel>;    
    std::unordered_map<int, std::unordered_map<typename Polyhedron::Vertex_handle, double>> dist_to_label;
#pragma omp parallel for
    for(int i = 0; i < all_labels.size(); i++)
    {
        int label = all_labels[i];
        std::unordered_map<typename Polyhedron::Vertex_handle, double> distances;
        HeatMethod heat_method(gum_mesh);
        for(auto hv : CGAL::vertices(gum_mesh))
        {
            if(hv->_label == label)
                heat_method.add_source(hv);
        }
        heat_method.estimate_geodesic_distances(boost::make_assoc_property_map(distances));
#pragma omp critical
{
        std::cout << "geo dist " << label << std::endl;
        dist_to_label.insert({label, std::move(distances)});
}
    }

    std::vector<std::array<std::pair<int, float>, 3>> weights_table;
    for(auto hv : CGAL::vertices(gum_mesh))
    {
        std::array<std::pair<int, float>, 3> v_weights;
        std::fill(v_weights.begin(), v_weights.end(), std::pair<int, float>(0, 0.f));
        if(hv->_label == 1)
        {
            v_weights[0] = {1, 1.f};
            v_weights[1] = {1, 0.f};
            v_weights[2] = {1, 0.f};
        }
        else if(hv->_label > 10)
        {
            v_weights[0] = {hv->_label, 1.f};
            v_weights[1] = {hv->_label, 0.f};
            v_weights[2] = {hv->_label, 0.f};
        }
        else
        {
            std::vector<std::pair<int, float>> label_dist_pairs;
            for(int label = 1; label < 50; label++)
            {
                if(dist_to_label.count(label) == 0 || dist_to_label.at(label).count(hv) == 0)
                    continue;
                label_dist_pairs.push_back({label, static_cast<float>(dist_to_label.at(label).at(hv))});
            }
            std::sort(label_dist_pairs.begin(), label_dist_pairs.end(), [](auto& lh, auto& rh){ return lh.second < rh.second; });
            label_dist_pairs.resize(v_weights.size());
            float w_sum = 0.f;
            for(int i = 0; i < v_weights.size(); i++)
            {
                v_weights[i].first = label_dist_pairs[i].first;
                v_weights[i].second = std::exp(-label_dist_pairs[i].second * label_dist_pairs[i].second / 10.f);
                w_sum += v_weights[i].second;
            }
            for(auto& [_, w] : v_weights)
                w /= w_sum;
        }
        weights_table.push_back(v_weights);
    }
//#define DEBUG_WEIGHTS
#ifdef DEBUG_WEIGHTS
    std::ofstream ofs("./weighted_v.obj");
    int i = 0;
    for(auto hv : CGAL::vertices(gum_mesh))
    {
        auto& weights = weights_table[i];
        auto p = hv->point();
        Eigen::Vector3f c = LabelColorMap(weights[0].first] * weights[0].second;
        c += LabelColorMap(weights[1].first] * weights[1].second;
        c += LabelColorMap(weights[2].first] * weights[2].second;
        ofs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << ' ' << c.x() << ' ' << c.y() << ' ' << c.z() << '\n';
        i++;
    }
#endif
    return weights_table;
}

void GenerateGum(std::string output_gum, std::string crown_frame, bool upper, Polyhedron& mesh)
{
    std::cout << "Creating gum..." << std::endl;
    // Here we remove tooth mesh one by one, in this way it creates smaller hole which are easier to close.
    int start_label = upper ? 11 : 31;
    int end_label = upper ? 29 : 49;
    for(int label = start_label; label < end_label; label++)
    {
        std::unordered_set<Polyhedron::Facet_handle> face_to_remove;
        for(auto hf : CGAL::faces(mesh))
        {
            for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
            {
                if(hv->_label == label)
                {
                    face_to_remove.insert(hf);
                    break;
                }
            }
        }
        if(face_to_remove.empty())
        {
            continue;
        }

        for(auto hf : face_to_remove)
        {
            if(hf != nullptr && hf->halfedge() != nullptr && !hf->halfedge()->is_border())
            {
                mesh.erase_facet(hf->halfedge());
            }
        }

        auto [V, F] = mesh.ToVerticesTriangles();
        std::vector<std::pair<std::vector<Polyhedron::Vertex_handle>, std::vector<Polyhedron::Facet_handle>>> patches;
        FixMeshWithLabel(V, F, mesh.WriteLabels(), mesh, true, 999, false, false, 99999999, 99999999.0f, true, 10, &patches);
        for(auto& [patch_vertices, patch_faces] : patches)
        {
            for(auto hf : patch_faces)
                hf->_label = 2;
            for(auto hv : patch_vertices)
                hv->_label = 2;
        }
    }
    std::vector<typename Polyhedron::Vertex_handle> vertex_to_fair;
    for(auto hv : CGAL::vertices(mesh))
        if(hv->_label == 2)
            vertex_to_fair.push_back(hv);

    CGAL::Polygon_mesh_processing::fair(mesh, vertex_to_fair);

    std::unordered_set<typename Polyhedron::Facet_handle> filter_faces;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label == 2)
            for(auto nei : CGAL::faces_around_target(hv->halfedge(), mesh))
                filter_faces.insert(nei);
    }
    CGAL::Face_filtered_graph<Polyhedron> filtered_graph(mesh, filter_faces);
    Polyhedron filtered_mesh;
    std::unordered_map<typename Polyhedron::Vertex_handle, typename Polyhedron::Vertex_handle> v_to_v_map;
    std::unordered_map<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle> f_to_f_map;
    CGAL::copy_face_graph(filtered_graph, filtered_mesh, 
        CGAL::parameters::vertex_to_vertex_map(boost::make_assoc_property_map(v_to_v_map)).face_to_face_map(boost::make_assoc_property_map(f_to_f_map)));
    for(auto& [src, tar] : v_to_v_map)
        tar->_label = src->_label;
    for(auto& [src, tar] : f_to_f_map)
        tar->_label = src->_label;

    using LA = CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<typename CGAL::Eigen_sparse_matrix<double>::EigenType>>;
    using HeatMethod = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Polyhedron, CGAL::Heat_method_3::Direct,
                                                                            typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::const_type, LA, Polyhedron::Traits::Kernel>;    
    std::unordered_map<typename Polyhedron::Vertex_handle, double> distances;
    HeatMethod heat_method(filtered_mesh);
    for(auto hv : CGAL::vertices(filtered_mesh))
        if(hv->_label != 2)
            heat_method.add_source(hv);
    heat_method.estimate_geodesic_distances(boost::make_assoc_property_map(distances));

    CrownFrames<KernelEpick> crown_frames(crown_frame);
    std::vector<typename Polyhedron::Vertex_handle> roi_vertices;
    for(auto hv : CGAL::vertices(mesh))
        if(hv->_label == 2)
            roi_vertices.push_back(hv);

    std::unordered_map<typename Polyhedron::Vertex_handle, typename Polyhedron::Point_3> control_vertices_map;
    std::vector<typename Polyhedron::Vertex_handle> control_vertices;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label == 2)
        {
            double min_dist = std::numeric_limits<double>::max();
            std::vector<std::pair<int, double>> frame_weights;
            for(const auto& frame : crown_frames.Frames())
            {
                if(upper == (frame.first >= 11 && frame.first <= 29))
                {
                    double dist = std::sqrt(CGAL::squared_distance(KernelEpick::Line_3(frame.second.pos, frame.second.up), hv->point()));
                    if(dist == 0.0)
                        dist = 1e-9; // avoid divide by 0
                    frame_weights.push_back({frame.first, dist});
                }
            }
            std::sort(frame_weights.begin(), frame_weights.end(), [](auto& l, auto& r){ return l.second < r.second; });
            if(frame_weights.size() > 1 && std::abs(frame_weights[0].second - frame_weights[1].second) < 0.8)
            {
                hv->_label = 3;
                continue;
            }
            auto dir = crown_frames.GetFrame(frame_weights.front().first).up;
            if(frame_weights.size() > 1)
            {
                double w0 = std::exp(-frame_weights[0].second * frame_weights[0].second);
                double w1 = std::exp(-frame_weights[1].second * frame_weights[1].second);

                dir = dir * w0 + crown_frames.GetFrame(frame_weights[1].first).up * w1;
                dir /= (w0 + w1);
                dir /= std::sqrt(dir.squared_length());
            }
            if(distances[v_to_v_map[hv]] > 0.0)
            {
                double x = distances[v_to_v_map[hv]];
                double dist = 2 * (1.0 - std::exp(-x * x / 4.0));
                control_vertices_map[hv] = hv->point() - dir * dist;
                control_vertices.push_back(hv);
                hv->_label = frame_weights[0].first;
            }
        }
    }
    
    CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default> deformation(mesh);
    deformation.insert_roi_vertices(roi_vertices.begin(), roi_vertices.end());
    deformation.insert_control_vertices(control_vertices.begin(), control_vertices.end());
    for(auto& [hv, pos] : control_vertices_map)
    {
        deformation.set_target_position(hv, pos);
    }
    deformation.deform(10, 1e-4);

    mesh.WriteOBJ(output_gum);
    mesh.WriteLabels(output_gum.substr(0, output_gum.rfind('.')) + ".json");

    std::cout << "Compute weights..." << std::endl;
    auto weights = ComputeDeformWeightsGeoDist(mesh, crown_frames, upper);
    {
        std::ofstream ofs;
        std::string weight_name;
        if( output_gum.rfind('.') == std::string::npos)
            ofs.open(output_gum + "_w.json");
        else 
            ofs.open(output_gum.substr(0, output_gum.rfind('.')) + "_w.json");
        nlohmann::json json;
        json["weights"] = weights;
        ofs << json;
    }

    std::cout << "Done." << std::endl;
}

#ifndef FOUND_PYBIND11
int main(int argc, char *argv[])
{
    argparse::ArgumentParser parser;
    parser.add_argument("--input_file", "-i").required().help("specify the input mesh.");
    parser.add_argument("--input_label", "-l").required().help("specify the input labels.");
    parser.add_argument("--output_file", "-o").help("specify the output file.");
    parser.add_argument("--output_label", "-ol").help("specify the output labels.");
    parser.add_argument("--output_gum", "-g").help("specify the output gum mesh file.");
    parser.add_argument("--crown_frame", "-c").help("specify the crown frame file.");
    try
    {
        parser.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }

    Polyhedron mesh;
    std::string input_file = parser.get("-i");
    std::string label_file = parser.get("-l");
    if (!CGAL::IO::read_polygon_mesh(input_file, mesh, CGAL::parameters::verbose(true)))
    {
        try
        {
            printf("possible invalid mesh, try fixing...\n");
            std::vector<typename Polyhedron::Traits::Point_3> vertices;
            std::vector<TTriangle<size_t>> faces;
            if(input_file.ends_with(".obj"))
            {
                LoadVFObj<typename Polyhedron::Traits, size_t>(input_file, vertices, faces);
            }
            else
            {
                // TODO: add custom file loading function for more formats, since assimp cannot keep vertex order sometimes.
                LoadVFAssimp<typename Polyhedron::Traits, size_t>(input_file, vertices, faces);
            }
            std::vector<int> labels = LoadLabels(label_file);
            FixMeshWithLabel(vertices, faces, labels, mesh, true, 9999, false, false, 0, 0, false, 10);
        }
        catch(const std::exception&)
        {
            throw IOError("Cannot read mesh file or mesh invalid: " + input_file);
        }
    }
    else
    {
        mesh.LoadLabels(label_file);
        auto [v, f] = mesh.ToVerticesTriangles();
        FixMeshWithLabel(v, f, mesh.WriteLabels(), mesh, true, 1000, false, false, 0, 0, false, 10);
        //mesh.UpdateFaceLabels2();
    }
    LabelProcessing(mesh);
    SegClean(mesh);

    bool upper = true;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label >= 11 && hv->_label <= 29)
            break;
        else if(hv->_label >= 31 && hv->_label <= 49)
        {
            upper = false;
            break;
        }
    }
    try
    {
        std::cout << "Optimizing...";
        Optimize(mesh);
        auto [v, f] = mesh.ToVerticesTriangles();
        FixMeshWithLabel(v, f, mesh.WriteLabels(), mesh, true, 1000, false, false, 0, 0, false, 10);
        std::cout << "Done." << std::endl;
#ifdef DEBUG_ORTHOSCANBASE
        mesh.WriteOBJ("optimized.obj");
#endif
        std::cout << "Generating...";
        GenerateBase2(mesh);
        std::cout << "Done." << std::endl;
        if(parser.present("-o").has_value()) {
            mesh.WriteOBJ(parser.get("-o"));
        }
        if(parser.present("-ol").has_value()) {
            mesh.WriteLabels(parser.get("-ol"));
        }

        if(parser.present("-g").has_value() && parser.present("-c").has_value())
        {
            GenerateGum(parser.present("-g").value(), parser.present("-c").value(), upper, mesh);
        }
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    return 0;
}

#else
int GenerateGumApi(std::string input_file, std::string input_label, std::string crown_frame, std::string output_gum)
{
    Polyhedron mesh;
    if (!CGAL::IO::read_polygon_mesh(input_file, mesh, CGAL::parameters::verbose(true)))
    {
        try
        {
            printf("possible invalid mesh, try fixing...\n");
            std::vector<typename Polyhedron::Traits::Point_3> vertices;
            std::vector<TTriangle<size_t>> faces;
            if(input_file.ends_with(".obj"))
            {
                LoadVFObj<typename Polyhedron::Traits, size_t>(input_file, vertices, faces);
            }
            else
            {
                // TODO: add custom file loading function for more formats, since assimp cannot keep vertex order sometimes.
                LoadVFAssimp<typename Polyhedron::Traits, size_t>(input_file, vertices, faces);
            }
            std::vector<int> labels = LoadLabels(input_label);
            FixMeshWithLabel(vertices, faces, labels, mesh, true, 9999, false, false, 0, 0, false, 10);
        }
        catch(const std::exception&)
        {
            throw IOError("Cannot read mesh file or mesh invalid: " + input_file);
        }
    }
    else
    {
        mesh.LoadLabels(input_label);
        //mesh.UpdateFaceLabels2();
    }
    LabelProcessing(mesh);
    SegClean(mesh);

    bool upper = true;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label >= 11 && hv->_label <= 29)
            break;
        else if(hv->_label >= 31 && hv->_label <= 49)
        {
            upper = false;
            break;
        }
    }
    try
    {
        std::cout << "Optimizing...";
        Optimize(mesh);
        std::cout << "Done." << std::endl;
        std::cout << "Generating...";
        GenerateBase2(mesh);
        std::cout << "Done." << std::endl;

        GenerateGum(output_gum, crown_frame, upper, mesh);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    return 0;
}

#endif