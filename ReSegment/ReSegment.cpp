#include <iostream>
#include <format>
#include <stack>
#include <string>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/version.h>
#include <nlohmann/json.hpp>
#include "../Polyhedron.h"

#define DEBUG_OUTPUT
namespace
{
    using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;
    using AABBPrimitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using AABBTraits = CGAL::AABB_traits<KernelEpick, AABBPrimitive>;
    using AABBTree = CGAL::AABB_tree<AABBTraits>;
    using hHalfedge = Polyhedron::Halfedge_handle;
    using hVertex = Polyhedron::Vertex_handle;
    using hFacet = Polyhedron::Facet_handle;
    using Halfedge = Polyhedron::Halfedge;
    using CVertex = Polyhedron::Vertex;
    using Facet = Polyhedron::Facet;
    using iHalfedge = Polyhedron::Halfedge_iterator;
    using iVertex = Polyhedron::Vertex_iterator;
    using iFacet = Polyhedron::Facet_iterator;
    using Point_3 = Polyhedron::Point_3;
    using Vec3 = Polyhedron::Traits::Vector_3;
    using Triangle = Polyhedron::Triangle;
    using Edge = Polyhedron::Edge;

    std::vector<hVertex> ConnectedComponent(hVertex hv, const Polyhedron &mesh, std::unordered_set<hVertex> &processed_set)
    {
        std::unordered_set<hVertex> vertices;

        vertices.insert(hv);
        processed_set.insert(hv);

        bool found = true;
        while (found)
        {
            found = false;
            auto temp = vertices;
            for (auto v : temp)
            {
                for (auto nei : CGAL::vertices_around_target(v, mesh))
                {
                    if (!processed_set.contains(nei))
                    {
                        vertices.insert(nei);
                        processed_set.insert(nei);
                        found = true;
                    }
                }
            }
        }

        return std::vector<hVertex>(vertices.begin(), vertices.end());
    }

    std::vector<Point_3> Resample(const std::vector<Point_3> &split_points, double threshold)
    {
        std::vector<Point_3> new_points;
        for (int i = 0; i < split_points.size(); i++)
        {
            Point_3 p0 = split_points[i];
            Point_3 p1 = split_points[(i + 1) % split_points.size()];
            new_points.push_back(p0);
            double dist2 = CGAL::squared_distance(p0, p1);
            if (dist2 > threshold * threshold)
            {
                double dist = std::sqrt(dist2);
                int nb_insert = static_cast<int>(dist / threshold);
                for (int j = 1; j < nb_insert; j++)
                {
                    new_points.push_back(p0 + (p1 - p0) * j / (nb_insert + 1));
                }
            }
            new_points.push_back(p1);
        }

        return new_points;
    }

    template <typename It>
    bool IsConnectedComponent(It begin, It end, const Polyhedron &mesh)
    {
        std::unordered_set<hFacet> all_faces(begin, end);
        hFacet start = *begin;
        std::unordered_set<hFacet> search_faces;
        search_faces.insert(start);
        bool found = true;
        while (found)
        {
            found = false;
            auto temp = search_faces;
            for (hFacet hf : temp)
            {
                for (hFacet nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                {
                    if (all_faces.count(nei) != 0 && search_faces.count(nei) == 0)
                    {
                        search_faces.insert(nei);
                        found = true;
                    }
                }
            }
        }

        return search_faces.size() == all_faces.size();
    }

    std::vector<KernelEpick::Segment_3> Subdivide(const KernelEpick::Segment_3 &seg, int num)
    {
        std::vector<KernelEpick::Segment_3> segments;
        auto p = seg.point(0);
        auto d = seg.point(1) - seg.point(0);
        for (int i = 0; i < num; i++)
        {
            segments.emplace_back(p + d * (double)(i) / num, p + d * (double)(i + 1) / num);
        }
        return segments;
    }

    std::vector<KernelEpick::Segment_3> Subdivide(const std::vector<KernelEpick::Segment_3> &seg, int num)
    {
        std::vector<KernelEpick::Segment_3> segments;
        for (const auto &s : seg)
        {
            auto cur = Subdivide(s, num);
            segments.insert(segments.end(), cur.begin(), cur.end());
        }
        return segments;
    }

    bool Backtracking(hFacet to, const Polyhedron &mesh, std::vector<hFacet> &path, std::unordered_set<hFacet> &processed_faces, size_t max_len)
    {
        if (path.back() == to)
            return true;
        if (path.size() >= max_len)
        {
            return false;
        }

        for (auto nei : CGAL::faces_around_face(path.back()->halfedge(), mesh))
        {
            if (processed_faces.contains(nei))
                continue;
            processed_faces.insert(nei);
            path.push_back(nei);
            if (Backtracking(to, mesh, path, processed_faces, max_len))
                return true;
            path.pop_back();
            processed_faces.erase(nei);
        }

        return false;
    }

    std::vector<hFacet> FindPath(hFacet from, hFacet to, const Polyhedron &mesh)
    {
        if (from == to)
            return std::vector<hFacet>{from};
        int dfs_range = 0;
        std::unordered_set<hFacet> processed_faces;
        std::queue<hFacet> q;
        q.push(from);
        processed_faces.insert(from);
        bool found = false;
        while (!q.empty())
        {
            size_t layer_size = q.size();
            for (int i = 0; i < layer_size; i++)
            {
                hFacet front_facet = q.front();
                q.pop();
                for (auto nei : CGAL::faces_around_face(front_facet->halfedge(), mesh))
                {
                    if(nei == nullptr || processed_faces.contains(nei))
                    {
                        continue;
                    }
                    q.push(nei);
                    processed_faces.insert(nei);
                    if (nei == to)
                        found = true;
                }
            }
            dfs_range++;
            if (found)
            {
                break;
            }
        }

        Polyhedron submesh;
        CGAL::Face_filtered_graph<Polyhedron> filtered_graph(mesh, processed_faces);
        std::unordered_map<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle> f_map;
        CGAL::copy_face_graph(filtered_graph, submesh, CGAL::parameters::face_to_face_map(boost::make_assoc_property_map(f_map)));
        CGAL::set_halfedgeds_items_id(submesh);
        std::unordered_map<typename Polyhedron::Facet_handle, typename Polyhedron::Facet_handle> inv_f_map;
        for(auto [src, tar] : f_map)
        {
            inv_f_map.insert({tar, src});
        }
        using SPTraits = CGAL::Surface_mesh_shortest_path_traits<typename Polyhedron::Traits::Kernel, Polyhedron>;
        struct Sequence_collector
        {
            std::vector<typename Polyhedron::Facet_handle> sequence;
            Polyhedron& mesh;

            Sequence_collector(Polyhedron& _mesh) : mesh(_mesh) {};

            void operator()(typename Polyhedron::Halfedge_handle he, double alpha)
            {
                if(he->facet() != nullptr)
                    sequence.push_back(he->facet());
                if(he->opposite()->facet() != nullptr)
                    sequence.push_back(he->opposite()->facet());
            }
            void operator()(typename Polyhedron::Vertex_handle v)
            {
                for(auto nei : CGAL::faces_around_target(v->halfedge(), mesh))
                {
                    if(nei != nullptr)
                        sequence.push_back(nei);
                }
            }
            void operator()(typename Polyhedron::Facet_handle f, SPTraits::Barycentric_coordinates alpha)
            {
                sequence.push_back(f);
            }
        };
        CGAL::Surface_mesh_shortest_path<SPTraits> shortest_path(submesh);
        shortest_path.add_source_point(f_map.at(from), SPTraits::Barycentric_coordinate({0.33, 0.33, 0.34}));
        Sequence_collector sc(submesh);
        shortest_path.shortest_path_sequence_to_source_points(f_map.at(to), SPTraits::Barycentric_coordinate({0.33, 0.33, 0.34}), sc);

        std::vector<hFacet> path;
        for(auto hf : sc.sequence)
        {
            path.push_back(inv_f_map[hf]);
        }
        return path;
    }

}

void ReSegmentOneLabel(const Polyhedron &mesh, const AABBTree &aabb_tree, const std::vector<std::vector<Point_3>> &split_points_list,
                        std::vector<int> &output_labels, int label, double intersection_width, int cutface_orit_smooth)
{
    std::unordered_set<hFacet> all_intersect_faces;
    for (auto &split_points : split_points_list)
    {
        std::vector<KernelEpick::Point_3> project_points;
        std::vector<KernelEpick::Vector_3> normals;
        std::vector<hFacet> project_facets;
        for (auto &pt : split_points)
        {
            auto [p, hf] = aabb_tree.closest_point_and_primitive(pt);
            project_points.push_back(p);
            normals.push_back(FaceNormal(hf));
            project_facets.push_back(hf);
        }

        // smooth the normals
        // for (int i = 0; i < cutface_orit_smooth; i++)
        // {
        //     std::vector<KernelEpick::Vector_3> new_normals(normals.size());
        //     for (size_t j = 0; j < normals.size(); j++)
        //     {
        //         size_t j0 = j > 0 ? (j - 1) : (normals.size() - 1);
        //         size_t j1 = j;
        //         size_t j2 = (j + 1) % normals.size();
        //         new_normals[j1] = (normals[j0] + normals[j2]) * 0.5;
        //         new_normals[j1] /= std::sqrt(new_normals[j1].squared_length());
        //     }
        //     normals = std::move(new_normals);
        // }

#ifdef DEBUG_OUTPUT
        std::vector<KernelEpick::Triangle_3> cut_faces;
#endif

        for (size_t i = 0; i < split_points.size(); i++)
        {
            size_t i0 = i;
            size_t i1 = (i + 1) % split_points.size();

            auto p0 = project_points[i0];
            auto p1 = project_points[i1];
            auto n0 = normals[i0];
            auto n1 = normals[i1];
            auto hf0 = project_facets[i0];
            auto hf1 = project_facets[i1];
            double depth = intersection_width;

            KernelEpick::Triangle_3 t0(p0 - n0 * depth, p1 - n1 * depth, p0 + n0 * depth);
            KernelEpick::Triangle_3 t1(p1 - n1 * depth, p0 + n0 * depth, p1 + n1 * depth);

            std::vector<hFacet> this_intersect_faces;
            aabb_tree.all_intersected_primitives(t0, std::back_inserter(this_intersect_faces));
            aabb_tree.all_intersected_primitives(t1, std::back_inserter(this_intersect_faces));

            if (!IsConnectedComponent(this_intersect_faces.begin(), this_intersect_faces.end(), mesh))
            {
                auto path = FindPath(hf0, hf1, mesh);
                for (auto hf : path)
                {
                    all_intersect_faces.insert(hf);
                }
#ifdef DEBUG_OUTPUT
                cut_faces.push_back(t0);
                cut_faces.push_back(t1);
#endif
            }
            else
            {
                all_intersect_faces.insert(this_intersect_faces.begin(), this_intersect_faces.end());
#ifdef DEBUG_OUTPUT
                cut_faces.push_back(t0);
                cut_faces.push_back(t1);
#endif
            }
        }
        std::cout << "Intersection: " << all_intersect_faces.size() << std::endl;
        #ifdef DEBUG_OUTPUT
            {
                std::ofstream ofs("cutfaces.obj");
                int cnt = 1;
                for (const auto &tri : cut_faces)
                {
                    for (int i = 0; i < 3; i++)
                        ofs << std::format("v {} {} {}\n", tri.vertex(i).x(), tri.vertex(i).y(), tri.vertex(i).z());
                    ofs << std::format("f {} {} {}\n", cnt, cnt + 1, cnt + 2);
                    cnt += 3;
                }
                ofs.close();
            }
        #endif
    }

    std::unordered_set<hVertex> processed_set;
#pragma omp critical
    {
        for (auto hf : all_intersect_faces)
        {
            for (auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
            {
                processed_set.insert(hv);
                output_labels[hv->id()] = label;
            }
        }
    }

    std::vector<std::vector<hVertex>> connected_components;
    for (auto hv : CGAL::vertices(mesh))
    {
        if (!processed_set.contains(hv))
        {
            connected_components.push_back(ConnectedComponent(hv, mesh, processed_set));
        }
    }

    std::sort(connected_components.begin(), connected_components.end(), [](auto &lhs, auto &rhs)
              { return lhs.size() < rhs.size(); });
    connected_components.pop_back();
#pragma omp critical
    {
        for (auto &teeth : connected_components)
        {
            std::cout << "Tooth part: " << teeth.size() << std::endl;
            for (auto hv : teeth)
            {
                output_labels[hv->id()] = label;
            }
        }
    }
}

void ReSegmentOneLabel(const Polyhedron &mesh, const AABBTree &aabb_tree, const std::vector<std::vector<double>> &split_points,
                        std::vector<int> &output_labels, int label, double intersection_width, int cutface_orit_smooth)
{
    std::vector<std::vector<Point_3>> points;
    for (auto &line : split_points)
    {
        std::vector<Point_3> line_points;
        for (int i = 0; i < line.size() / 3; i++)
            line_points.emplace_back(line[i * 3], line[i * 3 + 1], line[i * 3 + 2]);
        points.emplace_back(std::move(line_points));
    }
    ReSegmentOneLabel(mesh, aabb_tree, points, output_labels, label, intersection_width, cutface_orit_smooth);
}

bool ReSegmentLabels(
    std::string input_mesh,
    const std::vector<std::vector<std::vector<double>>> &splitlines,
    const std::vector<int> &splitline_labels,
    std::string output_json,
    double intersection_width,
    int cutface_orit_smooth,
    bool upper)
{
    Polyhedron mesh;
    if(CGAL::IO::read_polygon_mesh(input_mesh, mesh, CGAL::parameters::verbose(true)))
    {
        printf("Load mesh: V = %zd, F = %zd\n", mesh.size_of_vertices(), mesh.size_of_facets());
    }
    else
    {
        std::cout << "Error: failed to read mesh: " << input_mesh << std::endl;
        return false;
    }
    if(!mesh.is_valid(false))
    {
        std::cout << "Error: input mesh not valid:" << std::endl;
        mesh.is_valid(true);
        return false;
    }
    if(!mesh.is_pure_triangle())
    {
        std::cout << "Error: input mesh has non triangle face." << std::endl;
        return false;
    }
    std::cout << "Setting id" << std::endl; 
    CGAL::set_halfedgeds_items_id(mesh);
    std::cout << "Constructing aabb tree" << std::endl;
    AABBTree aabb_tree(mesh.facets_begin(), mesh.facets_end(), mesh);
    if(aabb_tree.empty())
    {
        std::cout << "Error: failed to create AABB tree." << std::endl;
        return false;
    }
    std::vector<int> output_labels(mesh.size_of_vertices(), 0);

    std::vector<int> indices_to_process;
    for (int i = 0; i < splitline_labels.size(); i++)
    {
        if (upper && splitline_labels[i] < 30)
            indices_to_process.push_back(i);
        if (!upper && splitline_labels[i] > 30)
            indices_to_process.push_back(i);
    }
//#pragma omp parallel for
    for (int i = 0; i < indices_to_process.size(); i++)
    {   
        int idx = indices_to_process[i];
        std::cout << "try resegment label " << splitline_labels[idx] << std::endl;
        ReSegmentOneLabel(mesh, aabb_tree, splitlines[idx], output_labels, splitline_labels[idx], intersection_width, cutface_orit_smooth);
    }

    nlohmann::json json;
    json["labels"] = output_labels;
    std::ofstream ofs(output_json);
    if(ofs.fail())
    {
        return false;
    }
    ofs << json;
    return true;
}
