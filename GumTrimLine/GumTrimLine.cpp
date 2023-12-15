#include <deque>
#include <iostream>
#include <filesystem>
#include <unordered_map>
#include "../Polyhedron.h"
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include "../MeshFix/MeshFix.h"

namespace
{
    using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;
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

    std::vector<hVertex> ConnectedComponents(hVertex hv, Polyhedron &mesh)
    {
        std::vector<hVertex> vertices;
        std::deque<hVertex> front;
        front.push_back(hv);
        vertices.push_back(hv);
        hv->_processed = true;

        while (!front.empty())
        {
            hVertex v = front.front();
            front.pop_front();
            for (auto nei : CGAL::vertices_around_target(v, mesh))
            {
                if (!nei->_processed && nei->_label != 0)
                {
                    front.push_back(nei);
                    vertices.push_back(nei);
                    nei->_processed = true;
                }
            }
        }

        return vertices;
    }

    std::vector<hFacet> FacialConnectedComponents(hFacet hf, Polyhedron &mesh)
    {
        std::vector<hFacet> faces;
        std::deque<hFacet> front;

        front.push_back(hf);
        faces.push_back(hf);

        hf->_processed = true;

        while (!front.empty())
        {
            hFacet f = front.front();
            front.pop_front();
            for (auto nei : CGAL::faces_around_face(f->halfedge(), mesh))
            {
                if (nei != nullptr && !nei->_processed && nei->_label != 0)
                {
                    front.push_back(nei);
                    faces.push_back(nei);
                    nei->_processed = true;
                }
            }
        }

        return faces;
    }

    std::vector<hHalfedge> GetBorderCycle(hHalfedge border, Polyhedron &mesh)
    {
        std::vector<hHalfedge> borders;
        for(auto hh : CGAL::halfedges_around_face(border, mesh))
        {
            borders.push_back(hh);
        }
        return borders;
    }

    bool WriteHalfEdges(const std::vector<hHalfedge> &edges, const std::string &path)
    {
        std::ofstream ofs(path);
        if (ofs.fail())
        {
            return false;
        }
        for (size_t i = 0; i < edges.size(); i++)
        {
            auto p0 = edges[i]->vertex()->point();
            auto p1 = edges[(i + 1) % edges.size()]->vertex()->point();
            ofs << "v " << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
            ofs << "v " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
            ofs << "l " << i * 2 + 1 << " " << i * 2 + 2 << "\n";
        }
        ofs.close();
        if (ofs.fail())
        {
            return false;
        }
        return true;
    }

    bool WriteHalfEdges(const std::vector<std::vector<hHalfedge>> &edges, const std::string &path)
    {
        std::ofstream ofs(path);
        if (ofs.fail())
        {
            return false;
        }
        int idx = 1;
        for (const std::vector<hHalfedge> &arr : edges)
        {
            for (size_t i = 0; i < arr.size(); i++)
            {
                auto p0 = arr[i]->vertex()->point();
                auto p1 = arr[(i + 1) % arr.size()]->vertex()->point();
                ofs << "v " << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
                ofs << "v " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
                ofs << "l " << idx++ << " " << idx++ << "\n";
            }
        }
        ofs.close();
        if (ofs.fail())
        {
            return false;
        }
        return true;
    }

    bool WritePoints(const std::vector<Point_3> &points, const std::string &path)
    {
        std::ofstream ofs(path);
        if (ofs.fail())
        {
            return false;
        }
        for (size_t i = 0; i < points.size(); i++)
        {
            auto p0 = points[i];
            auto p1 = points[(i + 1) % points.size()];
            ofs << "v " << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
            ofs << "v " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
            ofs << "l " << i * 2 + 1 << " " << i * 2 + 2 << "\n";
        }
        ofs.close();
        if (ofs.fail())
        {
            return false;
        }
        return true;
    }

    bool WritePoints(const std::vector<std::vector<Point_3>> &points, const std::string &path)
    {
        std::ofstream ofs(path);
        if (ofs.fail())
        {
            return false;
        }
        int idx = 1;
        for (const std::vector<Point_3> &arr : points)
        {
            for (size_t i = 0; i < arr.size(); i++)
            {
                auto p0 = arr[i];
                auto p1 = arr[(i + 1) % arr.size()];
                ofs << "v " << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
                ofs << "v " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
                ofs << "l " << idx++ << " " << idx++ << "\n";
            }
        }
        ofs.close();
        if (ofs.fail())
        {
            return false;
        }
        return true;
    }

    void CloseHoles(Polyhedron &mesh)
    {
        std::vector<hHalfedge> border_halfedges;
        CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_halfedges));
        for (auto hh : border_halfedges)
        {
            std::vector<hFacet> out_faces;
            CGAL::Polygon_mesh_processing::triangulate_hole(mesh, hh, std::back_inserter(out_faces));
        }
    }

    class Curve
    {
    public:
        void AddPoint(Point_3 p, int label)
        {
            _points.push_back(p);
            _labels.push_back(label);
        }
        Point_3& operator[](size_t i)
        {
            return _points[i];
        }
        const Point_3& operator[](size_t i) const
        {
            return _points[i];
        }
        int& Label(size_t i)
        {
            return _labels[i];
        }
        const int& Label(size_t i) const
        {
            return _labels[i];
        }
        void UpdateData()
        {
            std::unordered_map<int, int> counts;
            
            for(size_t i = 0; i < _points.size(); i++)
            {
                counts[_labels[i]]++;
                if(_centroids.count(_labels[i]) == 0)
                {
                    _centroids.insert({_labels[i], KernelEpick::Vector_3(0.0, 0.0, 0.0)});
                }
                _centroids[_labels[i]] += _points[i] - CGAL::ORIGIN;
            }
            for(auto& [label, c] : _centroids)
            {
                c /= counts[label];
            }

            
        }

    protected:
        std::vector<Point_3> _points;
        std::vector<int> _labels;
        std::unordered_map<int, KernelEpick::Vector_3> _centroids;
        std::unordered_map<int, Vec3> _upwards;
    };

    std::vector<Point_3> Merge(const std::vector<Point_3> &points0, const std::vector<Point_3> &points1)
    {
        std::pair<size_t, size_t> closest_pair;
        double min_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < points0.size(); i++)
        {
            for (size_t j = 0; j < points1.size(); j++)
            {
                double dist = CGAL::squared_distance(points0[i], points1[j]);
                if (dist < min_dist)
                {
                    closest_pair = {i, j};
                    min_dist = dist;
                }
            }
        }

        KernelEpick::Line_3 line(points0[closest_pair.first], points1[closest_pair.second]);
        KernelEpick::Vector_3 line_dir = points1[closest_pair.second] - points0[closest_pair.first];

        double w0 = 1.0f;
        double w1 = 0.0f;
        std::pair<size_t, size_t> end_pair;
        std::pair<size_t, size_t> start_pair;
        double max_score = -999.0;
        for (size_t i = 0; i < points0.size(); i++)
        {
            size_t idx = (closest_pair.first + i) % points0.size();
            KernelEpick::Vector_3 vec = points0[(idx + 1) % points0.size()] - points0[idx];
            double distance_score = std::sqrt(CGAL::squared_distance(line.projection(points0[idx]), points0[idx]));
            double direction_score = CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * distance_score;
            if(score >= max_score)
            {
                max_score = score;
                end_pair.first = idx;
            }
            else
            {
                break;
            }
        }
        max_score = -999.0;
        for (size_t i = 0; i < points1.size(); i++)
        {
            size_t idx = (closest_pair.first + i) % points1.size();
            KernelEpick::Vector_3 vec = points1[(idx + 1) % points1.size()] - points1[idx];
            double distance_score = std::sqrt(CGAL::squared_distance(line.projection(points1[idx]), points1[idx]));
            double direction_score = 180.0 - CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * distance_score;
            if(score >= max_score)
            {
                max_score = score;
                end_pair.second = idx;
            }
            else
            {
                break;
            }
        }
        max_score = -999.0;
        for (size_t i = 0; i < points0.size(); i++)
        {
            size_t idx = (closest_pair.first + points0.size() - i - 1) % points0.size();
            KernelEpick::Vector_3 vec = points0[idx] - points0[(idx + 1) % points0.size()];
            double distance_score = std::sqrt(CGAL::squared_distance(line.projection(points0[idx]), points0[idx]));
            double direction_score = CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * distance_score;
            if(score >= max_score)
            {
                max_score = score;
                start_pair.first = idx;
            }
            else
            {
                break;
            }
        }
        max_score = -999.0;
        for (size_t i = 0; i < points1.size(); i++)
        {
            size_t idx = (closest_pair.second + points1.size() - i - 1) % points1.size();
            KernelEpick::Vector_3 vec = points1[idx] - points1[(idx + 1) % points1.size()];
            double distance_score = std::sqrt(CGAL::squared_distance(line.projection(points1[idx]), points1[idx]));
            double direction_score = 180.0 - CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * distance_score;
            if(score >= max_score)
            {
                max_score = score;
                start_pair.second = idx;
            }
            else
            {
                break;
            }
        }
        
        std::vector<Point_3> result;
        size_t idx = end_pair.first;
        while(idx != start_pair.first)
        {
            result.push_back(points0[idx]);
            idx++;
            if(idx == points0.size())
            {
                idx = 0;
            }
        }
        for(int i = 0; i < 20; i++)
        {
            result.push_back(points0[start_pair.first] + (points1[end_pair.second] - points0[start_pair.first]) * (double)i / 20.0);
        }
        idx = end_pair.second;
        while(idx != start_pair.second)
        {
            result.push_back(points1[idx]);
            idx++;
            if(idx == points1.size())
            {
                idx = 0;
            }
        }
        for(int i = 0; i < 20; i++)
        {
            result.push_back(points1[start_pair.second] + (points0[end_pair.first] - points1[start_pair.second]) * (double)i / 20.0);
        }
        std::ofstream ofs("./merge" + std::to_string(result.size()) + ".obj");
        for (size_t i = 0; i < result.size(); i++)
        {
            ofs << "v " << result[i].x() << ' ' << result[i].y() << ' ' << result[i].z() << '\n';
        }
        return result;
    }
}

bool FixMeshWithLabel(
    const std::vector<KernelEpick::Point_3>& input_vertices,
    const std::vector<TTriangle<size_t>>& input_faces,
    const std::vector<int>& input_labels,
    Polyhedron& output_mesh,
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry
);

bool GumTrimLine(std::string input_file, std::string label_file, std::string output_file, int smooth)
{
    Polyhedron mesh;
    if (!CGAL::IO::read_polygon_mesh(input_file, mesh, CGAL::parameters::verbose(true)))
    {
        throw IOError("Cannot read from file " + input_file);
    }
    if (!mesh.is_valid(false))
    {
        mesh.is_valid(true);
        throw MeshError("Input mesh not valid: " + input_file);
    }
    if (!mesh.is_pure_triangle())
    {
        throw MeshError("Input mesh has non triangle face: " + input_file);
    }
    CloseHoles(mesh);
    CGAL::set_halfedgeds_items_id(mesh);
    printf("Load mesh: V = %zd, F = %zd.\n", mesh.size_of_vertices(), mesh.size_of_facets());
    mesh.LoadLabels(label_file);
    printf("Loaded labels.\n");
    // for(auto hf : CGAL::faces(mesh))
    // {
    //     int l0 = hf->halfedge()->vertex()->_label;
    //     int l1 = hf->halfedge()->next()->vertex()->_label;
    //     int l2 = hf->halfedge()->prev()->vertex()->_label;
    //     if(l0 == 0 || l1 == 0 || l2 == 0)
    //     {
    //         hf->_label = 0;
    //     }
    //     else
    //     {
    //         hf->_label = std::max(l0, std::max(l1, l2));
    //     }
    // }

    using SubMesh = std::vector<hFacet>;
    std::vector<SubMesh> components;

    for (auto hf : CGAL::faces(mesh))
    {
        if (!hf->_processed && hf->_label != 0)
        {
            components.emplace_back(FacialConnectedComponents(hf, mesh));
        }
    }

    if (components.empty())
    {
        throw AlgError("Cannot find gum part");
    }
    printf("Found %zd gum part by label\n", components.size());

    std::sort(components.begin(), components.end(),
        [](auto &lh, auto &rh) 
        {
            int maxl = 0;
            int minl = 100;
            int maxr = 0;
            int minr = 100;
            for(auto hf : lh)
            {
                if(hf->_label != 0)
                {
                    int l = hf->_label;
                    if(l <= 18 && l >= 11)
                    {
                        l = 29 - l;
                    }
                    else if(l >= 31 && l <= 38)
                    {
                        l = 69 - l;
                    }
                    maxl = std::max(maxl, l);
                    minl = std::min(minl, l);
                }
            }
            for(auto hf : rh)
            {
                if(hf->_label != 0)
                {
                    int l = hf->_label;
                    if(l <= 18 && l >= 11)
                    {
                        l = 29 - l;
                    }
                    else if(l >= 31 && l <= 38)
                    {
                        l = 69 - l;
                    }
                    maxr = std::max(maxr, l);
                    minr = std::min(minr, l);
                }
            }
            if(minl != minr)
            {
                return minl < minr;
            }
            return maxl < maxr;
        });

    using AABBPrimitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using AABBTraits = CGAL::AABB_traits<KernelEpick, AABBPrimitive>;
    using AABBTree = CGAL::AABB_tree<AABBTraits>;
    AABBTree aabb_tree(mesh.facets_begin(), mesh.facets_end(), mesh);
    if (aabb_tree.empty())
    {
        throw AlgError("Failed to build AABB tree.");
    }
    std::vector<std::vector<Point_3>> trim_points;
    for (auto &comp : components)
    {
        if(comp.size() < 100)
        {
            continue;
        }
        CGAL::Face_filtered_graph<Polyhedron> filtered_graph(mesh, comp);
        if (!filtered_graph.is_selection_valid())
        {
            throw AlgError("Invalid part selection!");
        }
        Polyhedron gum_mesh;
        std::vector<std::pair<hVertex, hVertex>> vtx_map;
        CGAL::copy_face_graph(filtered_graph, gum_mesh, CGAL::parameters::vertex_to_vertex_output_iterator(std::back_inserter(vtx_map)));
        for(auto& [source, target] : vtx_map)
        {
            target->_label = source->_label;
        }
        auto [gum_mesh_vertices, gum_mesh_faces] = gum_mesh.ToVerticesTriangles();
        FixMeshWithLabel(gum_mesh_vertices, gum_mesh_faces, gum_mesh.WriteLabels(), gum_mesh, false, 0, false, true, 100, 200, false, 10);
        if (gum_mesh.is_empty() || !gum_mesh.is_valid())
        {
            throw AlgError("Cannot find gum part");
        }
        printf("Mesh valid. F = %zd\n", gum_mesh.size_of_facets());
        std::vector<hHalfedge> border_halfedges;
        CGAL::Polygon_mesh_processing::extract_boundary_cycles(gum_mesh, std::back_inserter(border_halfedges));
        if (border_halfedges.empty())
        {
            throw AlgError("Cannot find trim line");
        }
        printf("Found %zd borders\n", border_halfedges.size());
        std::vector<std::vector<hHalfedge>> border_cycles;
        for (hHalfedge hh : border_halfedges)
        {
            border_cycles.emplace_back(GetBorderCycle(hh, gum_mesh));
        }

        printf("Found %zd possible trim line. ", border_cycles.size());
        border_cycles.erase(std::remove_if(border_cycles.begin(), border_cycles.end(), [](std::vector<hHalfedge> &edges)
                                           { return edges.size() <= 10; }), border_cycles.end());
        printf("Use %zd of them after removing small ones.\n", border_cycles.size());
        if (border_cycles.empty())
        {
            throw AlgError("No valid trim line after filtering.");
        }
        for (std::vector<hHalfedge> &trimline : border_cycles)
        {
            std::vector<Point_3> points;
            for (auto hh : trimline)
            {
                points.push_back(hh->vertex()->point());
            }
            for (size_t iteration = 0; iteration < smooth; iteration++)
            {
                std::vector<Point_3> new_points = points;
                for (size_t i = 0; i < new_points.size(); i++)
                {
                    size_t prev = i == 0 ? new_points.size() - 1 : i - 1;
                    size_t next = i == new_points.size() - 1 ? 0 : i + 1;
                    new_points[i] = CGAL::midpoint(points[prev], points[next]);
                    new_points[i] = aabb_tree.closest_point(new_points[i]);
                }
                points = new_points;
            }
            trim_points.push_back(std::move(points));
        }
    }

    if (trim_points.size() >= 2)
    {
        //Merge(trim_points[0], trim_points[1]);
        for(int i = 1; i < trim_points.size(); i++)
        {
            trim_points[0] = Merge(trim_points[0], trim_points[i]);
        }
    }
    trim_points.resize(1);

    return WritePoints(trim_points, output_file);
}

#ifndef FOUND_PYBIND11
int main(int argc, char *argv[])
{
    std::filesystem::current_path(R"(D:\dev\Ortho\OrthoMeshTools\test\GumTrimLine)");
    std::string input_file = "7.obj";
    std::string label_file = "7.json";
    std::string output_file = "7gumline.obj";
    int smooth = 20;
    for (int i = 1; i < argc; i++)
    {
        if (std::strcmp(argv[i], "-i") == 0)
        {
            input_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-o") == 0)
        {
            output_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-l") == 0)
        {
            label_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-s") == 0)
        {
            smooth = std::atoi(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-h") == 0)
        {
            std::cout << "Extract gum trim line and output it as an .obj file.\n"
                         "-i : file path of input mesh\n"
                         "-o : file path of output mesh\n"
                         "-l : file path of label file.\n"
                         "-s : a non-nagetive integer that specifies the iteration number of trim line smoothing."
                      << std::endl;
            return 0;
        }
    }
    if (input_file.empty() || output_file.empty() || label_file.empty() || smooth < 0)
    {
        std::cout << "Invalid paramters. Use -h for help." << std::endl;
        return -1;
    }
    try
    {
        GumTrimLine(input_file, label_file, output_file, smooth);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}
#endif
