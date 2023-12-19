#include <deque>
#include <iostream>
#include <filesystem>
#include <unordered_map>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include "../MeshFix/MeshFix.h"
#include "GumTrimLine.h"
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

}

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
    //CloseHoles(mesh);
    CGAL::set_halfedgeds_items_id(mesh);
    printf("Load mesh: V = %zd, F = %zd.\n", mesh.size_of_vertices(), mesh.size_of_facets());
    mesh.LoadLabels(label_file);
    mesh.UpdateFaceLabels2();
    printf("Loaded labels.\n");

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
    std::vector<internal::Curve<Polyhedron::Traits>> trim_points;
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
        std::vector<std::pair<hFacet, hFacet>> facet_map;
        CGAL::copy_face_graph(filtered_graph, gum_mesh,
         CGAL::parameters::vertex_to_vertex_output_iterator(std::back_inserter(vtx_map)).face_to_face_output_iterator(std::back_inserter(facet_map)));
        for(auto& [source, target] : vtx_map)
        {
            target->_label = source->_label;
        }
        for(auto& [source, target] : facet_map)
        {
            target->_label = source->_label;
        }
        auto [gum_mesh_vertices, gum_mesh_faces] = gum_mesh.ToVerticesTriangles();
        FixMeshWithLabel(gum_mesh_vertices, gum_mesh_faces, gum_mesh.WriteLabels(), gum_mesh, false, 0, false, true, 100, 200, false, 10);
        if (gum_mesh.is_empty() || !gum_mesh.is_valid())
        {
            throw AlgError("Cannot find gum part");
        }
        gum_mesh.UpdateFaceLabels2();
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
            internal::Curve<Polyhedron::Traits> curve;
            for (auto hh : trimline)
            {
                curve.AddPoint(hh->vertex()->point(), hh->opposite()->facet()->_label);
            }

            for (size_t iteration = 0; iteration < smooth; iteration++)
            {
                std::vector<Point_3> new_points = curve.GetPoints();
                for (size_t i = 0; i < new_points.size(); i++)
                {
                    size_t prev = i == 0 ? new_points.size() - 1 : i - 1;
                    size_t next = i == new_points.size() - 1 ? 0 : i + 1;
                    new_points[i] = CGAL::midpoint(curve[prev], curve[next]);
                    new_points[i] = aabb_tree.closest_point(new_points[i]);
                }
                for(size_t i = 0; i < curve.size(); i++)
                {
                    curve[i] = new_points[i];
                }
            }
            curve.UpdateData();
            //curve.WriteOBJ("curve" + std::to_string(curve.size()) + ".obj");
            trim_points.push_back(curve);
        }
    }

    if (trim_points.size() >= 2)
    {
        for(int i = 1; i < trim_points.size(); i++)
        {
            trim_points[0] = Merge(trim_points[0], trim_points[i]);
        }
    }

    for(size_t i = 0; i < trim_points[0].size(); i++)
    {
        trim_points[0][i] = aabb_tree.closest_point(trim_points[0][i]);
    }

    return WritePoints(trim_points[0].GetPoints(), output_file);
}

#ifndef FOUND_PYBIND11
int main(int argc, char *argv[])
{
    std::filesystem::current_path(R"(D:\dev\Ortho\OrthoMeshTools\test\GumTrimLine)");
    std::string input_file = "0.obj";
    std::string label_file = "0.json";
    std::string output_file = "0gumline.obj";
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
        GumTrimLine("0.obj", "0.json", "0gumline.obj", smooth);
        GumTrimLine("1.obj", "1.json", "1gumline.obj", smooth);
        GumTrimLine("2.obj", "2.json", "2gumline.obj", smooth);
        GumTrimLine("3.obj", "3.json", "3gumline.obj", smooth);
        GumTrimLine("4.ply", "4.json", "4gumline.obj", smooth);
        GumTrimLine("5.ply", "5.json", "5gumline.obj", smooth);
        GumTrimLine("6.obj", "6.json", "6gumline.obj", smooth);
        GumTrimLine("7.obj", "7.json", "7gumline.obj", smooth);

    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}
#endif
