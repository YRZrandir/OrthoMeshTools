#include <deque>
#include <iostream>
#include <filesystem>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include "../Polyhedron.h"

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
                if (!nei->_processed && nei->_label == hv->_label)
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
                if (nei != nullptr && !nei->_processed && nei->_label == hf->_label)
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
        hHalfedge hh = border;
        std::vector<hHalfedge> borders;
        do
        {
            borders.push_back(hh);
            for (auto nei : CGAL::halfedges_around_source(hh->vertex(), mesh))
            {
                if (nei->is_border())
                {
                    hh = nei;
                    break;
                }
            }
        } while (hh != border);

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

    bool WriteHalfEdges(const std::vector<std::vector<hHalfedge>>& edges, const std::string& path)
    {
        std::ofstream ofs(path);
        if (ofs.fail())
        {
            return false;
        }
        int idx = 1;
        for(const std::vector<hHalfedge>& arr : edges)
        {
            for (size_t i = 0; i < arr.size(); i++)
            {
                auto p0= arr[i]->vertex()->point();
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

    bool WritePoints(const std::vector<std::vector<Point_3>>& points, const std::string& path)
    {
        std::ofstream ofs(path);
        if (ofs.fail())
        {
            return false;
        }
        int idx = 1;
        for(const std::vector<Point_3>& arr : points)
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

    void CloseHoles(Polyhedron& mesh)
    {
        std::vector<hHalfedge> border_halfedges;
        CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_halfedges));
        for(auto hh : border_halfedges)
        {
            std::vector<hFacet> out_faces;
            CGAL::Polygon_mesh_processing::triangulate_hole(mesh, hh, std::back_inserter(out_faces));
        }       
    }

    std::vector<Point_3> Merge( const std::vector<Point_3>& point0, const std::vector<Point_3>& point1 )
    {

    }
}

bool GumTrimLine(std::string input_file, std::string label_file, std::string output_file, int smooth)
{
    Polyhedron mesh;
    if (!CGAL::IO::read_polygon_mesh(input_file, mesh, CGAL::parameters::verbose(true)))
    {
        printf_s("Failed to read mesh: %s\n", input_file.c_str());
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
    CloseHoles(mesh);
    CGAL::set_halfedgeds_items_id(mesh);
    printf_s("Load mesh: V = %zd, F = %zd.\n", mesh.size_of_vertices(), mesh.size_of_facets());
    if (!mesh.LoadLabels(label_file))
    {
        printf_s("Failed to read labels: %s\n", label_file.c_str());
        return false;
    }
    else
    {
        printf_s("Loaded labels.\n");
    }
    for(auto hf : CGAL::faces(mesh))
    {
        int l0 = hf->halfedge()->vertex()->_label;
        int l1 = hf->halfedge()->next()->vertex()->_label;
        int l2 = hf->halfedge()->prev()->vertex()->_label;
        if(l0 == 0 || l1 == 0 || l2 == 0)
        {
            hf->_label = 0;
        }
        else
        {
            hf->_label = std::max(l0, std::max(l1, l2));
        }
    }

    using SubMesh = std::vector<hFacet>;
    std::vector<SubMesh> components;

    for (auto hf : CGAL::faces(mesh))
    {
        if (!hf->_processed && hf->_label == 0)
        {
            components.emplace_back(FacialConnectedComponents(hf, mesh));
        }
    }

    if (components.empty())
    {
        printf_s("Cannot find gum part.\n");
        return false;
    }
    printf_s("Found %zd gum part by label\n", components.size());

    std::sort(components.begin(), components.end(), [](auto &lh, auto &rh)
              { return lh.size() > rh.size(); });
    CGAL::Face_filtered_graph<Polyhedron::Base> filtered_graph(mesh, components[0]);
    Polyhedron gum_mesh;
    CGAL::copy_face_graph(filtered_graph, gum_mesh);
    if(gum_mesh.is_empty() || !gum_mesh.is_valid())
    {
        printf_s("Error: Cannot find gum part\n");
        return false;
    }
    std::vector<hHalfedge> border_halfedges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(gum_mesh, std::back_inserter(border_halfedges));
    if (border_halfedges.empty())
    {
        printf_s("Error: Cannot find trim line\n");
        return false;
    }
    std::vector<std::vector<hHalfedge>> border_cycles;
    for (hHalfedge hh : border_halfedges)
    {
        border_cycles.emplace_back(GetBorderCycle(hh, gum_mesh));
    }

    printf_s("Found %zd possible trim line. ", border_cycles.size());
    border_cycles.erase(std::remove_if(border_cycles.begin(), border_cycles.end(), [](std::vector<hHalfedge>& edges){ return edges.size() <= 10; }), border_cycles.end());
    printf_s("Use %zd of them after removing small ones.\n", border_cycles.size());
    if(border_cycles.empty())
    {
        printf_s("Error: No valid trim line after filtering.\n");
        return false;
    }

    using AABBPrimitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using AABBTraits = CGAL::AABB_traits<KernelEpick, AABBPrimitive>;
    using AABBTree = CGAL::AABB_tree<AABBTraits>;
    AABBTree aabb_tree(mesh.facets_begin(), mesh.facets_end(), mesh);
    if(aabb_tree.empty())
    {
        printf_s("Error: Failed to build AABB tree.\n");
        return false;
    }
    std::vector<std::vector<Point_3>> trim_points;
    for(std::vector<hHalfedge>& trimline : border_cycles)
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
            "-s : a non-nagetive integer that specifies the iteration number of trim line smoothing." << std::endl;
            return 0;
        }
    }
    if(input_file.empty() || output_file.empty() || label_file.empty() || smooth < 0)
    {
        std::cout << "Invalid paramters. Use -h for help." << std::endl;
        return -1;
    }
    if (!GumTrimLine(input_file, label_file, output_file, smooth))
    {
        printf_s("GumTrimLine Failed.\n");
        return -1;
    }
    return 0;
}
#endif
