#include <atomic>
#include <cctype>
#include <fstream>
#include <filesystem>
#include <functional>
#include <iostream>
#include <list>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <omp.h>
#include "../Polyhedron.h"
#include "MeshFix.h"
#ifdef FOUND_PYBIND11
#include <pybind11/pybind11.h>
#endif
// TODO: remove this
bool gVerbose = true;
namespace 
{
using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;
using Triangle = Polyhedron::Triangle;
}

bool FixMeshFile(
    std::string input_mesh,
    std::string output_mesh, 
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry)
{
    try
    {
        std::vector<KernelEpick::Point_3> vertices;
        std::vector<Triangle> faces;
        LoadVFAssimp<KernelEpick, Triangle::size_type>(input_mesh, vertices, faces);
        if(gVerbose)
        {
            printf("Load mesh: V = %zd, F = %zd\n", vertices.size(), faces.size());
        }
        Polyhedron result;
        FixMesh<Polyhedron>(vertices, faces, result, keep_largest_connected_component,
        large_cc_threshold, fix_self_intersection,
        filter_small_holes, max_hole_edges, max_hole_diam, refine, max_retry);
        result.WriteAssimp(output_mesh);
        if(gVerbose)
        {
            printf("Output V = %zd, F = %zd.\n", result.size_of_vertices(), result.size_of_facets());
        }
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
    return true;
}

bool FixMeshFileWithLabel(
    std::string input_mesh,
    std::string output_mesh,
    std::string input_label,
    std::string output_label,
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry
)
{
    std::vector<KernelEpick::Point_3> vertices;
    std::vector<Triangle> faces;
    LoadVFAssimp<KernelEpick, Triangle::size_type>(input_mesh, vertices, faces);
    std::vector<int> labels = LoadLabels(input_label);
    if(gVerbose)
    {
        printf("Load mesh: V = %zd, F = %zd\n", vertices.size(), faces.size());
    }
    Polyhedron m;
    FixMeshWithLabel<Polyhedron>(vertices, faces, labels, m, keep_largest_connected_component, large_cc_threshold,
     fix_self_intersection, filter_small_holes, max_hole_edges, max_hole_diam, refine, max_retry);
    
    m.WriteAssimp(output_mesh);
    if(gVerbose)
    {
        printf("Output V = %zd, F = %zd.\n", m.size_of_vertices(), m.size_of_facets());
    }
    m.WriteLabels(output_label, input_label);
    return true;
}

bool FixMeshFileWithColor(
    std::string input_mesh,
    std::string output_mesh,
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry
)
{
    std::vector<KernelEpick::Point_3> vertices;
    std::vector<Triangle> faces;
    std::vector<Eigen::Vector3f> colors;
    LoadVCFAssimp<KernelEpick, Triangle::size_type>(input_mesh, vertices, colors, faces);

    auto color_to_int = [](Eigen::Vector3f c)->int
    {
        unsigned int r = static_cast<unsigned int>(c.x() * 255);
        unsigned int g = static_cast<unsigned int>(c.y() * 255);
        unsigned int b = static_cast<unsigned int>(c.z() * 255);
        int value = 0x00000000;
        value |= r;
        value = value << 8;
        value |= g;
        value = value << 8;
        value |= b;
        return value;
    };

    auto int_to_color = [](int v)->Eigen::Vector3f
    {
        unsigned int b = v & 0x000000FF;
        v = v >> 8;
        unsigned int g = v & 0x000000FF;
        v = v >> 8;
        unsigned int r = v & 0x000000FF;
        return Eigen::Vector3f((float)r / 255.f, (float)g / 255.f, (float)b / 255.f);
    };

    std::vector<int> labels;
    for(int i = 0; i < vertices.size(); i++)
    {
        labels.push_back(color_to_int(colors[i]));
    }
    if(gVerbose)
    {
        printf("Load mesh: V = %zd, F = %zd\n", vertices.size(), faces.size());
    }
    Polyhedron m;
    FixMeshWithLabel<Polyhedron>(vertices, faces, labels, m, keep_largest_connected_component, large_cc_threshold,
     fix_self_intersection, filter_small_holes, max_hole_edges, max_hole_diam, refine, max_retry);
    
    auto [out_v, out_f] = m.ToVerticesTriangles();
    std::vector<Eigen::Vector3f> out_c;
    for(auto hv : CGAL::vertices(m))
    {
        out_c.push_back(int_to_color(hv->_label));
    }

    WriteVCFAssimp<typename Polyhedron::Traits::Kernel, size_t>(output_mesh, out_v, out_c, out_f);
    if(gVerbose)
    {
        printf("Output V = %zd, F = %zd.\n", m.size_of_vertices(), m.size_of_facets());
    }
    return true;
}