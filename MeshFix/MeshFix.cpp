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
namespace 
{
using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;
using Triangle = Polyhedron::Triangle;
}

bool FixMeshFile(
    std::string input_mesh,
    std::string output_mesh, 
    const MeshFixConfig& cfg)
{
    try
    {
        std::vector<KernelEpick::Point_3> vertices;
        std::vector<Triangle> faces;
        LoadVFAssimp<KernelEpick, Triangle::size_type>(input_mesh, vertices, faces);
        if(cfg.verbosity > 0)
        {
            printf("Load mesh: V = %zd, F = %zd\n", vertices.size(), faces.size());
        }
        Polyhedron result;
        FixMesh<Polyhedron>(vertices, faces, result, cfg);
        result.WriteAssimp(output_mesh);
        if(cfg.verbosity > 0)
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
    const MeshFixConfig& cfg)
{
    std::vector<KernelEpick::Point_3> vertices;
    std::vector<Triangle> faces;
    LoadVFAssimp<KernelEpick, Triangle::size_type>(input_mesh, vertices, faces);
    std::vector<int> labels = LoadLabels(input_label);
    if(cfg.verbosity > 0)
    {
        printf("Load mesh: V = %zd, F = %zd\n", vertices.size(), faces.size());
    }
    Polyhedron m;
    FixMeshWithLabel<Polyhedron>(vertices, faces, labels, m, cfg);
    
    m.WriteAssimp(output_mesh);
    if(cfg.verbosity > 0)
    {
        printf("Output V = %zd, F = %zd.\n", m.size_of_vertices(), m.size_of_facets());
    }
    m.WriteLabels(output_label, input_label);
    return true;
}

bool FixMeshFileWithColor(
    std::string input_mesh,
    std::string output_mesh,
    const MeshFixConfig& cfg
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
    if(cfg.verbosity > 0)
    {
        printf("Load mesh: V = %zd, F = %zd\n", vertices.size(), faces.size());
    }
    Polyhedron m;
    FixMeshWithLabel<Polyhedron>(vertices, faces, labels, m, cfg);
    
    auto [out_v, out_f] = m.ToVerticesTriangles();
    std::vector<Eigen::Vector3f> out_c;
    for(auto hv : CGAL::vertices(m))
    {
        out_c.push_back(int_to_color(hv->_label));
    }

    WriteVCFAssimp<typename Polyhedron::Traits::Kernel, size_t>(output_mesh, out_v, out_c, out_f);
    if(cfg.verbosity > 0)
    {
        printf("Output V = %zd, F = %zd.\n", m.size_of_vertices(), m.size_of_facets());
    }
    return true;
}