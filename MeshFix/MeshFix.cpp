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