#include <deque>
#include <iostream>
#include <vector>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/version.h>
#include <nlohmann/json.hpp>
#include "../Polyhedron.h"
#ifdef FOUND_PYBIND11
#include <pybind11/pybind11.h>
#endif

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
    using Triangle = Polyhedron::Triangle;
    using Edge = Polyhedron::Edge;
    using PairHashUnordered = Polyhedron::PairHashUnordered;
    using PairPredUnordered = Polyhedron::PairPredUnordered;
    using PairHash = Polyhedron::PairHash;
    using PairPred = Polyhedron::PairPred;
        
    struct Config
    {
        std::string input;
        std::string output;
        std::string labels;
        std::string outmesh;
        int size_threshold;
    };

    Config LoadConfig(int argc, char *argv[])
    {
        Config cfg;
        for (int i = 1; i < argc; i++)
        {
            if (std::strcmp(argv[i], "-i") == 0)
            {
                cfg.input = std::string(argv[i + 1]);
            }
            else if (std::strcmp(argv[i], "-o") == 0)
            {
                cfg.output = std::string(argv[i + 1]);
            }
            else if (std::strcmp(argv[i], "-l") == 0)
            {
                cfg.labels = std::string(argv[i + 1]);
            }
            else if (std::strcmp(argv[i], "-m") == 0)
            {
                cfg.outmesh = std::string(argv[i + 1]);
            }
            else if (std::strcmp(argv[i], "-t") == 0)
            {
                cfg.size_threshold = std::atoi(argv[i + 1]);
            }
        }
        return cfg;
    }

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

    void CleanSmallComponents(const std::vector<hVertex> &component, Polyhedron &mesh)
    {
        std::vector<hVertex> surroundings;

        for (auto hv : component)
        {
            for (auto nei : CGAL::vertices_around_target(hv, mesh))
            {
                if (nei->_label != hv->_label)
                    surroundings.push_back(nei);
            }
        }

        if (surroundings.empty())
        {
            return;
        }

        for (auto hv : component)
        {
            auto p = hv->point();
            auto nearest = *std::min_element(surroundings.begin(), surroundings.end(),
                                             [&](hVertex nei0, hVertex nei1)
                                             { return CGAL::squared_distance(p, nei0->point()) < CGAL::squared_distance(p, nei1->point()); });

            hv->_label = nearest->_label;
        }
    }
}

bool SegClean(std::string input_mesh, std::string input_labels, std::string output_labels, int size_threshold)
{
    // Config cfg = LoadConfig(argc, argv);
    Polyhedron mesh;
    if(CGAL::IO::read_polygon_mesh(input_mesh, mesh))
    {
        printf_s("Load mesh: V = %zd, F = %zd\n", mesh.size_of_vertices(), mesh.size_of_facets());
    }
    else
    {
        printf_s("Error: failed to read mesh: %s\n", input_mesh.c_str());
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
    if(!mesh.LoadLabels(input_labels))
    {
        printf_s("Error: failed to load labels from: %s\n", input_labels.c_str());
        return false;
    }

    std::cout << "Find connected components...";
    std::vector<std::pair<int, std::vector<hVertex>>> connected_components;
    for (auto hv : CGAL::vertices(mesh))
    {
        if (!hv->_processed)
        {
            connected_components.emplace_back(hv->_label, ConnectedComponents(hv, mesh));
        }
    }
    if(connected_components.empty())
    {
        printf_s("Error: Failed to find connected components.\n");
        return false;
    }
    else
    {
        printf_s("Found %zd connected components by label\n", connected_components.size());
    }

    std::unordered_map<int, size_t> max_label_component_sizes;
    for (const auto &[label, vertices] : connected_components)
    {
        if (max_label_component_sizes.count(label) == 0)
            max_label_component_sizes[label] = vertices.size();
        else
            max_label_component_sizes[label] = std::max(vertices.size(), max_label_component_sizes[label]);
    }

    // Set Component Size Thresholds
    for (auto &[label, size] : max_label_component_sizes)
    {
        size = std::min(size, (size_t)size_threshold);
    }

    int cnt = 0;
    for (const auto &[label, vertices] : connected_components)
    {
        if (vertices.size() < max_label_component_sizes[label])
        {
            cnt++;
            CleanSmallComponents(vertices, mesh);
        }
    }
    printf_s("Cleaned %d small components.\n", cnt);
    return mesh.WriteLabels(output_labels, input_labels);
}

#ifndef FOUND_PYBIND11
int main(int argc, char *argv[])
{
    Config cfg = LoadConfig(argc, argv);
    SegClean(cfg.input, cfg.labels, cfg.output, cfg.size_threshold);
    return 0;
}
#endif