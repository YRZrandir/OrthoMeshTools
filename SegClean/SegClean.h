#ifndef SEGCLEAN_H
#define SEGCLEAN_H
#include <string>
#include <queue>
#include "../Polyhedron.h"

bool SegCleanF(std::string input_mesh, std::string input_labels, std::string output_labels, int size_threshold);

template <typename Polyhedron>
std::vector<typename Polyhedron::Vertex_handle> ConnectedComponents(typename Polyhedron::Vertex_handle hv, Polyhedron &mesh)
{
    using hVertex = typename Polyhedron::Vertex_handle;
    std::vector<hVertex> vertices;
    std::queue<hVertex> front;
    front.push(hv);
    vertices.push_back(hv);
    hv->_processed = true;

    while (!front.empty())
    {
        hVertex v = front.front();
        front.pop();
        for (auto nei : CGAL::vertices_around_target(v, mesh))
        {
            if (!nei->_processed && nei->_label == hv->_label)
            {
                front.push(nei);
                vertices.push_back(nei);
                nei->_processed = true;
            }
        }
    }

    return vertices;
}

template <typename Polyhedron>
void CleanSmallComponents(const std::vector<typename Polyhedron::Vertex_handle> &component, Polyhedron &mesh, 
                            const std::unordered_set<typename Polyhedron::Vertex_handle>& all_vertices_to_clean)
{
    using hVertex = typename Polyhedron::Vertex_handle;
    std::vector<hVertex> surroundings;
    for (auto hv : component)
    {
        for (auto nei : CGAL::vertices_around_target(hv, mesh))
        {
            if (nei->_label != hv->_label && all_vertices_to_clean.count(nei) == 0)
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

template <typename Polyhedron>
void SegClean(Polyhedron& mesh, size_t size_threshold = 1000)
{
    using hVertex = typename Polyhedron::Vertex_handle;
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
        printf("Error: Failed to find connected components.\n");
        return;
    }
    else
    {
        printf("Found %zd connected components by label\n", connected_components.size());
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
        size = std::min(size, size_threshold);
    }

    std::vector<const std::vector<hVertex>*> components_to_clean;
    std::unordered_set<hVertex> vertices_to_clean;
    int cnt = 0;
    for (const auto &[label, vertices] : connected_components)
    {
        if (vertices.size() < max_label_component_sizes[label])
        {
            cnt++;
            components_to_clean.push_back(&vertices);
            vertices_to_clean.insert(vertices.begin(), vertices.end());
        }
    }
    
    for(auto vertices : components_to_clean)
    {
        CleanSmallComponents(*vertices, mesh, vertices_to_clean);
    }
    printf("Cleaned %d small components.\n", cnt);
}
#endif