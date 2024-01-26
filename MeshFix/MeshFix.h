#ifndef MESH_FIX_H
#define MESH_FIX_H
#include <string>
#include <vector>
#include <unordered_map>
#include "../Polyhedron.h"
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

// TODO: remove this
extern bool gVerbose;
namespace internal
{
template <typename Kernel, typename SizeType>
std::vector<TTriangle<SizeType>> FixRoundingOrder(const std::vector<typename Kernel::Point_3>& vertices, const std::vector<TTriangle<SizeType>>& faces )
{
    std::unordered_map<std::pair<SizeType, SizeType>, std::vector<SizeType>, TPairHash<SizeType>, TPairPred<SizeType>> edgemap;
    std::vector<SizeType> faces_to_remove;
    for(SizeType i = 0; i < faces.size(); i++)
    {
        auto& f = faces[i];
        edgemap[{f[0], f[1]}].push_back(i);
        edgemap[{f[1], f[2]}].push_back(i);
        edgemap[{f[2], f[0]}].push_back(i);
    }

    std::unordered_set<std::pair<SizeType, SizeType>, TPairHash<SizeType>, TPairPred<SizeType>> problematic_edges;
    for(auto it = edgemap.begin(); it != edgemap.end(); it++)
    {
        if( it->second.size() > 1)
        {
            problematic_edges.insert(it->first);
        }
    }

    std::vector<TTriangle<SizeType>> new_faces;
    for(SizeType i = 0; i < faces.size(); i++)
    {
        auto& f = faces[i];
        if(problematic_edges.count({f[0], f[1]}) > 0 ||
        problematic_edges.count({f[1], f[2]}) > 0 ||
        problematic_edges.count({f[2], f[0]}) > 0)
        {
            continue;
        }
        else
        {
            new_faces.push_back(f);
        }
    }

    return new_faces;
}

template <typename Poly>
#if BOOST_CXX_VERSION >= 202002L
    requires std::derived_from<typename Poly::Items, ItemsWithLabelFlag>
#endif
void FixSelfIntersection( Poly& m, int max_retry )
{
    static_assert(std::is_base_of_v<ItemsWithLabelFlag, typename Poly::Items>);
    std::vector<std::pair<typename Poly::Facet_handle, typename Poly::Facet_handle>> intersect_faces;
    CGAL::Polygon_mesh_processing::self_intersections<CGAL::Parallel_if_available_tag>(m, std::back_inserter(intersect_faces));
    std::unordered_set<typename Poly::Facet_handle> face_to_remove;
    for(auto [f1, f2] : intersect_faces)
    {
        face_to_remove.insert(f1);  
        face_to_remove.insert(f2);
    }
    for(auto& hf : face_to_remove)
    {
        m.erase_facet(hf->halfedge());
    }

    auto [vertices, triangles] = m.ToVerticesTriangles();
    std::vector<int> labels;
    for(auto hv : CGAL::vertices(m))
    {
        labels.push_back(hv->_label);
    }
    size_t nb_removed_faces = 0;
    int cnt = 0;
    do
    {
        triangles = RemoveNonManifold<typename Poly::Vertex::Point_3::R, typename Poly::Face::size_type>(vertices, triangles, &nb_removed_faces);
        if(cnt++ > max_retry)
            break;
    } while(nb_removed_faces != 0);
    
    try
    {
        m.BuildFromVerticesFaces(vertices, triangles);
    }
    catch(const MeshError& e)
    {
        throw AlgError("Failed to fix self intersection: " + std::string(e.what()));
    }
    
    auto hv = m.vertices_begin();
    for(size_t i = 0; i < vertices.size(); i++)
    {
        hv->_label = labels[i];
        hv++;
    }
}

template <typename Kernel, typename SizeType>
std::vector<TTriangle<SizeType>> RemoveNonManifold(const std::vector<typename Kernel::Point_3>& vertices, const std::vector<TTriangle<SizeType>>& faces, size_t* nb_removed_face)
{
    using size_type = SizeType;
    std::vector<std::pair<TTriangle<SizeType>, bool>> faceflags;
    for(auto& f : faces)
    {
        faceflags.push_back(std::make_pair(f, true));
    }

    std::unordered_map<std::pair<SizeType, SizeType>, TEdge<SizeType>, TPairHashUnordered<SizeType>, TPairPredUnordered<SizeType>> edges;
    for(size_t i = 0; i < faceflags.size(); i++)
    {
        const auto& f = faceflags[i].first;
        auto ie0 = edges.find(std::make_pair(f[0], f[1]));
        if(ie0 == edges.end())
        {
            edges[{f[0], f[1]}] = TEdge<SizeType>(f[0], f[1]);
            edges[{f[0], f[1]}]._faces.push_back(i);
        }
        else
        {
            edges[{f[0], f[1]}]._faces.push_back(i);
        }

        auto ie1 = edges.find({f[1], f[2]});
        if(ie1 == edges.end())
        {
            edges[{f[1], f[2]}] = TEdge<SizeType>(f[1], f[2]);
            edges[{f[1], f[2]}]._faces.push_back(i);
        }
        else
        {
            edges[{f[1], f[2]}]._faces.push_back(i);
        }

        auto ie2 = edges.find({f[2], f[0]});
        if(ie2 == edges.end())
        {
            edges[{f[2], f[0]}] = TEdge<SizeType>(f[2], f[0]);
            edges[{f[2], f[0]}]._faces.push_back(i);
        }
        else
        {
            edges[{f[2], f[0]}]._faces.push_back(i);
        }
    }
    
    std::vector<SizeType> problematic_vertices;
    size_t nb_nm_edges = 0;
    for(auto it = edges.begin(); it != edges.end(); it++)
    {
        if(it->second._faces.size() <= 2)
        {
            continue;
        }

        problematic_vertices.push_back(it->first.first);
        problematic_vertices.push_back(it->first.second);

        for(const auto& hf : it->second._faces)
        {
            nb_nm_edges++;
            faceflags[hf].second = false;
        }
    }

    std::vector<std::vector<SizeType>> vneighbors;
    vneighbors.resize(vertices.size());
    for(size_t i = 0; i < faceflags.size(); i++)
    {
        if(faceflags[i].second)
        {
            vneighbors[faceflags[i].first[0]].push_back(i);
            vneighbors[faceflags[i].first[1]].push_back(i);
            vneighbors[faceflags[i].first[2]].push_back(i);
        }
    }

    for(auto pv : problematic_vertices)
    {
        for(auto f : vneighbors[pv])
        {
            faceflags[f].second = false;
        }
    }

    std::atomic_int nb_nm_vertices = 0;
#pragma omp parallel for
    for(int iv = 0; iv < vneighbors.size(); iv++)
    {
        auto& neighbors = vneighbors[iv];
        std::list<std::pair<SizeType, SizeType>> sur_edges;
        size_t nb_connect_faces = neighbors.size();
        std::vector<int> sampled(nb_connect_faces, 0);
        size_t nb_cluster = 0;
        for(size_t i = 0; i < nb_connect_faces; i++)
        {
            if(sampled[i] == 1)
                continue;
            std::list<size_t> cluster;
            cluster.push_back(i);
            sampled[i] = 1;
            do
            {
                auto e0 = faceflags[neighbors[cluster.front()]].first.GetEdge(0);
                auto e1 = faceflags[neighbors[cluster.front()]].first.GetEdge(1);
                auto e2 = faceflags[neighbors[cluster.front()]].first.GetEdge(2);

                for(size_t j = 0; j < nb_connect_faces; j++)
                {
                    if(j != cluster.front() && sampled[j] != 1)
                    {
                        auto e3 = faceflags[neighbors[j]].first.GetEdge(0);
                        auto e4 = faceflags[neighbors[j]].first.GetEdge(1);
                        auto e5 = faceflags[neighbors[j]].first.GetEdge(2);

                        if(TPairPredUnordered<SizeType>()(e0, e3) || TPairPredUnordered<SizeType>()(e0, e4) || TPairPredUnordered<SizeType>()(e0, e5) ||
                        TPairPredUnordered<SizeType>()(e1, e3) || TPairPredUnordered<SizeType>()(e1, e4) || TPairPredUnordered<SizeType>()(e1, e5) ||
                        TPairPredUnordered<SizeType>()(e2, e3) || TPairPredUnordered<SizeType>()(e2, e4) || TPairPredUnordered<SizeType>()(e2, e5))
                        {
                            cluster.push_back(j);
                            sampled[j] = 1;
                        }
                    }
                }
                cluster.pop_front();
            } while(!cluster.empty());
            nb_cluster++;
        }

        if(nb_cluster > 1)
        {
            nb_nm_vertices++;
            for(size_t hf : neighbors)
            {
                faceflags[hf].second = false;
            }
        }
    }

    std::vector<TTriangle<SizeType>> result_faces;
    for(const auto& [face, flag] : faceflags)
    {
        if(flag)
        {
            result_faces.push_back(face);
        }
    }

    if(gVerbose)
    {
        std::cout << "Find " << nb_nm_edges << " non-manifold edges and " << nb_nm_vertices << " non-manifold vertices." << std::endl;
        std::cout << "After remove non-manifold: " << result_faces.size() << " faces." << std::endl;
    }
    *nb_removed_face = faces.size() - result_faces.size();
    return result_faces;
}
}

bool FixMeshFile(
    std::string path,
    std::string output_path,
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry);

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
    int max_retry);

template <typename Polyhedron>
void FixMesh(
    const std::vector<typename Polyhedron::K::Point_3>& input_vertices,
    const std::vector<TTriangle<typename Polyhedron::Vertex::size_type>>& input_faces,
    Polyhedron& output_mesh,
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry,
    bool fair = false,
    std::vector<std::pair<std::vector<typename Polyhedron::Vertex_handle>, std::vector<typename Polyhedron::Facet_handle>>>* patch = nullptr
)
{
    using Kernel = typename Polyhedron::K;
    using Triangle = TTriangle<typename Polyhedron::Vertex::size_type>;
    auto faces = internal::FixRoundingOrder<Kernel, typename Triangle::size_type>(input_vertices, input_faces);
    std::cout << "After fix rounding F = " << faces.size() << std::endl;

    size_t nb_removed_faces = 0;
    int cnt = 0;
    do
    {
        faces = internal::RemoveNonManifold<Kernel, typename Triangle::size_type>(input_vertices, faces, &nb_removed_faces);
        if(cnt++ >= max_retry)
            break;
    } while (nb_removed_faces != 0);
    if(cnt >= max_retry)
    {
        throw AlgError("Cannot remove non-manifold parts. Try increasing retry times.");
    }

    Polyhedron& m = output_mesh;
    m.clear();
    m.BuildFromVerticesFaces(input_vertices, faces);
    
    CGAL::Polygon_mesh_processing::remove_isolated_vertices(m);

    if(fix_self_intersection)
    {
        internal::FixSelfIntersection(m, max_retry);
    }

    if(keep_largest_connected_component)
    {
        size_t num = CGAL::Polygon_mesh_processing::keep_large_connected_components(m, large_cc_threshold);
        if(gVerbose)
        {
            std::cout << "Remove " << num << " small connected components." << std::endl;
        }
    }

    std::vector<typename Polyhedron::Halfedge_handle> border_edges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m, std::back_inserter(border_edges));
    for(typename Polyhedron::Halfedge_handle hh : border_edges)
    {
        if(!filter_small_holes || (filter_small_holes && m.IsSmallHole(hh, max_hole_edges, max_hole_diam)))
        {
            if(refine && fair)
            {
                std::vector<typename Polyhedron::Vertex_handle> patch_vertices;
                std::vector<typename Polyhedron::Facet_handle> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(m, hh,
                 CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)).vertex_output_iterator(std::back_inserter(patch_vertices)).do_not_use_cubic_algorithm(true));
                bool success = !patch_faces.empty();
                if(!success)
                {
                    // relax the border edges and try again
                    std::cout << "failed to close hole, smoothing & trying again." << std::endl;
                    std::vector<typename Polyhedron::Halfedge_handle> borders(CGAL::halfedges_around_face(hh, m).begin(), CGAL::halfedges_around_face(hh, m).end());
                    std::vector<typename Polyhedron::Point_3> points;
                    for(auto edge : borders)
                    {
                        points.push_back(CGAL::midpoint(edge->next()->vertex()->point(), edge->prev()->vertex()->point()));
                    }
                    for(int i = 0; i < points.size(); i++)
                    {
                        borders[i]->vertex()->point() = points[i];
                    }
                    CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(m, hh,
                        CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)).vertex_output_iterator(std::back_inserter(patch_vertices)).do_not_use_cubic_algorithm(true));
                    success = !patch_faces.empty();
                }
                if(!success)
                {
                    std::cout << "Warning: failed to close hole." << std::endl;
                }
                if(patch != nullptr)
                {
                    patch->emplace_back(std::move(patch_vertices), std::move(patch_faces));
                }
            }
            else if(refine)
            {
                std::vector<typename Polyhedron::Vertex_handle> patch_vertices;
                std::vector<typename Polyhedron::Facet_handle> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(m, hh, 
                CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)).vertex_output_iterator(std::back_inserter(patch_vertices)).do_not_use_cubic_algorithm(true));
                if(patch != nullptr)
                {
                    patch->emplace_back(std::move(patch_vertices), std::move(patch_faces));
                }
            }
            else
            {
                std::vector<typename Polyhedron::Facet_handle> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_hole(m, hh,
                    CGAL::parameters::do_not_use_cubic_algorithm(true).face_output_iterator(std::back_inserter(patch_faces)));
            }
        }
    }
}

template <typename Polyhedron>
void FixMeshWithLabel(
    const std::vector<typename Polyhedron::Traits::Point_3> &input_vertices,
    const std::vector<TTriangle<typename Polyhedron::Vertex::size_type>> &input_faces,
    const std::vector<int> &input_labels,
    Polyhedron &output_mesh,
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry,
    std::vector<std::pair<std::vector<typename Polyhedron::Vertex_handle>, std::vector<typename Polyhedron::Facet_handle>>>* patch = nullptr)
{
    using Kernel = typename Polyhedron::Traits;
    using Triangle = TTriangle<typename Polyhedron::Vertex::size_type>;
    auto faces = internal::FixRoundingOrder<Kernel, typename Triangle::size_type>(input_vertices, input_faces);
    std::cout << "After fix rounding F=" << faces.size() << std::endl;
    size_t nb_removed_faces = 0;
    int cnt = 0;
    do
    {
        faces = internal::RemoveNonManifold<Kernel, typename Triangle::size_type>(input_vertices, faces, &nb_removed_faces);
        if (cnt++ >= max_retry)
            break;
    } while (nb_removed_faces != 0);
    if (cnt >= max_retry)
    {
        throw AlgError("Cannot remove non-manifold parts. Try increasing retry times.");
    }
    Polyhedron& m = output_mesh;
    m.clear();
    m.BuildFromVerticesFaces(input_vertices, faces);

    m.LoadLabels(input_labels);

    CGAL::Polygon_mesh_processing::remove_isolated_vertices(m);

    if (fix_self_intersection)
    {
        internal::FixSelfIntersection(m, max_retry);
    }

    if (keep_largest_connected_component)
    {
        size_t num = CGAL::Polygon_mesh_processing::keep_large_connected_components(m, large_cc_threshold);
        if (gVerbose)
        {
            std::cout << "Remove " << num << " small connected components." << std::endl;
        }
    }

    std::vector<typename Polyhedron::Halfedge_handle> border_edges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m, std::back_inserter(border_edges));
    for (typename Polyhedron::Halfedge_handle hh : border_edges)
    {
        if (!filter_small_holes || (filter_small_holes && m.IsSmallHole(hh, max_hole_edges, max_hole_diam)))
        {
            if(refine)
            {
                std::vector<typename Polyhedron::Vertex_handle> patch_vertices;
                std::vector<typename Polyhedron::Facet_handle> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(m, hh, 
                    CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)).vertex_output_iterator(std::back_inserter(patch_vertices)));

                // labels of patch vertices
                std::unordered_set<typename Polyhedron::Vertex_handle> patch_vset(patch_vertices.begin(), patch_vertices.end());
                std::unordered_set<typename Polyhedron::Vertex_handle> neighbor_vset;
                for(auto hv : patch_vset)
                {
                    for(auto nei : CGAL::vertices_around_target(hv, m))
                    {
                        if(patch_vset.count(nei) == 0)
                        {
                            neighbor_vset.insert(nei);
                        }
                    }
                }
                for(auto hv : patch_vset)
                {
                    auto nearest = *std::min_element(neighbor_vset.begin(), neighbor_vset.end(), [&](auto& lh, auto& rh) 
                        { return CGAL::squared_distance(lh->point(), hv->point()) < CGAL::squared_distance(rh->point(), hv->point()); });
                    hv->_label = nearest->_label;
                }
                if(patch != nullptr)
                {
                    patch->emplace_back(std::move(patch_vertices), std::move(patch_faces));
                }
            }
            else
            {
                std::vector<typename Polyhedron::Facet_handle> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_hole(m, hh, CGAL::parameters::face_output_iterator(std::back_inserter(patch_faces)));
                if(patch != nullptr)
                {
                    patch->emplace_back(std::vector<typename Polyhedron::Vertex_handle>(), std::move(patch_faces));
                }
            }
        }
    }
}

#endif
