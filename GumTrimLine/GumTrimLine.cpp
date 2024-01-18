#include <algorithm>
#include <deque>
#include <iostream>
#include <filesystem>
#include <unordered_map>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/bounding_box.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include "../MeshFix/MeshFix.h"
#include "GumTrimLine.h"
#include "../Ortho.h"

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

    std::vector<hFacet> FacialConnectedComponents(hFacet hf, Polyhedron &mesh, std::function<bool(hFacet)> pred)
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
                if (nei != nullptr && !nei->_processed && pred(nei))
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

    void LabelProcessing(Polyhedron& mesh)
    {
        auto aabb = CGAL::bbox_3(mesh.points_begin(), mesh.points_end());
        double threshold = std::max(aabb.x_span(), std::max(aabb.y_span(), aabb.z_span())) / 100.0;
        
        std::unordered_map<hVertex, int> new_label_set;
        for(auto hv : CGAL::vertices(mesh))
        {
            if(hv->_label != 0)
            {
                continue;
            }
            std::unordered_set<hVertex> neighbors;
            std::unordered_set<int> labels;
            std::queue<hVertex> q;
            q.push(hv);
            neighbors.insert(hv);
            labels.insert(hv->_label);
            while(!q.empty())
            {
                hVertex curr = q.front();
                q.pop();
                for(auto nei : CGAL::vertices_around_target(curr, mesh))
                {
                    if(neighbors.count(nei) == 0 && CGAL::squared_distance(nei->point(), hv->point()) < threshold * threshold)
                    {
                        neighbors.insert(nei);
                        labels.insert(nei->_label);
                        q.push(nei);
                    }
                }
            }
            if(labels.size() >= 3)
            {
                double nearest_dist = std::numeric_limits<double>::max();
                hVertex nearest_hv = nullptr;
                for(auto nei : neighbors)
                {
                    if(nei->_label == hv->_label)
                    {
                        continue;
                    }
                    double d = CGAL::squared_distance(nei->point(), hv->point());
                    if(d < nearest_dist)
                    {
                        nearest_dist = d;
                        nearest_hv = nei;
                    }
                }
                new_label_set.insert({hv, nearest_hv->_label});
            }
        }
    
        for(auto [hv, label] : new_label_set)
        {
            hv->_label = label;
        }
    }
}

bool GumTrimLine(std::string input_file, std::string label_file, std::string frame_file, std::string output_file, int smooth, double fix_factor)
{
    Polyhedron mesh;
    if (!CGAL::IO::read_polygon_mesh(input_file, mesh, CGAL::parameters::verbose(true)))
    {
        try
        {
            printf("possible invalid mesh, try fixing...\n");
            std::vector<typename Polyhedron::Traits::Point_3> vertices;
            std::vector<TTriangle<size_t>> faces;
            if(input_file.ends_with(".obj"))
            {
                LoadVFObj<typename Polyhedron::Traits, size_t>(input_file, vertices, faces);
            }
            else
            {
                // TODO: add custom file loading function for more formats, since assimp cannot keep vertex order sometimes.
                LoadVFAssimp<typename Polyhedron::Traits, size_t>(input_file, vertices, faces);
            }
            std::vector<int> labels = LoadLabels(label_file);
            FixMeshWithLabel(vertices, faces, labels, mesh, false, 0, false, false, 0, 0, false, 10);
        }
        catch(const std::exception&)
        {
            throw IOError("Cannot read mesh file or mesh invalid: " + input_file);
        }
    }
    else
    {
        mesh.LoadLabels(label_file);
        //mesh.UpdateFaceLabels2();
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

    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label == 1)
            hv->_label = 0;
    }
    std::unique_ptr<CrownFrames<typename Polyhedron::Traits>> crown_frames = nullptr;
    if(!frame_file.empty())
    {
        crown_frames = std::make_unique<CrownFrames<typename Polyhedron::Traits>>(frame_file);
    }
    CGAL::set_halfedgeds_items_id(mesh);
    printf("Load ortho scan mesh: V = %zd, F = %zd.\n", mesh.size_of_vertices(), mesh.size_of_facets());
    
    LabelProcessing(mesh);
    mesh.UpdateFaceLabels2();
    // mesh.WriteTriSoup("processed_mesh" + std::to_string(mesh.size_of_facets()) + ".obj");

    using SubMesh = std::vector<hFacet>;
    std::vector<SubMesh> components;
    for (auto hf : CGAL::faces(mesh))
    {
        if (!hf->_processed && hf->_label != 0)
        {
            components.emplace_back(FacialConnectedComponents(hf, mesh, [](hFacet nei){ return nei->_label != 0; }));
        }
        if(!hf->_processed && hf->_label == 0)
        {
            components.emplace_back(FacialConnectedComponents(hf, mesh, [](hFacet nei){ return nei->_label == 0; }));
        }
    }
    if (components.empty())
    {
        throw AlgError("Cannot find gum part");
    }

    // clean components
    std::vector<std::unordered_map<int, int>> comp_label_counts;
    std::unordered_map<int, int> max_label_size;
    std::array<bool, 50> label_exist;
    std::fill(label_exist.begin(), label_exist.end(), false);
    for(auto& comp : components)
    {
        comp_label_counts.emplace_back();
        for(auto hf : comp)
        {
            comp_label_counts.back()[hf->_label]++;
            label_exist[hf->_label] = true;
        }
    }
    for(auto& counts : comp_label_counts)
    {
        for(auto& pair : counts)
        {
            max_label_size[pair.first] = std::max(max_label_size[pair.first], pair.second);
        }
    }
    for(int i = 0; i < components.size(); i++)
    {
        std::unordered_set<int> label_to_remove;
        for(auto& [label, cnt] : comp_label_counts[i])
        {
            if(cnt != max_label_size[label])
            {
                label_to_remove.insert(label);
            }
        }
        std::vector<hFacet> face_to_relabel;
        for(auto hf : components[i])
        {
            if(label_to_remove.count(hf->_label))
            {
                face_to_relabel.push_back(hf);
            }
        }
        if(face_to_relabel.empty())
        {
            continue;
        }
        std::unordered_set<hFacet> surroundings;
        for (auto hf : face_to_relabel)
        {
            for (auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
            {
                if(label_to_remove.count(nei->_label) == 0)
                {
                    surroundings.insert(nei);
                }
            }
        }
        if(!surroundings.empty())
        {
            for (auto hf : face_to_relabel)
            {
                auto p = CGAL::centroid(hf->halfedge()->vertex()->point(), hf->halfedge()->next()->vertex()->point(), hf->halfedge()->prev()->vertex()->point());
                auto nearest = *std::min_element(surroundings.begin(), surroundings.end(), [&](hFacet nei0, hFacet nei1) {
                    auto nei0_p = CGAL::centroid(nei0->halfedge()->vertex()->point(), nei0->halfedge()->next()->vertex()->point(), nei0->halfedge()->prev()->vertex()->point());
                    auto nei1_p = CGAL::centroid(nei1->halfedge()->vertex()->point(), nei1->halfedge()->next()->vertex()->point(), nei1->halfedge()->prev()->vertex()->point());
                        return CGAL::squared_distance(p, nei0_p) < CGAL::squared_distance(p, nei1_p); });
                hf->_label = nearest->_label;
                hf->halfedge()->vertex()->_label = nearest->_label;
                hf->halfedge()->next()->vertex()->_label = nearest->_label;
                hf->halfedge()->prev()->vertex()->_label = nearest->_label;
            }
        }
    }

    // map face labels to vertex labels
    for(auto hv : CGAL::vertices(mesh))
    {
        int label = 0;
        for(auto hf : CGAL::faces_around_target(hv->halfedge(), mesh))
        {
            if(hf != nullptr)
            {
                label = std::max(label, hf->_label);
            }
        }
        hv->_label = label;
    }

    // Recompute label components
    for(auto hf : CGAL::faces(mesh))
    {
        hf->_processed = false;
    }
    components.clear();
    for (auto hf : CGAL::faces(mesh))
    {
        if (!hf->_processed && hf->_label != 0)
        {
            components.emplace_back(FacialConnectedComponents(hf, mesh, [](hFacet nei){ return nei->_label != 0; }));
        }
    }
    if (components.empty())
    {
        throw AlgError("Cannot find gum part");
    }
    printf("Divided ortho mesh into %zd submesh by label.\n", components.size());

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
            printf("Skip component with %zd faces.\n", comp.size());
            continue;
        }
        printf("Processing component with %zd faces...\n", comp.size());
        CGAL::Face_filtered_graph<Polyhedron> filtered_graph(mesh, comp);
        if (!filtered_graph.is_selection_valid())
        {
            throw AlgError("Invalid part selection!");
        }
        /* Extract sub mesh */
        Polyhedron part_mesh;
        std::vector<std::pair<hVertex, hVertex>> vtx_map;
        std::vector<std::pair<hFacet, hFacet>> facet_map;
        CGAL::copy_face_graph(filtered_graph, part_mesh,
        CGAL::parameters::vertex_to_vertex_output_iterator(std::back_inserter(vtx_map)).face_to_face_output_iterator(std::back_inserter(facet_map)));
        for(auto& [source, target] : vtx_map)
        {
            target->_label = source->_label;
        }
        for(auto& [source, target] : facet_map)
        {
            target->_label = source->_label;
        }
        
        /* Fix non-manifold */
        auto [gum_mesh_vertices, gum_mesh_faces] = part_mesh.ToVerticesTriangles();
        FixMeshWithLabel(gum_mesh_vertices, gum_mesh_faces, part_mesh.WriteLabels(), part_mesh, true, 0, false, true, 0, 0, false, 10);
        if (part_mesh.is_empty() || !part_mesh.is_valid())
        {
            throw AlgError("Cannot find trim line");
        }
        part_mesh.UpdateFaceLabels2();
        printf("SubMesh valid. F = %zd\n", part_mesh.size_of_facets());
        // part_mesh.WriteOBJ("part_mesh" + std::to_string(part_mesh.size_of_vertices()) + ".obj");
        /* Extract borders */
        std::vector<hHalfedge> border_halfedges;
        CGAL::Polygon_mesh_processing::extract_boundary_cycles(part_mesh, std::back_inserter(border_halfedges));
        std::vector<std::vector<hHalfedge>> border_cycles;
        std::transform(border_halfedges.begin(), border_halfedges.end(), std::back_inserter(border_cycles), [&part_mesh](auto& hf) { return GetBorderCycle(hf, part_mesh);});
        border_cycles.erase(std::remove_if(border_cycles.begin(), border_cycles.end(), [](std::vector<hHalfedge> &edges)
                                           { return edges.size() <= 10; }), border_cycles.end());
        if (border_cycles.empty())
        {
            throw AlgError("No valid trim line.");
        }
        for(auto& trimline : border_cycles)
        {
            if(trimline.size() < 10)
            {
                continue;
            }
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
                    new_points[i] = new_points[i] + 0.5 * (CGAL::midpoint(curve[prev], curve[next]) - curve[i]);
                    new_points[i] = aabb_tree.closest_point(new_points[i]);
                }
                for(size_t i = 0; i < curve.size(); i++)
                {
                    curve[i] = new_points[i];
                }
            }
            curve.UpdateData();
            if(crown_frames != nullptr)
            {
                curve.LoadCrownFrame(*crown_frames);
            }
            //curve.WriteOBJ("curve" + std::to_string(curve.size()) + ".obj");
            trim_points.push_back(curve);
            printf("Added curve of %zd points.\n", curve.size());
        }
    }

    std::sort(trim_points.begin(), trim_points.end(), [](auto& lh, auto& rh){
            int maxl = lh.MaxLabel();
            int minl = lh.MinLabel();
            int maxr = rh.MaxLabel();
            int minr = rh.MinLabel();

            if(maxl <= 18 && maxl >= 11)
                maxl = 29 - maxl;
            else if(maxl >= 31 && maxl <= 38)
                maxl = 69 - maxl;
            if(minl <= 18 && minl >= 11)
                minl = 29 - minl;
            else if(minl >= 31 && minl <= 38)
                minl = 69 - minl;
            
            if(maxr <= 18 && maxr >= 11)
                maxr = 29 - maxr;
            else if(maxr >= 31 && maxr <= 38)
                maxr = 69 - maxr;
            if(minr <= 18 && minr >= 11)
                minr = 29 - minr;
            else if(minr >= 31 && minr <= 38)
                minr = 69 - minr;

            if(minl != minr)
            {
                return minl < minr;
            }
            return maxl < maxr;
    });

    auto& final_curve = trim_points[0];
    if (trim_points.size() >= 2)
    {
        for(int i = 1; i < trim_points.size(); i++)
        {
            trim_points[0] = Merge(trim_points[0], trim_points[i], aabb_tree);
            if(crown_frames != nullptr)
            {
                trim_points[0].LoadCrownFrame(*crown_frames);
            }
        }
        for(size_t i = 0; i < final_curve.size(); i++)
        {
            final_curve[i] = aabb_tree.closest_point(final_curve[i]);
        }
    }
    final_curve.UpdateData();
    if(crown_frames != nullptr)
    {
        final_curve.LoadCrownFrame(*crown_frames);
    }
    if(fix_factor != 0.0)
    {
        final_curve.FixAllCurve(aabb_tree, fix_factor);
    }
    for (size_t iteration = 0; iteration < 0; iteration++)
    {
        std::vector<Point_3> new_points = final_curve.GetPoints();
        for (size_t i = 0; i < new_points.size(); i++)
        {
            size_t prev = i == 0 ? new_points.size() - 1 : i - 1;
            size_t next = i == new_points.size() - 1 ? 0 : i + 1;
            new_points[i] = CGAL::midpoint(final_curve[prev], final_curve[next]);
            new_points[i] = aabb_tree.closest_point(new_points[i]);
        }
        for(size_t i = 0; i < final_curve.size(); i++)
        {
            final_curve[i] = new_points[i];
        }
    }
    for (size_t iteration = 0; iteration < 1; iteration++)
    {
        std::vector<Point_3> new_points = final_curve.GetPoints();
        for (size_t i = 0; i < new_points.size(); i++)
        {
            size_t prev = i == 0 ? new_points.size() - 1 : i - 1;
            size_t next = i == new_points.size() - 1 ? 0 : i + 1;
            new_points[i] = new_points[i] + 0.5 * (CGAL::midpoint(final_curve[prev], final_curve[next]) - final_curve[i]);
        }
        for(size_t i = 0; i < final_curve.size(); i++)
        {
            final_curve[i] = new_points[i];
        }
    }

    return WritePoints(final_curve.GetPoints(), output_file);
}
