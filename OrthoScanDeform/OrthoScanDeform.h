#ifndef ORTHO_SCAN_DEFORM_H
#define ORTHO_SCAN_DEFORM_H
#include "../Polyhedron.h"
#include <CGAL/bounding_box.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <eigen3/Eigen/Eigen>
#include "../EasyOBJ.h"

template <typename MeshType, int WNum>
class OrthoScanDeform
{
public:
    using Mesh = MeshType;
    using Kernel = typename Mesh::Traits::Kernel;
    struct DeformWeights
    {
        DeformWeights()
        {
            std::fill(label.begin(), label.end(), 0);
            std::fill(w.begin(), w.end(), 0.0);
        }
        std::array<int, WNum> label;
        std::array<double, WNum> w;
    };

    void SetMesh(const MeshType* mesh, const CrownFrames<Kernel>* crown_frames, bool upper)
    {
        this->_mesh = mesh;
        this->_crown_frames = crown_frames;
        this->_deform_weights.clear();
        if(crown_frames->Frames().size() < WNum)
        {
            throw AlgError("Invalid crown frame data.");
        }
        
        std::unordered_map<int, std::vector<typename MeshType::Vertex_handle>> borders;
        for(auto hh : CGAL::halfedges(*mesh))
        {
            if(hh->facet() != nullptr && hh->opposite() != nullptr && hh->opposite()->facet() != nullptr)
            {
                if(hh->facet()->_label == 0 && hh->opposite()->facet()->_label != 0)
                {
                    borders[hh->opposite()->facet()->_label].push_back(hh->vertex());
                }
            }
        }
        
        using LA = CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<typename CGAL::Eigen_sparse_matrix<double>::EigenType>>;
        using HeatMethod = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<MeshType, CGAL::Heat_method_3::Intrinsic_Delaunay, typename boost::property_map<MeshType, CGAL::vertex_point_t>::const_type, LA, Kernel>;

        std::unordered_map<int, std::unordered_map<typename MeshType::Vertex_handle, double>> distances;
        HeatMethod heat_method(*mesh);
        for(auto& [label, vertices] : borders)
        {
            heat_method.clear_sources();
            heat_method.add_sources(vertices);
            heat_method.estimate_geodesic_distances(boost::make_assoc_property_map(distances[label]));
        }

        for(auto hv : CGAL::vertices(*mesh))
        {
            DeformWeights dw;
            if(hv->_label != 0)
            {
                dw.label[0] = hv->_label;
                dw.w[0] = 1.0;
            }
            else
            {
                std::vector<std::pair<int, double>> dists_to_labels;
                for(auto& [label, map] : distances)
                {
                    dists_to_labels.push_back({label, map[hv]});
                }
                std::sort(dists_to_labels.begin(), dists_to_labels.end(), [](auto& lh, auto& rh) { return lh.second < rh.second; });
                for(int i = 0; i < WNum; i++)            
                {
                    dw.label[i] = dists_to_labels[i].first;
                    dw.w[i] = dists_to_labels[i].second;
                }
                double w_sum = 0.0;
                for(int i = 0; i < WNum; i++)
                {
                    if(dw.w[i] <= 1e-8)
                    {
                        std::fill(dw.w.begin(), dw.w.end(), 0.0);
                        dw.w[i] = 1.0;
                        w_sum = 1.0;
                        break;
                    }
                    else
                    {
                        dw.w[i] = std::exp(-dw.w[i] * dw.w[i] / 4.0);
                        w_sum += dw.w[i];
                    }
                }
                for(int i = 0; i < WNum; i++)
                {
                    dw.w[i] /= w_sum;
                }
            }
            _deform_weights.push_back(dw);
        }

        for(auto hv : CGAL::vertices(*mesh))
        {
            auto p = hv->point() - CGAL::ORIGIN;
            auto avg = Kernel::Vector_3(0.0, 0.0, 0.0);
            int count = 0;
            for(auto nei : CGAL::vertices_around_target(hv, *mesh))
            {
                avg += nei->point() - CGAL::ORIGIN;
                count++;
            }
            avg /= count;
            _constraints.push_back(p - avg);
        }
        
        EasyOBJ testout("./weighted_vertices.obj");
        int cnt = 0;
        for(auto hv : CGAL::vertices(*mesh))
        {
            std::array<double, 3> c = {0.0, 0.0, 0.0};
            auto& dw = _deform_weights[cnt++];
            for(int i = 0; i < WNum; i++)
            {
                auto cc = LabelColorMap(dw.label[i]);
                c[0] += cc[0] * dw.w[i];
                c[1] += cc[1] * dw.w[i];
                c[2] += cc[2] * dw.w[i];
                // c[0] += cc[0];
                // c[1] += cc[1];
                // c[2] += cc[2];
            }
            testout.AddV(hv->point(), c[0], c[1], c[2]);
        }
    }

    void Deform(const CrownFrames<Kernel>& frames)
    {
        MeshType mesh = *_mesh;
        int count = 0;
        for(auto hv : CGAL::vertices(mesh))
        {
            auto& weights = _deform_weights[count++];
            auto p = hv->point();
            typename Kernel::Point_3 result = CGAL::ORIGIN;
            for(int i = 0; i < WNum; i++)
            {
                int label = weights.label[i];
                if(label != 0 && label != 1)
                {
                    double w = weights.w[i];
                    result += w * (frames.GetFrame(label).LocalToWorld().transform(CGAL::ORIGIN + (p - _crown_frames->GetFrame(label).pos)) - CGAL::ORIGIN);
                }
                else if(label == 1)
                {
                    double w = weights.w[i];
                    result += w * (p - CGAL::ORIGIN);
                }
            }
            hv->point() = result;
        }
        for(int i = 0; i < 0; i++)
        {
            std::vector<typename Kernel::Point_3> new_points(mesh.points_begin(), mesh.points_end());
            int cnt = 0;
            for(auto hv : CGAL::vertices(mesh))
            {
                if(hv->_label != 0)
                {
                    cnt++;
                    continue;
                }
                typename Kernel::Vector_3 avg(0.0, 0.0, 0.0);
                int nei_count = 0;
                for(auto nei : CGAL::vertices_around_target(hv, mesh))
                {
                    avg += nei->point() - CGAL::ORIGIN;
                    nei_count++;
                }
                avg /= nei_count;
                typename Kernel::Point_3 target = CGAL::ORIGIN + avg + _constraints[cnt];
                new_points[cnt] = CGAL::midpoint(target, new_points[cnt]);
                cnt++;
            }
            cnt = 0;
            for(auto hv : CGAL::vertices(mesh))
            {
                hv->point() = new_points[cnt++];
            }
        }
        mesh.WriteOBJ("./deformed.obj");
    }

protected:
    const MeshType* _mesh = nullptr;
    const CrownFrames<Kernel>* _crown_frames = nullptr;
    std::vector<DeformWeights> _deform_weights;
    std::vector<typename Kernel::Vector_3> _constraints;
};

#endif
