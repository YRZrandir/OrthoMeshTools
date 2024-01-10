#ifndef ORTHO_SCAN_DEFORM_H
#define ORTHO_SCAN_DEFORM_H
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include "../Polyhedron.h"
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

    void SetMesh(const MeshType* mesh, const CrownFrames<Kernel>* crown_frames)
    {
        this->_mesh = mesh;
        this->_crown_frames = crown_frames;
        this->_deform_weights.clear();
        if(crown_frames->Frames().size() < WNum)
        {
            throw AlgError("Invalid crown frame data.");
        }
        // test upper or lower ortho mesh
        bool upper = true;
        for(auto hv : CGAL::vertices(*mesh))
        {
            if(hv->_label != 0)
            {
                if(hv->_label <= 30)
                {
                    upper = true;
                    break;
                }
                else
                {
                    upper = false;
                    break;
                }
            }
        }

        EasyOBJ testout("./weighted_vertices.obj");
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
                std::vector<std::pair<int, double>> dist_to_axis;
                for(auto pair = crown_frames->Begin(); pair != crown_frames->End(); pair++)
                {
                    if(upper && pair->first < 30 || !upper && pair->first > 30)
                    {
                        dist_to_axis.push_back({pair->first, CGAL::squared_distance(hv->point(), Kernel::Line_3(pair->second.pos, pair->second.up))});
                    }
                }
                std::sort(dist_to_axis.begin(), dist_to_axis.end(), [](auto& lh, auto& rh) { return lh.second < rh.second; });

                for(int i = 0; i < WNum; i++)                
                {
                    dw.label[i] = dist_to_axis[i].first;
                    dw.w[i] = dist_to_axis[i].second;
                }
                double w_sum = 0.0;
                for(int i = 0; i < WNum; i++)
                {
                    if(dw.w[i] <= 1e-6)
                    {
                        std::fill(dw.w.begin(), dw.w.end(), 0.0);
                        dw.w[i] = 1.0;
                        w_sum = 1.0;
                        break;
                    }
                    else
                    {
                        dw.w[i] = 1.0 / dw.w[i];
                        w_sum += dw.w[i];
                    }
                }
                for(int i = 0; i < WNum; i++)
                {
                    dw.w[i] /= w_sum;
                }
            }
            _deform_weights.push_back(dw);

            std::array<double, 3> c = {0.0, 0.0, 0.0};
            for(int i = 0; i < WNum; i++)
            {
                auto cc = LabelColorMap(dw.label[i]);
                c[0] += cc[0] * dw.w[i];
                c[1] += cc[1] * dw.w[i];
                c[2] += cc[2] * dw.w[i];
            }
            testout.AddV(hv->point(), c[0], c[1], c[2]);
        }
    }


protected:
    const MeshType* _mesh = nullptr;
    const CrownFrames<Kernel>* _crown_frames = nullptr;
    std::vector<DeformWeights> _deform_weights;
};

#endif
