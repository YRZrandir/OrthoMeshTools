#ifndef ORTHO_SCAN_DEFORM_H
#define ORTHO_SCAN_DEFORM_H
#include "../Polyhedron.h"
#include <CGAL/bounding_box.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>
#include <CGAL/Surface_mesh_deformation.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include "../EasyOBJ.h"
#include "../MathTypeConverter.h"
#include "../MeshFix/MeshFix.h"

template <typename MeshType, int WNum>
class OrthoScanDeform
{
public:
    using Mesh = MeshType;
    using Kernel = typename Mesh::Traits::Kernel;

    void SetMesh(const MeshType* mesh, const std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>& crown_frames, bool upper)
    {
        _mesh = mesh;
        _crown_frames = crown_frames;
    }

    void SetCbctRegis( const CBCTRegis<double>* cbct_regis )
    {
        _cbct_regis = cbct_regis;
    }

    void SetCbctCentroids( const std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>& frames )
    {
        _cbct_centroids = frames;
    }

    MeshType Deform(const std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>& frames)
    {
        printf("preprocessing...");
        MeshType mesh = *_mesh;
        CGAL::set_halfedgeds_items_id(mesh);

        CGAL::Surface_mesh_deformation<MeshType> deformation(mesh);
        std::vector<typename MeshType::Vertex_handle> roi_vertices;
        for(auto hv : CGAL::vertices(mesh))
        {
            roi_vertices.push_back(static_cast<typename MeshType::Vertex_handle>(hv));
        }
        deformation.insert_roi_vertices(roi_vertices.begin(), roi_vertices.end());

        std::vector<typename MeshType::Vertex_handle> control_vertices;
        for(auto hv : CGAL::vertices(mesh))
            if(hv->_label != 0)
                control_vertices.push_back(static_cast<typename MeshType::Vertex_handle>(hv));
        deformation.insert_control_vertices(control_vertices.begin(), control_vertices.end());

        deformation.preprocess();

        printf("calculate target pos...");
        int count = 0;
        for(auto hv : CGAL::vertices(mesh))
        {
            if(hv->_label != 0)
            {
                Eigen::Vector3d p(hv->point().x(), hv->point().y(), hv->point().z());
                int label = hv->_label;
                if(label != 1)
                {
                    Eigen::Vector3d pos = _cbct_regis->CBCT_to_IOS(label) * frames.at(label) * (_cbct_regis->IOS_to_CBCT(label) * p - _cbct_centroids.at(label).translation());
                    deformation.set_target_position(hv, CGAL::ORIGIN + ToCGAL<double, Kernel>(pos));
                }
                else
                {
                    deformation.set_target_position(hv, hv->point());
                }
            }
        }

        printf("deform...");
        deformation.deform(10, 1e-4);
        
        printf("optimize...");
        std::unordered_set<typename MeshType::Facet_handle> faces_to_remove;
        for(auto hf : CGAL::faces(mesh))
        {
            int non_zero_label = 0;
            int l0 = hf->halfedge()->vertex()->_label;
            int l1 = hf->halfedge()->next()->vertex()->_label;
            int l2 = hf->halfedge()->prev()->vertex()->_label;
            if(l0 != 0)
                non_zero_label++;
            if(l1 != 0)
                non_zero_label++;
            if(l2 != 0)
                non_zero_label++;
            if(non_zero_label <= 1)
                continue;
            if((l0 != 0 && l1 != 0 && l0 != l1) || (l1 != 0 && l2 != 0 && l1 != l2) || (l0 != 0 && l2 != 0 && l0 != l2))
            {
                faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(hf));
                for(auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                {
                    faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(nei));
                }
            }
        }
        
        for(auto hf : faces_to_remove)
        {
            mesh.erase_facet(hf->halfedge());
        }

        auto [vertices, faces] = mesh.ToVerticesTriangles();
        auto labels = mesh.WriteLabels();
        std::vector<std::pair<std::vector<typename MeshType::Vertex_handle>, std::vector<typename MeshType::Facet_handle>>> patch;
        FixMesh(vertices, faces, mesh, true, 1000, true, false, 100, 100, true, 10, false, &patch);

        for(auto& pair : patch)
        {
            std::unordered_set<typename MeshType::Vertex_handle> vertex_to_smooth;
            for(auto hf : pair.second)
            {
                vertex_to_smooth.insert(hf->halfedge()->vertex());
                vertex_to_smooth.insert(hf->halfedge()->next()->vertex());
                vertex_to_smooth.insert(hf->halfedge()->prev()->vertex());
            }
            CGAL::Polygon_mesh_processing::fair(
                mesh, std::vector<typename MeshType::Vertex_handle>(vertex_to_smooth.begin(), vertex_to_smooth.end()));
            // std::unordered_set<typename MeshType::Facet_handle> facet_to_smooth;
            // for(auto hf : pair.second)
            // {
            //     facet_to_smooth.insert(hf);
            //     for(auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
            //     {
            //         facet_to_smooth.insert(nei);
            //     }
            // }
            // CGAL::Polygon_mesh_processing::angle_and_area_smoothing<MeshType, std::vector<typename MeshType::Facet_handle>>(
            //     std::vector<typename MeshType::Facet_handle>(facet_to_smooth.begin(), facet_to_smooth.end()), mesh,
            //     CGAL::parameters::do_project(false).number_of_iterations(10));
        }

        printf("finishing...");
        return std::move(mesh);
    }

protected:
    const MeshType* _mesh = nullptr;
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> _crown_frames;
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> _cbct_centroids;
    const CBCTRegis<double>* _cbct_regis;
};

#endif
