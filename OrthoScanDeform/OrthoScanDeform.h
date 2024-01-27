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

    void SetMesh(const MeshType *mesh, const std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> &crown_frames, bool upper)
    {
        _mesh = *mesh;
        _crown_frames = crown_frames;
    }

    void SetCbctRegis(const CBCTRegis<double> *cbct_regis)
    {
        _cbct_regis = cbct_regis;
    }

    void SetCbctCentroids(const std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> &frames)
    {
        _cbct_centroids = frames;
    }

    MeshType Deform(const std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> &frames)
    {
        printf("preprocessing...");
        MeshType mesh = _mesh;
        CGAL::set_halfedgeds_items_id(mesh);

        CGAL::Surface_mesh_deformation<MeshType> deformation(mesh);
        std::vector<typename MeshType::Vertex_handle> roi_vertices;
        for (auto hv : CGAL::vertices(mesh))
        {
            roi_vertices.push_back(static_cast<typename MeshType::Vertex_handle>(hv));
        }
        deformation.insert_roi_vertices(roi_vertices.begin(), roi_vertices.end());

        std::vector<typename MeshType::Vertex_handle> control_vertices;
        for (auto hv : CGAL::vertices(mesh))
            if (hv->_label != 0)
                control_vertices.push_back(static_cast<typename MeshType::Vertex_handle>(hv));
        deformation.insert_control_vertices(control_vertices.begin(), control_vertices.end());

        deformation.preprocess();

        printf("calculate target pos...");
        int count = 0;
        for (auto hv : control_vertices)
        {
            if (hv->_label != 0)
            {
                Eigen::Vector3d p = ToEigen(hv->point());
                int label = hv->_label;
                if (label != 1)
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
        for (auto hf : CGAL::faces(mesh))
        {
            int non_zero_label = 0;
            int l0 = hf->halfedge()->vertex()->_label;
            int l1 = hf->halfedge()->next()->vertex()->_label;
            int l2 = hf->halfedge()->prev()->vertex()->_label;
            if (l0 != 0)
                non_zero_label++;
            if (l1 != 0)
                non_zero_label++;
            if (l2 != 0)
                non_zero_label++;
            if (non_zero_label <= 1)
                continue;
            if ((l0 != 0 && l1 != 0 && l0 != l1) || (l1 != 0 && l2 != 0 && l1 != l2) || (l0 != 0 && l2 != 0 && l0 != l2))
            {
                faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(hf));
                for (auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                {
                    faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(nei));
                }
            }
        }

        for (auto hf : faces_to_remove)
        {
            mesh.erase_facet(hf->halfedge());
        }

        auto [vertices, faces] = mesh.ToVerticesTriangles();
        auto labels = mesh.WriteLabels();
        std::vector<std::pair<std::vector<typename MeshType::Vertex_handle>, std::vector<typename MeshType::Facet_handle>>> patch;
        FixMeshWithLabel(vertices, faces, labels, mesh, true, 1000, true, false, 100, 100, true, 10, &patch);
        for (auto &pair : patch)
        {
            std::unordered_set<typename MeshType::Vertex_handle> vertex_to_smooth;
            for (auto hf : pair.second)
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

    void Deform(const std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> &frames)
    {
        if (frames.empty())
        {
            return;
        }

        MeshType mesh = _mesh;
        for (auto hv : CGAL::vertices(mesh))
        {
            hv->ori_pos = hv->point();
        }
        for (int step = 0; step < frames.size(); step++)
        {
            printf("preprocessing...");
            CGAL::set_halfedgeds_items_id(mesh);

            CGAL::Surface_mesh_deformation<MeshType, CGAL::Default, CGAL::Default, CGAL::Deformation_algorithm_tag::SRE_ARAP> deformation(mesh);
            deformation.set_sre_arap_alpha(100.0);
            std::vector<typename MeshType::Vertex_handle> roi_vertices;
            for (auto hv : CGAL::vertices(mesh))
            {
                if(hv->_label != 1)
                    roi_vertices.push_back(static_cast<typename MeshType::Vertex_handle>(hv));
            }
            deformation.insert_roi_vertices(roi_vertices.begin(), roi_vertices.end());

            std::vector<typename MeshType::Vertex_handle> control_vertices;
            for (auto hv : CGAL::vertices(mesh))
                if (hv->_label >= 11 && hv->_label <= 49 || hv->_label == 1)
                    control_vertices.push_back(static_cast<typename MeshType::Vertex_handle>(hv));

            // For better deformation, we want to control the positions of some gum vertices that are close to the gumline.
            // This can be very time-consuming. maybe need a better solution.
            // std::unordered_map<Polyhedron::Vertex_handle, int> gumv_label_map;
            // {
            //     using LA = CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<typename CGAL::Eigen_sparse_matrix<double>::EigenType>>;
            //     using HeatMethod = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Polyhedron, CGAL::Heat_method_3::Intrinsic_Delaunay,
            //                                                                             typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::const_type, LA, Polyhedron::Traits::Kernel>;    
            //     std::unordered_map<int, std::unordered_map<typename Polyhedron::Vertex_handle, double>> distances;
            //     std::unordered_map<int, std::vector<Polyhedron::Vertex_handle>> labelwise_borders;
            //     HeatMethod heat_method(mesh);
            //     // for(auto hh : CGAL::halfedges(mesh))
            //     // {
            //     //     if(hh->facet() && hh->opposite() && hh->opposite()->facet() && hh->facet()->_label != 0 && hh->opposite()->facet()->_label == 0)
            //     // }
            //     for(auto hv : CGAL::vertices(mesh))
            //     {
            //         if(hv->_label != 0 && hv->_label != 1)
            //             labelwise_borders[hv->_label].push_back(hv);
            //     }

            //     for(auto& [label, sources] : labelwise_borders)
            //     {
            //         heat_method.clear_sources();
            //         heat_method.add_sources(sources);
            //         heat_method.estimate_geodesic_distances(distances[label]);
            //     }

            //     for(auto hv : CGAL::vertices(mesh))
            //     {
            //         if(hv->_label != 0)
            //             continue;
            //         if(gumv_label_map.count(hv) != 0)
            //             continue;
                    
            //     }
            // }
            deformation.insert_control_vertices(control_vertices.begin(), control_vertices.end());

            if(!deformation.preprocess())
            {
                throw AlgError("Failed to deform step " + std::to_string(step));
            }

            printf("calculate target pos...");
            if (step == 0)
            {
                for (auto hv : CGAL::vertices(mesh))
                {
                    if (hv->_label != 0)
                    {
                        Eigen::Vector3d p = ToEigen(hv->ori_pos);
                        int label = hv->_label;
                        if (label >= 11)
                        {
                            Eigen::Vector3d pos = _cbct_regis->CBCT_to_IOS(label) * frames[step].at(label) * (_cbct_regis->IOS_to_CBCT(label) * p - _cbct_centroids.at(label).translation());
                            deformation.set_target_position(hv, CGAL::ORIGIN + ToCGAL<double, Kernel>(pos));
                        }
                        else if(label == 1)
                        {
                            deformation.set_target_position(hv, hv->point());
                        }
                    }
                }
            }
            else
            {
                for (auto hv : CGAL::vertices(mesh))
                {
                    if (hv->_label != 0)
                    {
                        int label = hv->_label;
                        if (label >= 11)
                        {
                            Eigen::Vector3d p = ToEigen(hv->point());
                            Eigen::Vector3d p_cbct = _cbct_regis->CBCT_to_IOS(label).inverse() * p;
                            Eigen::Vector3d pos = _cbct_regis->CBCT_to_IOS(label) * frames[step].at(label) * frames[step - 1].at(label).inverse() * p_cbct;
                            deformation.set_target_position(hv, CGAL::ORIGIN + ToCGAL<double, Kernel>(pos));
                        }
                        else if(label == 1)
                        {
                            deformation.set_target_position(hv, hv->point());
                        }
                    }
                }
            }
            printf("deform...");
            deformation.deform(10, 1e-4);

            printf("optimize...");
            std::unordered_set<typename MeshType::Facet_handle> faces_to_remove;
            for (auto hf : CGAL::faces(mesh))
            {
                int non_zero_label = 0;
                int l0 = hf->halfedge()->vertex()->_label;
                int l1 = hf->halfedge()->next()->vertex()->_label;
                int l2 = hf->halfedge()->prev()->vertex()->_label;
                if(l0 == 3 || l1 == 3 || l2 == 3)
                {
                    faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(hf));
                    for (auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                    {
                        faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(nei));
                    }
                    continue;
                }
                if (l0 != 0)
                    non_zero_label++;
                if (l1 != 0)
                    non_zero_label++;
                if (l2 != 0)
                    non_zero_label++;
                if (non_zero_label <= 1)
                    continue;
                if ((l0 != 0 && l1 != 0 && l0 != l1) || (l1 != 0 && l2 != 0 && l1 != l2) || (l0 != 0 && l2 != 0 && l0 != l2))
                {
                    faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(hf));
                    for (auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                    {
                        faces_to_remove.insert(static_cast<typename MeshType::Facet_handle>(nei));
                    }
                }
            }

            for (auto hf : faces_to_remove)
            {
                mesh.erase_facet(hf->halfedge());
            }

            //mesh.WriteOBJ("mid.obj");
            auto [vertices, faces] = mesh.ToVerticesTriangles();
            std::vector<std::pair<std::vector<typename MeshType::Vertex_handle>, std::vector<typename MeshType::Facet_handle>>> patch;
            FixMeshWithLabel(vertices, faces, mesh.WriteLabels(), mesh, true, 1000, true, false, 100, 100, true, 10, &patch);
            for (auto &pair : patch)
            {
                std::unordered_set<typename MeshType::Vertex_handle> vertex_to_smooth;
                for (auto hf : pair.second)
                {
                    vertex_to_smooth.insert(hf->halfedge()->vertex());
                    vertex_to_smooth.insert(hf->halfedge()->next()->vertex());
                    vertex_to_smooth.insert(hf->halfedge()->prev()->vertex());
                }
                CGAL::Polygon_mesh_processing::fair(
                    mesh, std::vector<typename MeshType::Vertex_handle>(vertex_to_smooth.begin(), vertex_to_smooth.end()));
            }
            mesh.UpdateFaceLabels();
            mesh.WriteOBJ("deformed_step" + std::to_string(step) + ".obj");
            mesh.WriteLabels("deformed_step" + std::to_string(step) + ".json");
            printf("finishing...");
        }
    }

protected:
    MeshType _mesh;
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> _crown_frames;
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> _cbct_centroids;
    const CBCTRegis<double> *_cbct_regis;
};

#endif
