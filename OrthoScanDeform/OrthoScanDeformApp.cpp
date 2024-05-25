#include "OrthoScanDeform.h"
#include <vector>
#include <unordered_map>
#include <CGAL/boost/graph/io.h>
#include <argparse/argparse.hpp>

template <class Refs, typename Tag, typename Point>
class VertexDeform : public VertexWithLabelFlag<Refs, Tag, Point>
{
public:
    VertexDeform() = default;
    explicit VertexDeform(const Point &p) : VertexWithLabelFlag<Refs, Tag, Point>(p) {}

public:
    Point ori_pos;
};

class ItemsDeform : public ItemsWithLabelFlag
{
public:
    template <class Refs, class Traits>
    struct Vertex_wrapper
    {
        using Point = typename Traits::Point_3;
        using Vertex = VertexDeform<Refs, CGAL::Tag_true, Point>;
    };
};

int main(int argc, char* argv[])
{
    Eigen::Matrix3f m0;
    using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Polyhedron = TPolyhedronWithLabel<ItemsDeform, KernelEpick>;
    argparse::ArgumentParser argparse("OrthoScanDeform");
    argparse.add_argument("--input_file", "-i").required();
    argparse.add_argument("--label_file", "-l").required();
    argparse.add_argument("--frame_file", "-f").required();
    argparse.add_argument("--path_file", "-p").required();
    argparse.add_argument("--cbct_regis_file", "-c").required();
    argparse.add_argument("--cbct_teeth", "-t").required();
    try
    {
        argparse.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    std::string input_file = argparse.get("-i");
    std::string label_file = argparse.get("-l");
    std::string frame_file = argparse.get("-f");
    std::string path_file = argparse.get("-p");
    std::string cbct_regis_file = argparse.get("-c");
    std::string cbct_teeth = argparse.get("-t");
   
#ifdef _DEBUG
    std::filesystem::current_path(R"(D:\dev\Ortho\OrthoMeshTools\test\MeshDeform)");
    input_file = "oral_scan_L.ply";
    label_file = "oral_scan_L.json";
    frame_file = "crown_frame.json";
    path_file = "path.json";
    cbct_regis_file = "registration.json";
    cbct_teeth = "teeth_fdi.glb";
#endif
    auto start_time = std::chrono::high_resolution_clock::now();
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
            MeshFixConfig cfg;
            cfg.keep_largest_connected_component = false;
            cfg.fix_self_intersection = false;
            cfg.filter_small_holes = false;
            cfg.refine = false;
            cfg.fair = false;
            FixMeshWithLabel(vertices, faces, labels, mesh, cfg);
        }
        catch(const std::exception&)
        {
            throw IOError("Cannot read mesh file or mesh invalid: " + input_file);
        }
    }
    else
    {
        mesh.LoadLabels(label_file);
    }
    auto [vertices, faces] = mesh.ToVerticesTriangles();
    MeshFixConfig cfg;
    cfg.keep_largest_connected_component = true;
    cfg.large_cc_threshold = 1000;
    cfg.fix_self_intersection = true;
    cfg.filter_small_holes = false;
    cfg.refine = true;
    cfg.fair = true;
    FixMeshWithLabel(vertices, faces, mesh.WriteLabels(), mesh, cfg);

    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label % 10 == 8)
        {
            hv->_label = 0;
        }
    }
    mesh.UpdateFaceLabels();

    printf("Load mesh: V=%zd, F=%zd\n", mesh.size_of_vertices(), mesh.size_of_facets());
    // test upper or lower ortho mesh
    bool upper = true;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label >= 11)
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

    // A temp solution to determine bottom part (which won't deform).
    // auto aabb = CGAL::bbox_3(mesh.points_begin(), mesh.points_end());
    // if(upper)
    // {
    //     for(auto hv : CGAL::vertices(mesh))
    //     {
    //         if(hv->point().z() > aabb.zmax() - aabb.z_span() * 0.1)
    //         {
    //             hv->_label = 1;
    //         }
    //     }
    // }
    // else
    // {
    //     for(auto hv : CGAL::vertices(mesh))
    //     {
    //         if(hv->point().z() < aabb.zmin() + aabb.z_span() * 0.1)
    //         {
    //             hv->_label = 1;
    //         }
    //     }
    // }
    // for(auto hf : CGAL::faces(mesh))
    // {
    //     if(hf->halfedge()->vertex()->_label == 1 || hf->halfedge()->next()->vertex()->_label == 1 || hf->halfedge()->prev()->vertex()->_label == 1)
    //     {
    //         hf->_label = 1;
    //     }
    // }

    try
    {
        printf("Loading crown frames...");
        auto crown_frames = LoadCrownFrameEigen(frame_file);
        printf("Loading paths...");
        auto paths = LoadPathsEigen(path_file);
        printf("Loading cbct registration...");
        CBCTRegis<double> cbct_regis(cbct_regis_file);
        printf("Loading cbct teeth...");
        auto cbct_frames = LoadCBCTTeethFrames(cbct_teeth);

        OrthoScanDeform<Polyhedron, 3> deformer;
        deformer.SetMesh(&mesh, crown_frames, upper);
        deformer.SetCbctRegis(&cbct_regis);
        deformer.SetCbctCentroids(cbct_frames);
        
        printf("Load %zd steps.\n", paths.size());
        deformer.Deform(paths);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}