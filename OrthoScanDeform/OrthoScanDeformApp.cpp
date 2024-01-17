#include "OrthoScanDeform.h"
#include <vector>
#include <unordered_map>
#include <argparse/argparse.hpp>
#include <CGAL/boost/graph/io.h>

template <typename Kernel>
std::vector<CrownFrames<Kernel>> LoadPaths( const std::string& path )
{
    nlohmann::json json = nlohmann::json::parse(std::ifstream(path));
    std::vector<CrownFrames<Kernel>> paths;
    int step = 0;
    while(json.find(std::to_string(step)) != json.end())
    {
        CrownFrames<Kernel> frames;
        auto step_data = json[std::to_string(step)];
        for(int label = 11; label < 50; label++)
        {
            if(step_data.find(std::to_string(label)) != step_data.end())
            {
                std::vector<double> mat = step_data[std::to_string(label)].get<std::vector<double>>();
                Frame<Kernel> frame( typename Kernel::Aff_transformation_3(
                    mat[0], mat[1], mat[2], mat[3],
                    mat[4], mat[5], mat[6], mat[7],
                    mat[8], mat[9], mat[10], mat[11],
                    1.0
                ));
                frames.Insert(label, frame);
            }
        }
        paths.push_back(frames);
        step++;
    }
    
    return paths;
}

std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> LoadPathsEigen( const std::string& path )
{
    nlohmann::json json = nlohmann::json::parse(std::ifstream(path));
    std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> paths;
    int step = 0;
    while(json.find(std::to_string(step)) != json.end())
    {
        std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> frames;
        auto step_data = json[std::to_string(step)];
        for(int label = 11; label < 50; label++)
        {
            if(step_data.find(std::to_string(label)) != step_data.end())
            {
                std::vector<double> data = step_data[std::to_string(label)].get<std::vector<double>>();
                Eigen::Matrix4d mat;
                for(int r = 0; r < 4; r++)
                    for(int c = 0; c < 4; c++)
                        mat(r, c) = data.at(r * 4 + c);
                frames[label] = Eigen::Transform<double, 3, Eigen::Affine>(mat);
            }
        }
        paths.push_back(frames);
        step++;
    }
    return paths;
}

std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> LoadCrownFrameEigen( const std::string& path )
{
    nlohmann::json json = nlohmann::json::parse(std::ifstream(path));
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> frames;
    for(int label = 11; label < 50; label++)
    {
        if(json.find(std::to_string(label)) != json.end())
        {
            std::vector<std::vector<double>> frame_data = json[std::to_string(label)].get<std::vector<std::vector<double>>>();
            Eigen::Matrix4d mat;
            for(int r = 0; r < 3; r++)
                for(int c = 0; c < 4; c++)
                    mat(r, c) = frame_data.at(c).at(r);
            mat(3, 0) = 0.0;
            mat(3, 1) = 0.0;
            mat(3, 2) = 0.0;
            mat(3, 3) = 1.0;
            frames[label] = Eigen::Transform<double, 3, Eigen::Affine>(mat);
        }
    }
    return frames;
}

std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> LoadCBCTTeethFrames( const std::string& path )
{
    Assimp::Importer importer;
    importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS,
     aiComponent_ANIMATIONS | aiComponent_COLORS | aiComponent_NORMALS | aiComponent_TEXCOORDS | aiComponent_TANGENTS_AND_BITANGENTS );
    const aiScene* scene = importer.ReadFile(path, aiProcess_JoinIdenticalVertices | aiProcess_RemoveComponent);
    if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode || scene->mNumMeshes < 1)
    {
        throw IOError("cannot read mesh: " + path);
    }
    printf("Load %d meshes\n", scene->mNumMeshes);
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> frames;
    for(unsigned i = 0; i < scene->mNumMeshes; i++)
    {
        const aiMesh* mesh = scene->mMeshes[i];
        
        int label = std::atoi(mesh->mName.C_Str());

        Eigen::Vector3d cent(0.0, 0.0, 0.0);
        double w_sum = 0.0;
        for(unsigned j = 0; j < mesh->mNumFaces; j++)
        {
            Eigen::Matrix3d m;
            m.col(0) = ToEigen(mesh->mVertices[mesh->mFaces[j].mIndices[0]]).cast<double>();
            m.col(1) = ToEigen(mesh->mVertices[mesh->mFaces[j].mIndices[1]]).cast<double>();
            m.col(2) = ToEigen(mesh->mVertices[mesh->mFaces[j].mIndices[2]]).cast<double>();
            double w = m.determinant();
            cent += (m.col(0) + m.col(1) + m.col(2)) / 4.0 * w;
            w_sum += w;
        }
        cent /= w_sum;

        frames.insert({label, Eigen::Transform<double, 3, Eigen::Affine>::Identity()});
        frames[label].pretranslate(cent);
    }
    printf("Loaded %zd cbct frames: ", frames.size());
    for(auto& [label, _] : frames)
    {
        printf("%d, ", label);
    }
    printf("\n");
    return frames;
}

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
    }
    auto [vertices, faces] = mesh.ToVerticesTriangles();
    FixMeshWithLabel(vertices, faces, mesh.WriteLabels(), mesh, true, 1000, true, false, 100, 100, true, 10);

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