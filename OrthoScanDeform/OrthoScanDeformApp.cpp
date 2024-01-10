#include "OrthoScanDeform.h"
#include <CGAL/boost/graph/io.h>
#include <vector>
#include <unordered_map>
#include "../MeshFix/MeshFix.h"

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
            for(int r = 0; r < 4; r++)
                for(int c = 0; c < 4; c++)
                    mat(r, c) = frame_data[r][c];
            frames[label] = Eigen::Transform<double, 3, Eigen::Affine>(mat);
        }
    }
    return frames;
}

int main(int argc, char* argv[])
{
    using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string input_file = "";
    std::string label_file = "";
    std::string frame_file = "";
    std::string path_file = "";
    std::string cbct_regis_file = "";
    for (int i = 1; i < argc; i++)
    {
        if (std::strcmp(argv[i], "-i") == 0)
        {
            input_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-l") == 0)
        {
            label_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-f") == 0)
        {
            frame_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-p") == 0)
        {
            path_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-c") == 0)
        {
            cbct_regis_file = std::string(argv[i + 1]);
        }
        else if (std::strcmp(argv[i], "-h") == 0)
        {
            return 0;
        }
    }
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
    auto aabb = CGAL::bbox_3(mesh.points_begin(), mesh.points_end());
    if(upper)
    {
        for(auto hv : CGAL::vertices(mesh))
        {
            if(hv->point().z() > aabb.zmax() - aabb.z_span() * 0.1)
            {
                hv->_label = 1;
            }
        }
    }
    else
    {
        for(auto hv : CGAL::vertices(mesh))
        {
            if(hv->point().z() < aabb.zmin() + aabb.z_span() * 0.1)
            {
                hv->_label = 1;
            }
        }
    }
    for(auto hf : CGAL::faces(mesh))
    {
        if(hf->halfedge()->vertex()->_label == 1 || hf->halfedge()->next()->vertex()->_label == 1 || hf->halfedge()->prev()->vertex()->_label == 1)
        {
            hf->_label = 1;
        }
    }

    // Load crown frames
    auto crown_frames = LoadCrownFrameEigen(frame_file);
    auto paths = LoadPathsEigen(path_file);

    // Load cbct registration
    CBCTRegis<double> cbct_regis(cbct_regis_file);

    OrthoScanDeform<Polyhedron, 3> deformer;
    deformer.SetMesh(&mesh, crown_frames, upper);
    deformer.SetCbctRegis(&cbct_regis);

    printf("Load %zd steps.\n", paths.size());
    try
    {
        deformer.Deform(paths[0]);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}