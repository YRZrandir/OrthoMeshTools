#include "Ortho.h"
#include "MathTypeConverter.h"
#include <assimp/scene.h>
#include <assimp/mesh.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>

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
