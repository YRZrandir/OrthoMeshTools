#include "OrthoScanDeform.h"
#include <CGAL/boost/graph/io.h>
#include "../MeshFix/MeshFix.h"

int main(int argc, char* argv[])
{
    using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string input_file = "";
    std::string label_file = "";
    std::string frame_file = "";
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
        //mesh.UpdateFaceLabels2();
    }
    CrownFrames<KernelEpick> crown_frames(frame_file);
    OrthoScanDeform<Polyhedron, 3> deformer;
    deformer.SetMesh(&mesh, &crown_frames);
    return 0;
}