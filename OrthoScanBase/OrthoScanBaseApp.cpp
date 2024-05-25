#include "OrthoScanBase.h"
#include <argparse/argparse.hpp>
#include "../Polyhedron.h"
#include "../MeshFix/MeshFix.h"
#include "../EasyOBJ.h"
#include "../SegClean/SegClean.h"

using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = TPolyhedronWithLabel<ItemsWithLabelFlag, KernelEpick>;

void LabelProcessing(Polyhedron& mesh);
void GenerateBase2(Polyhedron &mesh);
void Optimize(Polyhedron &mesh);
void GenerateGumImpl(std::string output_gum, std::string crown_frame, bool upper, Polyhedron& mesh);

int main(int argc, char *argv[])
{
    argparse::ArgumentParser parser;
    parser.add_argument("--input_file", "-i").required().help("specify the input mesh.");
    parser.add_argument("--input_label", "-l").required().help("specify the input labels.");
    parser.add_argument("--output_file", "-o").help("specify the output file.");
    parser.add_argument("--output_label", "-ol").help("specify the output labels.");
    parser.add_argument("--output_gum", "-g").help("specify the output gum mesh file.");
    parser.add_argument("--crown_frame", "-c").help("specify the crown frame file.");
    try
    {
        parser.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }

    Polyhedron mesh;
    std::string input_file = parser.get("-i");
    std::string label_file = parser.get("-l");
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
            cfg.keep_largest_connected_component = true;
            cfg.large_cc_threshold = 9999;
            cfg.fix_self_intersection = false;
            cfg.filter_small_holes = false;
            cfg.max_hole_diam = 0.f;
            cfg.max_hole_edges = 0;
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
        auto [v, f] = mesh.ToVerticesTriangles();
        MeshFixConfig cfg;
        cfg.keep_largest_connected_component = true;
        cfg.large_cc_threshold = 9999;
        cfg.fix_self_intersection = false;
        cfg.filter_small_holes = false;
        cfg.max_hole_diam = 0.f;
        cfg.max_hole_edges = 0;
        FixMeshWithLabel(v, f, mesh.WriteLabels(), mesh, cfg);
        //mesh.UpdateFaceLabels2();
    }
    LabelProcessing(mesh);
    SegClean(mesh);

    bool upper = true;
    for(auto hv : CGAL::vertices(mesh))
    {
        if(hv->_label >= 11 && hv->_label <= 29)
            break;
        else if(hv->_label >= 31 && hv->_label <= 49)
        {
            upper = false;
            break;
        }
    }
    try
    {
        std::cout << "Optimizing...";
        Optimize(mesh);
        auto [v, f] = mesh.ToVerticesTriangles();
        MeshFixConfig cfg;
        cfg.keep_largest_connected_component = true;
        cfg.large_cc_threshold = 9999;
        cfg.fix_self_intersection = false;
        cfg.filter_small_holes = false;
        cfg.max_hole_diam = 0.f;
        cfg.max_hole_edges = 0;
        FixMeshWithLabel(v, f, mesh.WriteLabels(), mesh, cfg);
        std::cout << "Done." << std::endl;
#ifdef DEBUG_ORTHOSCANBASE
        mesh.WriteOBJ("optimized.obj");
#endif
        std::cout << "Generating...";
        GenerateBase2(mesh);
        std::cout << "Done." << std::endl;
        if(parser.present("-o").has_value()) {
            mesh.WriteOBJ(parser.get("-o"));
        }
        if(parser.present("-ol").has_value()) {
            mesh.WriteLabels(parser.get("-ol"));
        }

        if(parser.present("-g").has_value() && parser.present("-c").has_value())
        {
            GenerateGumImpl(parser.present("-g").value(), parser.present("-c").value(), upper, mesh);
        }
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    return 0;
}