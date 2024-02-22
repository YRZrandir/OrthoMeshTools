#include "Polyhedron.h"
#include <argparse/argparse.hpp>
#include "OrthoScanPath/OrthoScanPath.h"
#include "MathTypeConverter.h"

using Kernel = CGAL::Simple_cartesian<double>;
using Polyhedron = TPolyhedron<CGAL::Polyhedron_items_with_id_3, Kernel>;

void MoveTeeth(
    std::vector<std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>>> input_path,
    const CBCTRegis<double>& cbct_regis,
    std::unordered_map<int, Eigen::Transform<double, 3, Eigen::Affine>> cbct_frames,
    int step
)
{
    auto paths = OrthoScanPath(input_path, cbct_regis, cbct_frames);
    if(paths.size() <= step)
        return;
    std::unordered_map<int, std::vector<typename Kernel::Point_3>> teeth_vertices;
    std::unordered_map<int, std::vector<TTriangle<size_t>>> teeth_faces;
    for(int label = 11; label < 50; label++)
    {
        std::ifstream ifs(std::to_string(label) + ".obj");
        if(ifs.good())
        {
            std::vector<typename Kernel::Point_3> vertices;
            std::vector<TTriangle<size_t>> faces;
            LoadVFAssimp<Kernel, size_t>(std::to_string(label) + ".obj", vertices, faces);
            teeth_vertices[label] = std::move(vertices);
            teeth_faces[label] = std::move(faces);
        }
    }

    auto& transforms = paths[step];
    for(auto& [label, vertices] : teeth_vertices)
    {
        if(transforms.count(label) == 0)
        {
            continue;
        }
        for(auto& v : vertices)
        {
            auto p = ToEigen(v);
            p = transforms[label] * p;
            v = Kernel::Point_3(p.x(), p.y(), p.z());
        }
        WriteVFObj<Kernel, size_t>("./" + std::to_string(step) + "/" + std::to_string(label) + ".obj", vertices, teeth_faces[label]);
    }
}

int main(int argc, char* argv[])
{
    argparse::ArgumentParser args("MoveTeeth");
    args.add_argument("-p").help("path file").required();
    args.add_argument("-r").help("cbct registration").required();
    args.add_argument("-c").help("cbct teeth").required();
    args.add_argument("-s").help("step").required().nargs(1).scan<'i', int>();
    try
    {
        args.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }
    MoveTeeth(LoadPathsEigen(args.get("-p")), CBCTRegis<double>(args.get("-r")), LoadCBCTTeethFrames(args.get("-c")), args.get<int>("-s"));
    return 0;
}