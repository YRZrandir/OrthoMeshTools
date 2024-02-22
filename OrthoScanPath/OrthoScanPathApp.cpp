#include "OrthoScanPath.h"
#include <argparse/argparse.hpp>

int main(int argc, char* argv[])
{
    argparse::ArgumentParser args("OrthoScanPath");
    args.add_argument("-i").help("input path file").required();
    args.add_argument("-o").help("output path file").required();
    args.add_argument("-r").help("cbct registration").required();
    args.add_argument("-c").help("cbct teeth").required();
    args.parse_args(argc, argv);

    auto paths = OrthoScanPath(LoadPathsEigen(args.get("-i")), CBCTRegis<double>(args.get("-r")), LoadCBCTTeethFrames(args.get("-c")));
    
    nlohmann::json json;
    for(int step = 0; step < paths.size(); step++)
    {
        const auto& step_data = paths[step];
        nlohmann::json step_json;
        for(int label = 11; label < 50; label++)
        {
            if(step_data.count(label) == 0)
                continue;
            Eigen::Matrix4d mat = step_data.at(label).matrix();
            std::vector<double> data(16, 0.0);
            for(int r = 0; r < 4; r++)
                for(int c = 0; c < 4; c++)
                    data[r * 4 + c] = mat(r, c);
            step_json[std::to_string(label)] = data;
        }
        json[std::to_string(step)] = step_json;
    }
    std::ofstream ofs(args.get("-o"));
    ofs << json;
    return 0;
}