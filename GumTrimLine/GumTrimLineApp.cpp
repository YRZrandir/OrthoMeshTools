#include "GumTrimLine.h"
#include <chrono>
#include <argparse/argparse.hpp>
                                                                                                                                                                                                                                                       
int main(int argc, char *argv[])
{
    argparse::ArgumentParser argparse("GumTrimLine");
    argparse.add_description("Extract gum trim line and output it as an .obj file.");
    argparse.add_argument("--input_file", "-i").required().help("specify the input mesh file.");
    argparse.add_argument("--label_file", "-l").required().help("specify the input label file.");
    argparse.add_argument("--frame_file", "-f").required().help("specify the crown frame file.");
    argparse.add_argument("--output_file", "-o").required().help("specify the output file.");
    argparse.add_argument("--smooth", "-s").default_value("10").help("a non-nagetive integer that specifies the iteration number of trim line smoothing");
    argparse.add_argument("--fix_factor", "-f").default_value(0.0).help("");
    try
    {
        argparse.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    
    try
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        GumTrimLine(argparse.get("input_file"), argparse.get("label_file"), argparse.get("frame_file"), argparse.get("output_file"), argparse.get<int>("smooth"), argparse.get<double>("fix_factor"));
        std::cout << "Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_time) << std::endl;
        std::cout << "===============================" << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}