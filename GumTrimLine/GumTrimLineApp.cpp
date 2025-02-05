#include "GumTrimLine.h"
#include <chrono>
#include <argparse/argparse.hpp>
                                                                                                                                                                                                                                                       
int main(int argc, char *argv[])
{
    argparse::ArgumentParser argparse("GumTrimLine");
    argparse.add_description("Extract gum trim line and output it as an .obj file.");
    argparse.add_argument("--input_file", "-i").required().help("specify the input mesh file.");
    argparse.add_argument("--label_file", "-l").required().help("specify the input label file.");
    argparse.add_argument("--frame_file", "-fr").default_value("").help("specify the crown frame file.");
    argparse.add_argument("--output_file", "-o").required().help("specify the output file.");
    argparse.add_argument("--smooth", "-s").scan<'i', int>().default_value(10).help("a non-nagetive integer that specifies the iteration number of trim line smoothing");
    argparse.add_argument("--fix", "-f").flag().help("Turn on convexity fixing");
    argparse.add_argument("--debug_output", "-d").flag().help("Output some middle results.");
    argparse.add_argument("--hvalue").scan<'f', float>();
    argparse.add_argument("--mu").scan<'f', float>();
    try
    {
        argparse.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cerr << "Invalid arguments: " << e.what() << '\n';
    }
    
    try
    {
        GumTrimLine(argparse.get("-i"), argparse.get("-l"), "", argparse.get("-o"), argparse.get<int>("-s"), argparse.get<bool>("-f"), argparse.get<bool>("-d"));

    }
    catch(const std::exception& e)
    {
        std::cerr << "GumTrimLine Error: " << e.what() << std::endl;
        return -1;
    }
    return 0;
}