#include "MeshFix.h"
#include <argparse/argparse.hpp>

extern bool gVerbose;

int main(int argc, char* argv[])
{
    auto print_help_msg = []()
    {
        std::cout << "usage:\n"
        "\t-i filename \tPath to input mesh.\n"
        "\t-o filename \tFile name of output mesh.\n"
        "\t-k threshold\tDelete connected components smaller than threshold (default=off)\n"
        "\t-s \tFix self intersection\n"
        "\t-f max_hole_edges max_hole_diam\t Do not fill holes that satisfiy (edge number > max_hole_edges) OR (AABB size > max_hole_diam)\n"
        "\t-r refine after filling holes.\n"
        "\t-m max_retry The program will repeatedly try to fix the mesh, this is the max retry time. (default=10)"
        "\t-v \tPrint debug messages" << std::endl;
    };
    if(argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
    {
        print_help_msg();
        return -1;
    }
    argparse::ArgumentParser argparse("MeshFix");
    argparse.add_argument("--input", "-i").required().help("specify the input mesh file");
    argparse.add_argument("--output", "-o").required().help("specify the output mesh file");
    argparse.add_argument("--input_label", "-li").help("specify the input label file");
    argparse.add_argument("--output_label", "-lo").help("specify the output label file");
    argparse.add_argument("--remove_small_components", "-k").help("").nargs(1).scan<'g', double>().default_value(0);
    argparse.add_argument("--fix_self_intersection").flag();
    argparse.add_argument("--filter_small_holes", "-f").help("").nargs(2).default_value(std::pair<int, double>(0, 0.0));
    argparse.add_argument("--refine", "-r").help("").flag();
    argparse.add_argument("--max_retry", "-m").help("").scan<'i', int>().default_value(10);
    
    std::string path;
    std::string output_path;
    std::string input_label;
    std::string output_label;
    bool keep_largest_connected_component = false;
    int large_cc_threshold = 100;
    bool fix_self_intersection = false;
    bool filter_small_holes = false;
    int max_hole_edges = std::numeric_limits<int>::max();
    float max_hole_diam = std::numeric_limits<float>::max();
    bool refine = false;
    int max_retry = 10;

    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-i") == 0)
        {
            path = std::string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-o") == 0)
        {
            output_path = std::string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            gVerbose = true;
        }
        else if (strcmp(argv[i], "-k") == 0)
        {
            keep_largest_connected_component = true;
            if(i < argc - 1 && std::atoi(argv[i+1]) != 0)
            {
                large_cc_threshold = std::atoi(argv[i + 1]);
            }
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            fix_self_intersection = true;
        }
        else if (strcmp(argv[i], "-f") == 0)
        {
            filter_small_holes = true;
            max_hole_edges = std::atoi(argv[i+1]);
            max_hole_diam = static_cast<float>(std::atof(argv[i+2]));
        }
        else if (strcmp(argv[i], "-r") == 0)
        {
            refine = true;
        }
        else if (strcmp(argv[i], "-m") == 0)
        {
            max_retry = std::atoi(argv[i+1]);
        }
        else if (strcmp(argv[i], "-li") == 0)
        {
            input_label = std::string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-lo") == 0)
        {
            output_label = std::string(argv[i + 1]);
        }
    }
    if(path.empty() || output_path.empty())
    {
        print_help_msg();
        return -1;
    }
    try
    {
        if(!input_label.empty() && !output_label.empty())
        {
            FixMeshFileWithLabel(
                path,
                output_path,
                input_label,
                output_label,
                keep_largest_connected_component,
                large_cc_threshold,
                fix_self_intersection, 
                filter_small_holes, 
                max_hole_edges, 
                max_hole_diam, 
                refine,
                max_retry
            );
        }
        else
        {
            FixMeshFile( 
                path,
                output_path,
                keep_largest_connected_component,
                large_cc_threshold,
                fix_self_intersection, 
                filter_small_holes, 
                max_hole_edges, 
                max_hole_diam, 
                refine,
                max_retry
            );
        }
    }
    catch( const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }
    return 0;
}