#include "MeshFix.h"

#ifndef FOUND_PYBIND11
#include <argparse/argparse.hpp>

extern bool gVerbose;

int main(int argc, char* argv[])
{
    argparse::ArgumentParser argparse("MeshFix");
    argparse.add_argument("--input", "-i").required().help("specify the input mesh file");
    argparse.add_argument("--output", "-o").required().help("specify the output mesh file");
    argparse.add_argument("--input_label", "-li").help("specify the input label file");
    argparse.add_argument("--output_label", "-lo").help("specify the output label file");
    argparse.add_argument("--small_component_threshold", "-sc").help("connected components whose size is smaller than the value are removed.").nargs(1).scan<'i', int>().default_value(0);
    argparse.add_argument("--fix_self_intersection").flag().help("detect and remove self-intersecting faces.");
    argparse.add_argument("--smallhole_edge_num", "-sh").help("holes whose edge number is smaller than the value are closed.").nargs(1).scan<'i', int>().default_value(0);
    argparse.add_argument("--smallhole_size", "-ss").help("holes whose edge bounding box smaller than the value are closed.").nargs(1).scan<'f', float>().default_value(0.0f);
    argparse.add_argument("--refine", "-r").help("refine the filled holes.").flag();
    argparse.add_argument("--max_retry", "-m").help("max retry number to fix the mesh.").scan<'i', int>().default_value(10);
    argparse.add_argument("--color", "-c").help("keep the colors").flag();
    try
    {
        argparse.parse_args(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }
    try
    {
        std::string path = argparse.get("-i");
        std::string output_path = argparse.get("-o");
        std::string input_label = argparse.present("-li").value_or("");
        std::string output_label = argparse.present("-lo").value_or("");
        int large_cc_threshold = argparse.get<int>("--small_component_threshold");
        bool keep_largest_connected_component = (large_cc_threshold != 0);
        bool fix_self_intersection = argparse.get<bool>("--fix_self_intersection");
        int smallhole_edge_num = argparse.get<int>("--smallhole_edge_num");
        float smallhole_size = argparse.get<float>("--smallhole_size");
        bool filter_small_holes = smallhole_edge_num <= 2 && smallhole_size <= 0.0;
        bool refine = argparse.get<bool>("--refine");
        int max_retry = argparse.get<int>("--max_retry");
        bool color = argparse.get<bool>("--color");
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
                smallhole_edge_num, 
                smallhole_size, 
                refine,
                max_retry
            );
        }
        else
        {
            if(color)
            {
                FixMeshFileWithColor( 
                    path,
                    output_path,
                    keep_largest_connected_component,
                    large_cc_threshold,
                    fix_self_intersection, 
                    filter_small_holes, 
                    smallhole_edge_num, 
                    smallhole_size, 
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
                    smallhole_edge_num, 
                    smallhole_size, 
                    refine,
                    max_retry
                );
            }
        }
    }
    catch( const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return -1;
    }
    return 0;
}
#endif